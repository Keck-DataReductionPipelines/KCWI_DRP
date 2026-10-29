"""Draw a recorded ``PlotFigure``/``PlotGrid`` (see ``kcwi_plotting``) into a
Matplotlib ``Figure``.  Private to ``kcwi_plotting``; call
``kcwi_plotting.render``.

Bokeh conventions kept (bokeh-plan 1.1); otherwise plots use Matplotlib's
own look, as agreed with the maintainer (same information, simpler code):

* sizes are in screen pixels at 100 dpi; one pixel is 0.72 point, for marker
  ``size`` and ``line_width``;
* default glyph colour ``#1f77b4``, no colour cycling; ``color`` sets both the
  outline and the fill of a marker, ``line_color``/``fill_color`` override;
* automatic ranges use Matplotlib's autoscaling with a 5 % margin at each
  end, which is Bokeh's default padding; ``range_padding = 0`` removes it;
  ``start``/``end`` override either end;
* glyphs are drawn in call order; one legend, entries merged by label, on the
  top-most axes; a left-aligned bold title;
* ``Spectral11`` is Bokeh's eleven discrete colours, mapped linearly from the
  array minimum to its maximum; images have row 0 at the bottom.

No ``pyplot``: every call builds or reuses its own ``Figure``, so renders on
several threads at once do not share Matplotlib state.
"""
import numpy as np

from matplotlib.backends.backend_agg import FigureCanvasAgg, RendererAgg
from matplotlib.colors import ListedColormap, to_rgba
from matplotlib.collections import PolyCollection
from matplotlib.figure import Figure
from matplotlib.path import Path

from kcwidrp.core import kcwi_plotting as kp

DPI = 100
PT_PER_PX = 72.0 / DPI          # 0.72
DEFAULT_COLOR = '#1f77b4'
MARGIN = 0.05                   # per end; Bokeh pads 10 % of the span in all
MIN_TITLE_PT = 7.0

SPECTRAL11 = ('#5e4fa2', '#3288bd', '#66c2a5', '#abdda4', '#e6f598',
              '#ffffbf', '#fee08b', '#fdae61', '#f46d43', '#d53e4f',
              '#9e0142')
# Built once and never modified afterwards, so it may be shared by threads.
_CMAPS = {'Spectral11': ListedColormap(list(SPECTRAL11), name='Spectral11')
          .with_extremes(bad='gray')}   # Bokeh's default nan_color

_DASHES = {'solid': '-', 'dashed': '--', 'dotted': ':', 'dashdot': '-.'}

# Bokeh's diamond: height = size, width = size / 1.5.  Matplotlib scales a
# custom path so its largest extent equals markersize.
_DIAMOND = Path([(0, 1), (1 / 1.5, 0), (0, -1), (-1 / 1.5, 0), (0, 1)],
                closed=True)
_MARKERS = {'circle': 'o', 'diamond': _DIAMOND, 'x': 'x'}


def _pt(px):
    return px * PT_PER_PX


def _inches(px, dpi):
    """px / dpi, nudged up so that Agg's truncation gives exactly px."""
    v = px / dpi
    while v * dpi < px:
        v = np.nextafter(v, np.inf)
    return v


# ---------------------------------------------------------------------------

def render(obj, fig=None):
    w, h = kp.pixel_size(obj)
    if fig is None:
        fig = Figure(dpi=DPI)
        FigureCanvasAgg(fig)
    else:
        fig.clear()
    dpi = getattr(fig, '_original_dpi', fig.dpi)
    fig.set_size_inches(_inches(w, dpi), _inches(h, dpi), forward=False)
    fig.set_layout_engine('constrained')
    fig.set_facecolor('white')

    if isinstance(obj, kp.PlotGrid):
        figs = obj._children
        nrows, ncols = obj._shape()
        cell_w = obj._cell_size()[0]
    else:
        figs = [obj]
        nrows = ncols = 1
        cell_w = w
    measure = RendererAgg(10, 10, fig.dpi)
    for i, pf in enumerate(figs):
        ax = fig.add_subplot(nrows, ncols, i + 1)
        _draw(ax, pf, cell_w, fig.dpi / dpi, measure)
    return fig


def _set_limits(ax, axis, rng, log):
    """Apply a recorded range to one axis of autoscaled axes.  An end that
    cannot be drawn (NaN or infinite, or not positive on a log axis) is left
    to autoscaling."""
    def usable(v):
        return v is not None and np.isfinite(v) and (v > 0 or not log)
    start = rng.start if usable(rng.start) else None
    end = rng.end if usable(rng.end) else None
    if start is not None or end is not None:
        set_lim = ax.set_xlim if axis == 'x' else ax.set_ylim
        set_lim(start, end)         # None keeps the autoscaled end


def _margin(rng):
    # A Range1d fixes both ends, so its margin never shows
    if isinstance(rng, kp.Range1d):
        return MARGIN
    return 0.0 if rng.range_padding == 0 else MARGIN


def _draw(ax, pf, cell_w, scale, measure):
    if pf._y_log:
        ax.set_yscale('log', nonpositive='mask')
    ax.grid(True)
    ax.set_axisbelow(True)
    if pf._xgrid.grid_line_color is None:
        ax.xaxis.grid(False)
    ax.set_xmargin(_margin(pf._x_range))
    ax.set_ymargin(_margin(pf._y_range))

    axes = {'default': ax}
    twin = None
    extra = pf._extra_x_ranges
    if len(extra) > 1:
        kp._unsupported("more than one extra_x_ranges entry")
    if extra:
        name, xrng2 = next(iter(extra.items()))
        twin = ax.twiny()
        axes[name] = twin
        axis = None
        for lay in pf._layouts:
            if lay.x_range_name == name and axis is None:
                axis = lay
            else:
                kp._unsupported("add_layout(%r) without a matching single "
                                "extra x range" % (lay,))
        if axis is None:
            twin.xaxis.set_visible(False)
        else:
            twin.set_xlabel(axis.axis_label, parse_math=False)
    else:
        for lay in pf._layouts:
            kp._unsupported("add_layout(%r) with no extra_x_ranges" % (lay,))

    legend = {}                         # label -> [artists], first-use order
    for i, g in enumerate(pf._glyphs):
        target = axes.get(g.get('x_range_name', 'default'))
        if target is None:
            kp._unsupported("%s(x_range_name=%r): no such range"
                            % (g['kind'], g.get('x_range_name')))
            continue
        try:
            artist = _GLYPHS[g['kind']](target, g)
        except Exception as e:
            kp._failed("drawing %s()" % g['kind'], e)
            continue
        artist.set_zorder(2 + i * 1e-4)     # Bokeh draws in call order
        label = g.get('legend_label')
        if label:
            legend.setdefault(label, []).append(artist)

    # Autoscale first, so a range with only one end set keeps the other
    # end autoscaled.
    for a in axes.values():
        a.autoscale_view()
    _set_limits(ax, 'x', pf._x_range, False)
    _set_limits(ax, 'y', pf._y_range, pf._y_log)
    if twin is not None:
        _set_limits(twin, 'x', xrng2, False)

    ax.set_xlabel(pf._x_axis_label, parse_math=False)
    ax.set_ylabel(pf._y_axis_label, parse_math=False)
    title = ax.set_title(pf._title.text, loc='left', fontweight='bold',
                         parse_math=False)
    _fit_title(title, cell_w, scale, measure)

    if legend:
        top = twin if twin is not None else ax
        handles = [a[0] if len(a) == 1 else tuple(a) for a in legend.values()]
        leg = top.legend(handles, list(legend.keys()),
                         loc=kp.LEGEND_LOCATIONS[pf._legend.location])
        for t in leg.get_texts():
            t.set_parse_math(False)


def _fit_title(title, cell_w, scale, measure):
    """Make a title that is wider than its cell (less room for the y-axis
    labels) fit, so constrained layout never squeezes the axes: shrink it
    down to MIN_TITLE_PT, then let Matplotlib wrap it at draw time (the
    text itself, ``get_title``, is unchanged)."""
    if not title.get_text():
        return
    avail = (cell_w - 80) * scale
    width = title.get_window_extent(renderer=measure).width
    if not width > avail > 0:
        return
    title.set_fontsize(max(MIN_TITLE_PT, title.get_fontsize() * avail / width))
    if title.get_window_extent(renderer=measure).width > avail:
        title.set_wrap(True)


# ---------------------------------------------------------------------------
# Glyphs: each draws one recorded glyph and returns its artist
# ---------------------------------------------------------------------------

def _first(g, *names, default=None):
    for n in names:
        if n in g:
            return g[n]
    return default


def _rgba(color, alpha):
    if color is None:
        return (0, 0, 0, 0)
    return to_rgba(color, alpha)


def _line(ax, g):
    lw = _pt(g.get('line_width', 1.0))
    color = _rgba(_first(g, 'line_color', 'color', default=DEFAULT_COLOR),
                  _first(g, 'line_alpha', 'alpha', default=1.0))
    line, = ax.plot(g['x'], g['y'], color=color, linewidth=lw,
                    linestyle=_DASHES[g.get('line_dash', 'solid')],
                    marker='None')
    return line


def _marker(ax, g):
    kind = g['kind']
    name = g.get('marker', 'circle') if kind == 'scatter' else kind
    edge_alpha = _first(g, 'line_alpha', 'alpha', default=1.0)
    edge = _rgba(_first(g, 'line_color', 'color', default=DEFAULT_COLOR),
                 edge_alpha)
    face = _rgba(_first(g, 'fill_color', 'color', default=DEFAULT_COLOR),
                 g.get('alpha', 1.0))
    mew = _pt(g.get('line_width', 1.0))
    if edge[3] == 0:
        edge, mew = 'none', 0.0     # line_alpha=0: filled, no outline
    line, = ax.plot(g['x'], g['y'], linestyle='None', marker=_MARKERS[name],
                    markersize=_pt(g.get('size', 4.0)),
                    markerfacecolor=face, markeredgecolor=edge,
                    markeredgewidth=mew)
    return line


def _quad(ax, g):
    left, right, bottom, top = g['left'], g['right'], g['bottom'], g['top']
    ok = np.isfinite(left) & np.isfinite(right) & np.isfinite(bottom) & \
        np.isfinite(top)
    verts = [((l, b), (r, b), (r, t), (l, t)) for l, r, b, t in
             zip(left[ok], right[ok], bottom[ok], top[ok])]
    face = _rgba(_first(g, 'fill_color', 'color', default=DEFAULT_COLOR),
                 g.get('alpha', 1.0))
    edge = _rgba(_first(g, 'line_color', 'color', default=DEFAULT_COLOR),
                 _first(g, 'line_alpha', 'alpha', default=1.0))
    pc = PolyCollection(verts, facecolors=[face], edgecolors=[edge],
                        linewidths=_pt(g.get('line_width', 1.0)))
    ax.add_collection(pc)
    return pc


def _image(ax, g):
    a = g['image']
    finite = a[np.isfinite(a)]
    if finite.size:
        vmin, vmax = float(finite.min()), float(finite.max())
    else:
        vmin, vmax = 0.0, 1.0
    x0, y0 = g['x'], g['y']
    x1, y1 = x0 + g['dw'], y0 + g['dh']
    im = ax.imshow(np.ma.masked_invalid(a), origin='lower',
                   extent=(x0, x1, y0, y1), cmap=_CMAPS[g['palette']],
                   vmin=vmin, vmax=vmax, interpolation='nearest',
                   aspect='auto', resample=False)
    return im


_GLYPHS = {'line': _line, 'circle': _marker, 'diamond': _marker,
           'scatter': _marker, 'quad': _quad, 'image': _image}
