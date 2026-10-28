"""Draw a recorded ``PlotFigure``/``PlotGrid`` (see ``kcwi_plotting``) into a
Matplotlib ``Figure``.  Private to ``kcwi_plotting``; call
``kcwi_plotting.render``.

Bokeh conventions reproduced here (bokeh-plan 1.1):

* sizes are in screen pixels at 100 dpi; one pixel is 0.72 point, for marker
  ``size``, ``line_width`` and font sizes;
* default glyph colour ``#1f77b4``, no colour cycling; ``color`` sets both the
  outline and the fill of a marker, ``line_color``/``fill_color`` override;
* automatic ranges are padded as ``DataRange1d`` pads them (``range_padding``
  is a fraction of the data span, half added at each end; in decades on a
  log axis); ``start``/``end`` override either end;
* glyphs are drawn in call order; one legend, entries merged by label, on the
  top-most axes;
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
import matplotlib as mpl

from kcwidrp.core import kcwi_plotting as kp

DPI = 100
PT_PER_PX = 72.0 / DPI          # 0.72
DEFAULT_COLOR = '#1f77b4'
GRID_COLOR = '#e5e5e5'
OUTLINE_COLOR = '#e5e5e5'
TITLE_PX = 13                   # Bokeh title: 13px bold
AXIS_LABEL_PX = 13              # Bokeh axis label: 13px italic
TICK_LABEL_PX = 11              # Bokeh major tick label: 11px
LEGEND_PX = 13                  # Bokeh legend label: 13px
MIN_TITLE_PT = 7.0

SPECTRAL11 = ('#5e4fa2', '#3288bd', '#66c2a5', '#abdda4', '#e6f598',
              '#ffffbf', '#fee08b', '#fdae61', '#f46d43', '#d53e4f',
              '#9e0142')
# Built once and never modified afterwards, so it may be shared by threads.
_CMAPS = {'Spectral11': ListedColormap(list(SPECTRAL11), name='Spectral11')
          .with_extremes(bad='gray')}   # Bokeh's default nan_color

# Bokeh dash patterns, in pixels (on, off, ...)
_DASHES = {'dashed': (6, 6), 'dotted': (2, 4), 'dashdot': (6, 4, 2, 4)}

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


class _Bounds:
    def __init__(self, log=False):
        self.lo = None
        self.hi = None
        self.log = log

    def add(self, *arrays):
        for a in arrays:
            a = np.asarray(a, dtype=float)
            ok = np.isfinite(a)
            if self.log:
                ok &= a > 0
            if not ok.any():
                continue
            lo, hi = float(a[ok].min()), float(a[ok].max())
            self.lo = lo if self.lo is None else min(self.lo, lo)
            self.hi = hi if self.hi is None else max(self.hi, hi)


def _auto_range(bounds, padding, log):
    """Bokeh DataRange1d: pad the data span by ``padding`` (total), or use
    the default span 2 (two decades on a log axis) for a single value."""
    if bounds.lo is None:
        return None
    lo, hi = bounds.lo, bounds.hi
    if log:
        lo, hi = np.log10(lo), np.log10(hi)
    if lo == hi:
        lo, hi = lo - 1.0, hi + 1.0
    else:
        c = 0.5 * (lo + hi)
        half = 0.5 * (hi - lo) * (1.0 + padding)
        lo, hi = c - half, c + half
    if log:
        lo, hi = 10.0 ** lo, 10.0 ** hi
    return lo, hi


def _limits(rng, bounds, log):
    """(start, end) for an axis.  An end that cannot be drawn (NaN or
    infinite, or not positive on a log axis) falls back to the data."""
    start, end = rng.start, rng.end
    # A Range1d has no padding of its own; 0.1 is used only for a fallback.
    padding = 0.1 if isinstance(rng, kp.Range1d) else rng.range_padding
    auto = _auto_range(bounds, padding, log)

    def usable(v):
        return v is not None and np.isfinite(v) and (v > 0 or not log)
    start = start if usable(start) else None
    end = end if usable(end) else None
    if start is None and auto is not None:
        start = auto[0]
    if end is None and auto is not None:
        end = auto[1]
    if start is None and end is None:
        return (1.0, 10.0) if log else (0.0, 1.0)
    if start is None:
        start = end / 10.0 if log else end - 1.0
    if end is None:
        end = start * 10.0 if log else start + 1.0
    return start, end


def _style_axes(ax, pf):
    # Limits are always set explicitly (Bokeh's padding), so Matplotlib's
    # autoscaling is never wanted; it would also warn on an empty log axis.
    ax.set_autoscale_on(False)
    ax.set_facecolor('white')
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_color(OUTLINE_COLOR)
    ax.grid(True, which='major', color=GRID_COLOR, linewidth=_pt(1))
    if pf._xgrid.grid_line_color is None:
        ax.xaxis.grid(False)
    ax.minorticks_on()
    ax.tick_params(which='both', direction='out', labelsize=_pt(TICK_LABEL_PX))
    ax.tick_params(which='major', length=_pt(6))
    ax.tick_params(which='minor', length=_pt(4))
    ax.tick_params(top=False, right=False, which='both')


def _label_kw():
    return dict(fontsize=_pt(AXIS_LABEL_PX), fontstyle='italic',
                parse_math=False)


def _draw(ax, pf, cell_w, scale, measure):
    if pf._y_log:
        ax.set_yscale('log', nonpositive='mask')
    _style_axes(ax, pf)

    axes = {'default': ax}
    twin = None
    extra = pf._extra_x_ranges
    if len(extra) > 1:
        kp._unsupported("more than one extra_x_ranges entry")
    if extra:
        name, xrng2 = next(iter(extra.items()))
        twin = ax.twiny()
        twin.set_autoscale_on(False)
        axes[name] = twin
        for side in ('left', 'right', 'bottom'):
            twin.spines[side].set_visible(False)
        twin.spines['top'].set_color('black')
        twin.minorticks_on()
        twin.tick_params(which='both', direction='out',
                         labelsize=_pt(TICK_LABEL_PX))
        twin.tick_params(which='major', length=_pt(6))
        twin.tick_params(which='minor', length=_pt(4))
        twin.grid(False)
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
            twin.set_xlabel(axis.axis_label, **_label_kw())
    else:
        for lay in pf._layouts:
            kp._unsupported("add_layout(%r) with no extra_x_ranges" % (lay,))

    xb = _Bounds()
    yb = _Bounds(log=pf._y_log)
    legend = {}                         # label -> [artists], first-use order
    for i, g in enumerate(pf._glyphs):
        target = axes.get(g.get('x_range_name', 'default'))
        if target is None:
            kp._unsupported("%s(x_range_name=%r): no such range"
                            % (g['kind'], g.get('x_range_name')))
            continue
        try:
            artist, xs, ys = _GLYPHS[g['kind']](target, g)
        except Exception as e:
            kp._failed("drawing %s()" % g['kind'], e)
            continue
        artist.set_zorder(2 + i * 1e-4)     # Bokeh draws in call order
        if target is ax:
            xb.add(*xs)
        yb.add(*ys)
        label = g.get('legend_label')
        if label:
            legend.setdefault(label, []).append(artist)

    ax.set_xlim(*_limits(pf._x_range, xb, False))
    ax.set_ylim(*_limits(pf._y_range, yb, pf._y_log))
    if twin is not None:
        twin.set_xlim(*_limits(xrng2, _Bounds(), False))

    ax.set_xlabel(pf._x_axis_label, **_label_kw())
    ax.set_ylabel(pf._y_axis_label, **_label_kw())
    title = ax.set_title(pf._title.text, loc='left', fontweight='bold',
                         fontsize=_pt(TITLE_PX), parse_math=False)
    _fit_title(title, cell_w, scale, measure)

    if legend:
        top = twin if twin is not None else ax
        handles = [a[0] if len(a) == 1 else tuple(a) for a in legend.values()]
        leg = top.legend(handles, list(legend.keys()),
                         loc=kp.LEGEND_LOCATIONS[pf._legend.location],
                         fontsize=_pt(LEGEND_PX), frameon=True,
                         fancybox=False, framealpha=0.95, facecolor='white',
                         edgecolor=OUTLINE_COLOR)
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
# Glyphs: each returns (artist, x arrays for bounds, y arrays for bounds)
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
    dash = g.get('line_dash', 'solid')
    if dash == 'solid' or lw <= 0:
        ls = '-'
    else:
        pattern = [_pt(d) for d in _DASHES[dash]]
        if mpl.rcParams['lines.scale_dashes']:
            pattern = [d / lw for d in pattern]
        ls = (0, tuple(pattern))
    line, = ax.plot(g['x'], g['y'], color=color, linewidth=lw, linestyle=ls,
                    marker='None', solid_capstyle='butt',
                    dash_capstyle='butt', solid_joinstyle='bevel',
                    dash_joinstyle='bevel')
    return line, (g['x'],), (g['y'],)


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
    return line, (g['x'],), (g['y'],)


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
    ax.add_collection(pc, autolim=False)
    return pc, (left, right), (bottom, top)


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
    return im, ([x0, x1],), ([y0, y1],)


_GLYPHS = {'line': _line, 'circle': _marker, 'diamond': _marker,
           'scatter': _marker, 'quad': _quad, 'image': _image}
