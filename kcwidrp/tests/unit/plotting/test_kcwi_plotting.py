"""Tests for the Matplotlib plotting wrapper ``kcwidrp.core.kcwi_plotting``
(bokeh plan task 1.2; upgrade plan s7 amendments 5, 8 and 10).

Written from the specification only -- docs/bokeh_port_plan.md ("Core design
decision", "1.1 Wrapper", "1.2 Wrapper tests", "Interfaces fixed for parallel
implementation") and docs/plot_inventory.md section 2 -- never from the
wrapper's source.  The AST coverage scan (the other half of task 1.2) is in
``test_plot_coverage_scan.py``.

How the tests look at a plot
----------------------------
``render(obj)`` returns a Matplotlib ``Figure``; the tests assert on its
axes (lines, markers, rectangles, images, limits, labels, title, legend).
They deliberately accept either of the usual ways to draw a Bokeh glyph
(``Line2D`` or a collection for markers; ``Rectangle`` patches or a
``PolyCollection`` for ``quad``) so that they pin the behaviour, not the
implementation.  Where Matplotlib offers no public accessor the tests use
one private attribute and say so (``Axes._left_title``, ``Legend._loc``).

Units and sources of expected values
------------------------------------
* ``PT_PER_PX`` -- Bokeh sizes are pixels, Matplotlib sizes are points; at
  100 dpi one pixel is 0.72 point.
  # source: analytic - 72 pt per inch / 100 px per inch
* Bokeh defaults (colour ``#1f77b4``, line width 1 px, circle size 4 px,
  10 % ``DataRange1d`` padding) are the values the plan's "Appearance
  defaults" paragraph names.
  # source: published - docs/bokeh_port_plan.md section 1.1, Bokeh 2.4 docs
* Every expected value below carries its own ``# source:`` comment.

All tests run with ``STRICT`` true unless they are in the tolerant group.
Nothing here needs a display, a network, a browser or Bokeh.
"""
import io
import json
import logging
import os
import re
import subprocess
import sys
import threading
from pathlib import Path

import matplotlib.colors as mcolors
import numpy as np
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.collections import PathCollection, PolyCollection
from matplotlib.legend import Legend
from matplotlib.markers import MarkerStyle
from matplotlib.patches import Rectangle
from PIL import Image

from kcwidrp.core import kcwi_plotting as kp

REPO = Path(__file__).resolve().parents[4]
THIS_FILE = os.path.basename(__file__)

# source: analytic - 72 points per inch, 100 pixels per inch (plan 1.1,
# "Marker size units" and "Line width units")
PT_PER_PX = 0.72

# source: published - Bokeh palettes.Spectral11 (ColorBrewer Spectral, 11
# classes, in Bokeh's order: first #5e4fa2, last #9e0142); plan 1.1 trap
# "Spectral11 is eleven discrete colours"
SPECTRAL11 = ['#5e4fa2', '#3288bd', '#66c2a5', '#abdda4', '#e6f598',
              '#ffffbf', '#fee08b', '#fdae61', '#f46d43', '#d53e4f',
              '#9e0142']

# source: published - plan 1.1 "Appearance defaults": default colour
DEFAULT_COLOUR = '#1f77b4'

X = np.arange(11.0)
Y = 10.0 * X


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(autouse=True)
def strict_mode(monkeypatch):
    """The whole suite runs with STRICT true (plan 1.2)."""
    monkeypatch.setattr(kp, 'STRICT', True)


@pytest.fixture
def tolerant(monkeypatch):
    monkeypatch.setattr(kp, 'STRICT', False)


@pytest.fixture
def in_tmp(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


@pytest.fixture
def kcwi_log(caplog, monkeypatch):
    """Capture the 'KCWI' logger (it may not propagate in a reduction)."""
    lg = logging.getLogger('KCWI')
    monkeypatch.setattr(lg, 'propagate', False)
    old_level = lg.level
    lg.setLevel(logging.DEBUG)
    lg.addHandler(caplog.handler)
    caplog.set_level(logging.DEBUG)
    yield caplog
    lg.removeHandler(caplog.handler)
    lg.setLevel(old_level)


def warnings_of(caplog):
    return [r for r in caplog.records
            if r.name == 'KCWI' and r.levelno >= logging.WARNING]


# ---------------------------------------------------------------------------
# builders and readers
# ---------------------------------------------------------------------------
def fig_(**kw):
    kw.setdefault('title', 'T')
    kw.setdefault('x_axis_label', 'X')
    kw.setdefault('y_axis_label', 'Y')
    kw.setdefault('plot_width', 400)
    kw.setdefault('plot_height', 300)
    return kp.figure(**kw)


def draw(f):
    """Force layout and drawing; returns the figure."""
    f.savefig(io.BytesIO(), format='png')
    return f


def rendered(obj, **kw):
    return draw(kp.render(obj, **kw))


def arr(a):
    return np.ma.filled(np.ma.asarray(a, dtype=float), np.nan)


def rgba(c, alpha=None):
    r = mcolors.to_rgba(c)
    return r if alpha is None else r[:3] + (alpha,)


def eff_alpha(artist, colour):
    a = artist.get_alpha()
    return colour[3] if a is None else a * colour[3]


def _no_marker(m):
    return m in (None, 'None', 'none', '', ' ')


def glyph_lines(ax):
    """Line2D artists that are drawn as lines (the ``line`` glyph)."""
    return [l for l in ax.lines if _no_marker(l.get_marker())
            and l.get_linestyle() not in ('None', 'none', '', ' ')]


def line_info(l):
    c = rgba(l.get_color())
    return dict(x=arr(l.get_xdata()), y=arr(l.get_ydata()), colour=c,
                alpha=eff_alpha(l, c), lw=l.get_linewidth(),
                ls=l.get_linestyle())


_KIND_CACHE = {}


def _kind_of_path(path):
    for m in ('o', 'x', 'D'):
        if m not in _KIND_CACHE:
            ms = MarkerStyle(m)
            _KIND_CACHE[m] = ms.get_path().transformed(ms.get_transform())
        ref = _KIND_CACHE[m]
        if ref.vertices.shape == path.vertices.shape and np.allclose(
                ref.vertices, path.vertices):
            return m
    return '?'


def markers(ax):
    """One dict per marker artist (Line2D markers or scatter collections)."""
    out = []
    for l in ax.lines:
        if _no_marker(l.get_marker()):
            continue
        m = l.get_marker()
        kind = {'d': 'D'}.get(m, m)
        line_c = rgba(l.get_color())
        mfc = l.get_markerfacecolor()
        mec = l.get_markeredgecolor()
        filled = MarkerStyle(m).is_filled()
        face = rgba(line_c if mfc == 'auto' and filled else
                    ('none' if mfc == 'auto' else mfc))
        edge = rgba(line_c if mec == 'auto' else mec)
        ew = l.get_markeredgewidth()
        a = l.get_alpha()
        out.append(dict(
            kind=kind, x=arr(l.get_xdata()), y=arr(l.get_ydata()),
            diameter=l.get_markersize(), face=face, edge=edge, ew=ew,
            alpha=a, has_edge=edge[3] > 0 and ew > 0,
            has_face=face[3] > 0 and filled))
    for c in ax.collections:
        if not isinstance(c, PathCollection):
            continue
        offs = np.asarray(c.get_offsets(), float).reshape(-1, 2)
        fc = c.get_facecolor()
        ec = c.get_edgecolor()
        face = tuple(fc[0]) if len(fc) else (0, 0, 0, 0)
        edge = tuple(ec[0]) if len(ec) else (0, 0, 0, 0)
        lw = float(np.atleast_1d(c.get_linewidth())[0])
        sizes = np.atleast_1d(c.get_sizes())
        out.append(dict(
            kind=_kind_of_path(c.get_paths()[0]), x=offs[:, 0], y=offs[:, 1],
            diameter=float(np.sqrt(sizes[0])), face=face, edge=edge, ew=lw,
            alpha=c.get_alpha(), has_edge=edge[3] > 0 and lw > 0,
            has_face=face[3] > 0))
    return out


def rects(ax):
    """Sorted (left, bottom, right, top) of every rectangle, with colours."""
    out = []
    for p in ax.patches:
        if isinstance(p, Rectangle):
            out.append(dict(
                box=(p.get_x(), p.get_y(), p.get_x() + p.get_width(),
                     p.get_y() + p.get_height()),
                face=rgba(p.get_facecolor()), edge=rgba(p.get_edgecolor()),
                alpha=p.get_alpha()))
    for c in ax.collections:
        if isinstance(c, PolyCollection):
            fcs, ecs = c.get_facecolor(), c.get_edgecolor()
            for i, path in enumerate(c.get_paths()):
                v = path.vertices
                out.append(dict(
                    box=(v[:, 0].min(), v[:, 1].min(), v[:, 0].max(),
                         v[:, 1].max()),
                    face=tuple(fcs[i % len(fcs)]),
                    edge=tuple(ecs[i % len(ecs)]) if len(ecs)
                    else (0, 0, 0, 0), alpha=c.get_alpha()))
    return sorted(out, key=lambda d: d['box'])


def legend_labels(ax):
    leg = ax.get_legend()
    return None if leg is None else [t.get_text() for t in leg.get_texts()]


def legend_loc(ax):
    # Legend has no public accessor for the location it was built with.
    return {v: k for k, v in Legend.codes.items()}[ax.get_legend()._loc]


def all_legends(f):
    return [a.get_legend() for a in f.axes if a.get_legend()] + list(f.legends)


def title_of(ax):
    return ax.get_title(loc='left')


def grid_visible(ax, which):
    axis = ax.xaxis if which == 'x' else ax.yaxis
    return any(g.get_visible() for g in axis.get_gridlines())


def png_size(path):
    with Image.open(path) as im:
        return im.size


def png_pixels(path):
    with Image.open(path) as im:
        return np.asarray(im.convert('RGBA'))


def describe(f):
    """A plain-data summary of a rendered figure, for equality checks."""
    out = []
    for ax in f.axes:
        out.append(dict(
            lines=[(tuple(np.nan_to_num(line_info(l)['x'], nan=-9e9)),
                    tuple(np.nan_to_num(line_info(l)['y'], nan=-9e9)),
                    line_info(l)['colour'], line_info(l)['ls'],
                    line_info(l)['lw']) for l in glyph_lines(ax)],
            nmark=len(markers(ax)), nrect=len(rects(ax)),
            xlim=tuple(ax.get_xlim()), ylim=tuple(ax.get_ylim()),
            title=title_of(ax), xlabel=ax.get_xlabel(),
            ylabel=ax.get_ylabel(), legend=legend_labels(ax),
            loc=legend_loc(ax) if ax.get_legend() else None,
            xgrid=grid_visible(ax, 'x'), ygrid=grid_visible(ax, 'y'),
            yscale=ax.get_yscale()))
    return out


def solvearcs_like(order='main_first'):
    """The SolveArcs twin-axis figure: 3 glyphs on the main x range, 5 on
    the 'pix' range, one legend of six labels ('Arc' and 'Atlas' twice)."""
    p = kp.figure(title='ATLAS/ARC LINE FITS', x_axis_label='Wavelength (A)',
                  y_axis_label='Flux', plot_width=1000, plot_height=600)
    p.x_range = kp.Range1d(start=4000., end=4100.)
    p.extra_x_ranges = {'pix': kp.Range1d(start=10., end=60.)}
    p.add_layout(kp.LinearAxis(x_range_name='pix', axis_label='CCD Y pix'),
                 'above')
    w = np.linspace(4000., 4100., 51)
    pix = np.linspace(10., 60., 51)
    f1 = np.sin(w / 7.) ** 2
    f2 = np.cos(w / 11.) ** 2

    def main():
        p.line(w, f1, color='black', legend_label='Arc')
        p.line(w, f2, color='red', legend_label='Atlas')
        p.circle(w[::5], f1[::5], color='blue', legend_label='Fit')

    def twin():
        p.line(pix, f1, color='black', legend_label='Arc',
               x_range_name='pix', line_dash='dashed')
        p.line(pix, f2, color='red', legend_label='Atlas',
               x_range_name='pix', line_dash='dashdot')
        p.circle(pix[::5], f1[::5], color='green', legend_label='Peaks',
                 x_range_name='pix')
        p.circle(pix[::7], f2[::7], color='orange', legend_label='Cat lines',
                 x_range_name='pix')
        p.circle(pix[::9], f2[::9], color='purple', legend_label='Resid',
                 x_range_name='pix')

    if order == 'main_first':
        main()
        twin()
    else:
        twin()
        main()
    return p


# ===========================================================================
# 1. Module surface and appearance defaults
# ===========================================================================
def test_module_exposes_the_fixed_interface():
    # source: published - plan 0.3 "Interfaces fixed for parallel
    # implementation": STRICT, figure/PlotFigure, Range1d, LinearAxis,
    # gridplot/PlotGrid, render, pixel_size, snapshot, plot_label,
    # save_plot, get_plot_lims, set_plot_lims, oplot_slices
    for name in ('STRICT', 'figure', 'PlotFigure', 'Range1d', 'LinearAxis',
                 'gridplot', 'PlotGrid', 'render', 'pixel_size', 'snapshot',
                 'plot_label', 'save_plot', 'get_plot_lims', 'set_plot_lims',
                 'oplot_slices'):
        assert hasattr(kp, name), name


def test_figure_returns_a_recorder_and_gridplot_a_grid():
    p = fig_()
    # source: published - plan 0.3 interfaces: figure -> PlotFigure,
    # gridplot -> PlotGrid
    assert isinstance(p, kp.PlotFigure)
    g = kp.gridplot([p], ncols=1, plot_width=100, plot_height=100)
    assert isinstance(g, kp.PlotGrid)


def test_figure_title_text_is_readable():
    # source: published - plan 0.3: "p.title.text is readable"
    assert fig_(title='Hello').title.text == 'Hello'


def test_title_labels_and_size_of_a_single_figure():
    p = fig_(title='My title', x_axis_label='lam', y_axis_label='flux',
             plot_width=400, plot_height=300)
    f = rendered(p)
    ax = f.axes[0]
    # source: published - plan 1.1 table, figure(title, x_axis_label,
    # y_axis_label, plot_width, plot_height); title is the left title
    assert len(f.axes) == 1
    assert title_of(ax) == 'My title'
    assert ax.get_title() == ''
    assert ax.get_xlabel() == 'lam'
    assert ax.get_ylabel() == 'flux'
    # source: published - plan 1.1: "size maps to pixels at 100 dpi"
    assert f.dpi == 100
    assert tuple(f.get_size_inches() * f.dpi) == pytest.approx((400, 300))
    assert kp.pixel_size(p) == (400, 300)


def test_title_is_left_aligned_and_bold():
    ax = rendered(fig_(title='abc')).axes[0]
    # source: published - plan 1.1 appearance defaults: "left-aligned bold
    # title".  _left_title is private: there is no public accessor for the
    # left title's font weight.
    w = ax._left_title.get_fontweight()
    assert w == 'bold' or (isinstance(w, (int, float)) and w >= 700)


def test_grid_is_light_grey_by_default():
    ax = rendered(fig_()).axes[0]
    ax.figure.canvas.draw()
    lines = [g for g in ax.xaxis.get_gridlines() + ax.yaxis.get_gridlines()
             if g.get_visible()]
    # source: published - plan 1.1 appearance defaults: "light grey grid"
    assert grid_visible(ax, 'x') and grid_visible(ax, 'y')
    for g in lines:
        r, gr, b, _ = rgba(g.get_color())
        assert abs(r - gr) < 0.02 and abs(gr - b) < 0.02   # grey
        # Tolerance: "light" is not numeric in the spec; > 0.7 of full
        # scale excludes Matplotlib's dark defaults and any black grid.
        assert r > 0.7


def test_tooltips_are_accepted_and_ignored():
    p = fig_(tooltips=[("x", "$x"), ("y", "$y")])
    p.line(X, Y)
    q = fig_()
    q.line(X, Y)
    # source: published - plan 1.1 table: "tooltips accepted and ignored"
    assert describe(rendered(p)) == describe(rendered(q))


# ===========================================================================
# 2. Glyphs, colours, units
# ===========================================================================
def test_line_data_and_default_colour():
    p = fig_()
    p.line(X, Y)
    ax = rendered(p).axes[0]
    ls = glyph_lines(ax)
    # source: published - plan 1.1 table (`line`) and appearance defaults
    assert len(ls) == 1
    assert np.array_equal(line_info(ls[0])['x'], X)
    assert np.array_equal(line_info(ls[0])['y'], Y)
    assert mcolors.to_hex(line_info(ls[0])['colour']) == DEFAULT_COLOUR


def test_no_colour_cycling():
    p = fig_()
    p.line(X, Y)
    p.line(X, Y + 5)
    p.line(X, Y + 10)
    p.circle(X, Y + 15)
    p.scatter(X, Y + 20)
    ax = rendered(p).axes[0]
    cols = [mcolors.to_hex(line_info(l)['colour']) for l in glyph_lines(ax)]
    cols += [mcolors.to_hex(m['face']) for m in markers(ax)]
    cols += [mcolors.to_hex(m['edge']) for m in markers(ax)]
    # source: published - plan 1.1 appearance defaults: "default glyph colour
    # #1f77b4 with no colour cycling"
    assert cols == [DEFAULT_COLOUR] * 7


@pytest.mark.parametrize("spec", ['red', '#00ff00', 'black', 'purple'])
def test_line_colour_keyword(spec):
    p = fig_()
    p.line(X, Y, color=spec)
    c = line_info(glyph_lines(rendered(p).axes[0])[0])['colour']
    # source: published - plan 1.1 table, glyph keyword `color` (colour names
    # and hex strings mean what they mean in CSS / Matplotlib)
    assert c == mcolors.to_rgba(spec)


def test_line_color_keyword_sets_the_line():
    p = fig_()
    p.line(X, Y, line_color='navy')
    c = line_info(glyph_lines(rendered(p).axes[0])[0])['colour']
    # source: published - plan 1.1 table, glyph keyword `line_color`
    assert c == mcolors.to_rgba('navy')


def test_trap_pixel_line_width():
    p = fig_()
    p.line(X, Y)
    p.line(X, Y + 1, line_width=2)
    p.line(X, Y + 2, line_width=4)
    ls = [line_info(l)['lw'] for l in glyph_lines(rendered(p).axes[0])]
    # source: analytic - Bokeh line_width is pixels (default 1);
    # Matplotlib points = pixels * 0.72 at 100 dpi (plan 1.1 trap "Line
    # width units")
    # Tolerance: single float multiplication, relative 1e-15.
    assert ls == pytest.approx([1 * PT_PER_PX, 2 * PT_PER_PX, 4 * PT_PER_PX],
                               rel=1e-9)


def test_line_alpha_keyword():
    p = fig_()
    p.line(X, Y, alpha=0.7)
    li = line_info(glyph_lines(rendered(p).axes[0])[0])
    # source: published - plan 1.1 table, glyph keyword `alpha` (used with
    # line_width=4 in MakeMasterBias)
    assert li['alpha'] == pytest.approx(0.7, rel=1e-12)


def test_line_alpha_keyword_line_alpha():
    p = fig_()
    p.line(X, Y, line_alpha=0.25)
    li = line_info(glyph_lines(rendered(p).axes[0])[0])
    # source: published - plan 1.1 table, glyph keyword `line_alpha`
    assert li['alpha'] == pytest.approx(0.25, rel=1e-12)


def dash_shape(line):
    """The on/off dash sequence of a Line2D relative to its first dash, or
    None for a solid line (independent of line width and units)."""
    seq = line._unscaled_dash_pattern[1]
    return None if not seq else [v / seq[0] for v in seq]


# source: published - Bokeh 2.4.3 dash patterns (core/property/
# visual.py DashPattern): dashed [6, 6], dotted [2, 4], dashdot [6, 4, 2, 4]
BOKEH_DASH = {'dashed': [1., 1.], 'dotted': [1., 2.],
              'dashdot': [1., 4 / 6, 2 / 6, 4 / 6]}


@pytest.mark.parametrize("name", ['dashed', 'dashdot', 'dotted'])
def test_trap_dash_names(name):
    p = fig_()
    p.line(X, Y)
    p.line(X, Y + 1, line_dash=name)
    got = [dash_shape(l) for l in glyph_lines(rendered(p).axes[0])]
    # source: published - plan 1.1 trap "Dash names": each name maps to its
    # own pattern; the shape is Bokeh's (a solid line has no dashes)
    # Tolerance: ratios of float64 lengths, relative 1e-9.
    assert got[0] is None
    assert got[1] == pytest.approx(BOKEH_DASH[name], rel=1e-9)


@pytest.mark.parametrize("glyph", ['circle', 'scatter', 'diamond'])
def test_marker_glyph_positions(glyph):
    p = fig_()
    getattr(p, glyph)(X, Y)
    ms = markers(rendered(p).axes[0])
    # source: published - plan 1.1 table: the six glyph methods
    assert len(ms) == 1
    assert np.array_equal(ms[0]['x'], X) and np.array_equal(ms[0]['y'], Y)


def test_marker_shapes():
    p = fig_()
    p.circle(X, Y)
    p.scatter(X, Y + 1)
    p.scatter(X, Y + 2, marker='x')
    p.diamond(X, Y + 3)
    ms = markers(rendered(p).axes[0])
    # source: published - plan 1.1 table: `scatter` with no marker is a
    # circle (FitCenter); 'x' is the only named marker
    assert [m['kind'] for m in ms][:3] == ['o', 'o', 'x']
    # source: published - Bokeh 2.4.3 diamond marker (markers/defs.ts): a
    # rhombus as tall as `size` and 2/3 as wide, not Matplotlib's square 'D'
    # Tolerance: path vertices are exact fractions, absolute 1e-12.
    verts = np.asarray(ms[3]['kind'].vertices)
    w = verts[:, 0].max() - verts[:, 0].min()
    h = verts[:, 1].max() - verts[:, 1].min()
    assert w / h == pytest.approx(2 / 3, abs=1e-12)
    unit = verts / (h / 2)
    for corner in [(0., 1.), (2 / 3, 0.), (0., -1.), (-2 / 3, 0.)]:
        assert np.min(np.hypot(*(unit - corner).T)) < 1e-12, corner


def test_circle_and_scatter_default_size_is_four_pixels():
    p = fig_()
    p.circle(X, Y)
    p.scatter(X, Y + 1)
    d = [m['diameter'] for m in markers(rendered(p).axes[0])]
    # source: published - plan 1.1 appearance defaults: "marker size
    # equivalent to Bokeh's 4 px"; analytic: 4 px = 2.88 pt
    # Tolerance: float multiplication / square root, relative 1e-12.
    assert d == pytest.approx([4 * PT_PER_PX] * 2, rel=1e-9)


@pytest.mark.parametrize("glyph, kw", [
    ('circle', {}), ('scatter', {}), ('scatter', {'marker': 'x'})])
def test_trap_one_pixel_markers(glyph, kw):
    p = fig_()
    getattr(p, glyph)(X, Y, size=1, **kw)
    m = markers(rendered(p).axes[0])[0]
    # source: analytic - plan 1.1 trap "Marker size units": Bokeh size is a
    # diameter in pixels; 1 px = 0.72 pt at 100 dpi, so size=1 stays one
    # pixel wide.  Tolerance: float multiplication / square root, 1e-12.
    assert m['diameter'] == pytest.approx(1 * PT_PER_PX, rel=1e-9)


def test_marker_size_is_linear_in_pixels():
    p = fig_()
    p.circle(X, Y, size=2)
    p.circle(X, Y + 1, size=8)
    p.diamond(X, Y + 2, size=8)
    p.diamond(X, Y + 3, size=16)
    d = [m['diameter'] for m in markers(rendered(p).axes[0])]
    # source: analytic - sizes 2 and 8 px -> 1.44 and 5.76 pt; diamond sizes
    # 8 and 16 differ by exactly a factor 2 (the absolute diamond scale is
    # a Matplotlib marker convention, so only the ratio is asserted)
    # Tolerance: float arithmetic, relative 1e-12.
    assert d[0] == pytest.approx(2 * PT_PER_PX, rel=1e-9)
    assert d[1] == pytest.approx(8 * PT_PER_PX, rel=1e-9)
    assert d[3] / d[2] == pytest.approx(2.0, rel=1e-9)


def test_filled_marker_with_zero_line_alpha_has_no_edge():
    p = fig_()
    p.circle(X, Y, size=1, line_alpha=0., fill_color='purple')
    m = markers(rendered(p).axes[0])[0]
    # source: published - plan 1.1 table: "line_alpha=0. with fill_color= on
    # circle ... filled marker with no outline ... not an invisible marker"
    assert m['has_face']
    assert mcolors.to_hex(m['face'][:3]) == mcolors.to_hex('purple')
    assert not m['has_edge']


def test_color_on_a_marker_sets_outline_and_fill():
    p = fig_()
    p.circle(X, Y, color='red')
    m = markers(rendered(p).axes[0])[0]
    # source: published - plan 1.1 table: "color= on a marker glyph sets
    # both outline and fill, as in Bokeh"
    assert mcolors.to_hex(m['face'][:3]) == '#ff0000'
    assert mcolors.to_hex(m['edge'][:3]) == '#ff0000'
    assert m['has_face'] and m['has_edge']


def test_line_color_and_fill_color_override_each_part():
    p = fig_()
    p.circle(X, Y, color='red', fill_color='blue')
    p.circle(X, Y + 1, color='red', line_color='green')
    a, b = markers(rendered(p).axes[0])
    # source: published - plan 1.1 table: "`line_color` and `fill_color`
    # override each part"
    assert mcolors.to_hex(a['face'][:3]) == '#0000ff'
    assert mcolors.to_hex(a['edge'][:3]) == '#ff0000'
    assert mcolors.to_hex(b['face'][:3]) == '#ff0000'
    assert mcolors.to_hex(b['edge'][:3]) == '#008000'


def test_circle_fill_color_only():
    p = fig_()
    p.circle(X, Y, fill_color='orange')
    m = markers(rendered(p).axes[0])[0]
    # source: published - inventory 2.2: circle(fill_color=...) (MakeMasterFlat)
    assert mcolors.to_hex(m['face'][:3]) == mcolors.to_hex('orange')


def test_marker_x_is_drawn_with_the_colour_given():
    p = fig_()
    p.scatter(X, Y, marker='x', color='red', legend_label='bad')
    m = markers(rendered(p).axes[0])[0]
    # source: published - inventory 2.2: scatter(marker='x', color=..)
    assert m['kind'] == 'x'
    assert mcolors.to_hex(m['edge'][:3]) == '#ff0000'


def test_quad_rectangles():
    left = np.array([0., 1., 2.])
    right = np.array([1., 2., 3.])
    top = np.array([1., 3., 2.])
    p = fig_()
    p.quad(top=top, bottom=0, left=left, right=right,
           fill_color='gray', line_color='black', alpha=0.5)
    ax = rendered(p).axes[0]
    rs = rects(ax)
    # source: published - plan 1.1 table: quad(top=, bottom=, left=, right=),
    # array or scalar per edge, rectangles with fill_color, line_color, alpha
    assert [r['box'] for r in rs] == [(0, 0, 1, 1), (1, 0, 2, 3), (2, 0, 3, 2)]
    assert mcolors.to_hex(rs[0]['face'][:3]) == mcolors.to_hex('gray')
    assert mcolors.to_hex(rs[0]['edge'][:3]) == '#000000'
    # source: published - Bokeh alpha applies to the fill (and line); the
    # RGBA of a Matplotlib patch already carries the alpha
    assert rs[0]['face'][3] == pytest.approx(0.5, rel=1e-12)
    # source: analytic - x data 0..3, y data 0..3, Bokeh's padding of 5 % of
    # the span 3 at each end (Bokeh 2.4.3 data_range1d.js)
    # Tolerance: float64 arithmetic, relative 1e-12.
    assert ax.get_xlim() == pytest.approx((-0.15, 3.15), rel=1e-9)
    assert ax.get_ylim() == pytest.approx((-0.15, 3.15), rel=1e-9)


def test_quad_scalar_edges_broadcast():
    p = fig_()
    p.quad(top=2., bottom=1., left=[0., 5.], right=[1., 6.])
    # source: published - plan 1.1 table: "array or scalar per edge"
    assert [r['box'] for r in rects(rendered(p).axes[0])] == [
        (0, 1, 1, 2), (5, 1, 6, 2)]


def test_quad_is_keyword_only():
    p = fig_()
    with pytest.raises(TypeError):
        # source: published - plan 1.1 table: "Keyword-only"
        p.quad([1.], [0.], [0.], [1.])


# ---------------------------------------------------------------------------
# image
# ---------------------------------------------------------------------------
def _image_fig(a, x=1., y=2., dw=4., dh=6., **kw):
    p = fig_(**kw)
    p.image([a], x=x, y=y, dw=dw, dh=dh, palette="Spectral11", level="image")
    return p


def image_value_at(im, xd, yd):
    """Value displayed at data position (xd, yd), from extent and origin."""
    l, r, b, t = im.get_extent()
    a = np.asarray(im.get_array())
    nrow, ncol = a.shape
    col = min(int((xd - l) / (r - l) * ncol), ncol - 1)
    frow = min(int((yd - b) / (t - b) * nrow), nrow - 1)
    row = frow if im.origin == 'lower' else nrow - 1 - frow
    return a[row, col]


def test_trap_image_orientation_and_extent():
    # an asymmetric array: row 0 is [1, 2, 3], the last row is [7, 8, 9]
    a = np.array([[1., 2., 3.], [4., 5., 6.], [7., 8., 9.]])
    ax = rendered(_image_fig(a, x=1., y=2., dw=4., dh=6.)).axes[0]
    im = ax.images[0]
    # source: published - plan 1.1 trap "Image orientation": Bokeh puts row 0
    # at the bottom; extent (x, x + dw, y, y + dh).
    assert len(ax.images) == 1
    assert tuple(im.get_extent()) == (1., 5., 2., 8.)
    # the cell at the bottom-left corner of the extent is a[0, 0], the
    # top-left is a[-1, 0], the bottom-right a[0, -1], the top-right a[-1,-1]
    # source: published - Bokeh image() semantics (row 0 at y)
    assert image_value_at(im, 1.1, 2.1) == a[0, 0]
    assert image_value_at(im, 1.1, 7.9) == a[-1, 0]
    assert image_value_at(im, 4.9, 2.1) == a[0, -1]
    assert image_value_at(im, 4.9, 7.9) == a[-1, -1]
    assert im.origin == 'lower'
    # source: published - plan 1.1 trap: "Use origin='lower'"


def test_image_fills_the_axes_and_can_have_zero_padding():
    a = np.arange(12.).reshape(3, 4)
    p = _image_fig(a, x=1., y=2., dw=4., dh=6.)
    p.x_range.range_padding = p.y_range.range_padding = 0
    ax = rendered(p).axes[0]
    # source: published - plan 1.1 table: range_padding = 0 gives zero
    # margins on that axis (MakeCube); the limits are the image extent
    assert ax.get_xlim() == pytest.approx((1., 5.), rel=1e-12)
    assert ax.get_ylim() == pytest.approx((2., 8.), rel=1e-12)
    # source: analytic - an image must fill the plot area, not keep square
    # pixels (dw != dh here)
    assert ax.get_aspect() == 'auto'


def test_trap_spectral11_endpoints_and_discreteness():
    a = np.linspace(-3., 7., 40).reshape(5, 8)
    im = rendered(_image_fig(a)).axes[0].images[0]
    cmap = im.cmap
    vmin, vmax = im.get_clim()
    # source: published - plan 1.1 trap "Spectral11": eleven discrete
    # colours, #5e4fa2 first and #9e0142 last, mapped linearly from the
    # array minimum to its maximum
    assert (vmin, vmax) == (a.min(), a.max())
    assert mcolors.to_hex(cmap(im.norm(a.min()))) == SPECTRAL11[0]
    assert mcolors.to_hex(cmap(im.norm(a.max()))) == SPECTRAL11[-1]
    ramp = cmap(np.linspace(0., 1., 2000))
    assert len({tuple(c) for c in ramp}) == 11
    # each of eleven equal-width bins of the data range has its own colour,
    # in Bokeh's order
    # source: published - Bokeh LinearColorMapper: 11 equal bins low..high
    centres = vmin + (np.arange(11) + 0.5) / 11. * (vmax - vmin)
    assert [mcolors.to_hex(cmap(im.norm(c))) for c in centres] == SPECTRAL11


def test_spectral11_is_not_matplotlibs_continuous_spectral():
    a = np.linspace(0., 1., 30).reshape(5, 6)
    cmap = rendered(_image_fig(a)).axes[0].images[0].cmap
    import matplotlib
    spectral = matplotlib.colormaps['Spectral']
    # source: published - plan 1.1 trap: Matplotlib's Spectral is
    # continuous and runs in the opposite direction
    assert mcolors.to_hex(cmap(0.0)) != mcolors.to_hex(spectral(0.0))
    assert mcolors.to_hex(cmap(0.0)) == SPECTRAL11[0]


def test_image_level_and_tooltips_accepted():
    p = fig_(tooltips=[("x", "$x"), ("y", "$y"), ("value", "@image")])
    p.image([np.arange(6.).reshape(2, 3)], x=0, y=0, dw=3, dh=2,
            palette="Spectral11", level="image")
    # source: published - plan 1.1 table: `level` accepted and ignored
    assert len(rendered(p).axes[0].images) == 1


# ===========================================================================
# 3. Ranges, limits, padding
# ===========================================================================
@pytest.mark.parametrize("glyph", ['line', 'circle', 'scatter', 'diamond'])
def test_default_limits_pad_ten_percent_like_bokeh(glyph):
    p = fig_()
    getattr(p, glyph)(X, Y)
    ax = rendered(p).axes[0]
    # source: published - Bokeh 2.4.3 DataRange1d (data_range1d.js):
    # span = (max - min) * (1 + range_padding), range_padding = 0.1, centred,
    # so 5 % of the span is added at each end.
    # x span 10 -> (-0.5, 10.5); y span 100 -> (-5, 105)
    # Tolerance: float64 arithmetic on values ~1e2, relative 1e-12.
    assert ax.get_xlim() == pytest.approx((-0.5, 10.5), rel=1e-9)
    assert ax.get_ylim() == pytest.approx((-5., 105.), rel=1e-9)


def test_limits_cover_all_glyphs_on_the_axes():
    p = fig_()
    p.line(X, Y)
    p.circle([20.], [-50.])
    ax = rendered(p).axes[0]
    # source: analytic - x data 0..20 (pad 1), y data -50..100 (pad 7.5),
    # 5 % of the span at each end as Bokeh's DataRange1d
    # Tolerance: float64 arithmetic, relative 1e-12.
    assert ax.get_xlim() == pytest.approx((-1., 21.), rel=1e-9)
    assert ax.get_ylim() == pytest.approx((-57.5, 107.5), rel=1e-9)


@pytest.mark.parametrize("make", [
    lambda: [2., 8.], lambda: (2., 8.), lambda: kp.Range1d(2., 8.),
    lambda: kp.Range1d(start=2., end=8.)],
    ids=['list', 'tuple', 'Range1d-positional', 'Range1d-keyword'])
def test_figure_y_range_forms(make):
    p = fig_(y_range=make())
    p.line(X, Y)
    ax = rendered(p).axes[0]
    # source: published - plan 1.1 table: figure(y_range=[lo, hi]) accepts a
    # list, tuple or Range1d and sets fixed y limits; x stays automatic
    assert ax.get_ylim() == (2., 8.)
    assert ax.get_xlim() == pytest.approx((-0.5, 10.5), rel=1e-9)


def test_range1d_both_forms():
    a = kp.Range1d(1., 5.)
    b = kp.Range1d(start=1., end=5.)
    # source: published - plan 1.1 table: "Range1d(a, b) and
    # Range1d(start=a, end=b)"
    assert (a.start, a.end) == (1., 5.)
    assert (b.start, b.end) == (1., 5.)


def test_assigned_ranges_replace_the_range_object():
    p = fig_()
    p.line(X, Y)
    p.x_range = kp.Range1d(3., 4.)
    p.y_range = kp.Range1d(start=-5., end=500.)
    ax = rendered(p).axes[0]
    # source: published - plan 1.1 table: "p.x_range = Range1d(...) replace
    # the range object"
    assert ax.get_xlim() == (3., 4.)
    assert ax.get_ylim() == (-5., 500.)


def test_range_start_and_end_are_settable_on_their_own():
    p = fig_()
    p.line(X, Y)
    p.y_range.start = 0
    p.x_range.end = 7
    ax = rendered(p).axes[0]
    # source: published - plan 1.1 table: each of x/y_range.start/.end
    # settable on its own, leaving the other end automatic (SolveArcs sets
    # only y_range.start = 0).  The automatic end keeps Bokeh's padding
    # (5 % of the data span at each end, Bokeh 2.4.3 data_range1d.js):
    # y end = 100 + 5 = 105; x start = 0 - 0.5 = -0.5.
    # Tolerance: float64 arithmetic, relative 1e-12.
    assert ax.get_ylim()[0] == 0
    assert ax.get_ylim()[1] == pytest.approx(105., rel=1e-9)
    assert ax.get_xlim()[1] == 7
    assert ax.get_xlim()[0] == pytest.approx(-0.5, rel=1e-9)


def test_all_four_range_ends_set_like_set_plot_lims():
    p = fig_()
    p.line(X, Y)
    p.x_range.start = 2.
    p.x_range.end = 4.
    p.y_range.start = 20.
    p.y_range.end = 40.
    ax = rendered(p).axes[0]
    # source: published - plan 1.1 table: set_plot_lims assigns all four
    assert ax.get_xlim() == (2., 4.)
    assert ax.get_ylim() == (20., 40.)


def test_replaced_range_can_be_modified():
    p = fig_()
    p.line(X, Y)
    p.x_range = kp.Range1d(0., 10.)
    p.x_range.end = 5.
    # source: analytic - a Range1d is a plain mutable pair
    assert rendered(p).axes[0].get_xlim() == (0., 5.)


def test_range_padding_zero_chained_assignment():
    p = fig_()
    p.line(X, Y)
    p.x_range.range_padding = p.y_range.range_padding = 0
    ax = rendered(p).axes[0]
    # source: published - plan 1.1 table: chained assignment in MakeCube,
    # "Zero margins on that axis"; limits are exactly the data range
    assert ax.get_xlim() == (0., 10.)
    assert ax.get_ylim() == (0., 100.)


def test_range_padding_zero_on_one_axis_only():
    p = fig_()
    p.line(X, Y)
    p.x_range.range_padding = 0
    ax = rendered(p).axes[0]
    # source: published - plan 1.1 table; y keeps Bokeh's default padding
    # (5 % of the span at each end, Bokeh 2.4.3 data_range1d.js)
    assert ax.get_xlim() == (0., 10.)
    assert ax.get_ylim() == pytest.approx((-5., 105.), rel=1e-9)


def test_set_plot_lims_sets_the_four_ends():
    p = fig_()
    p.line(X, Y)
    kp.set_plot_lims(p, xlim=[2., 6.], ylim=[10., 70.])
    ax = rendered(p).axes[0]
    # source: published - plan 0.3: set_plot_lims unchanged (inventory 2.5:
    # fig.x_range.start/end and fig.y_range.start/end are assigned)
    assert ax.get_xlim() == (2., 6.)
    assert ax.get_ylim() == (10., 70.)


def test_set_plot_lims_with_one_argument_leaves_the_other_automatic():
    p = fig_()
    p.line(X, Y)
    kp.set_plot_lims(p, ylim=[10., 70.])
    ax = rendered(p).axes[0]
    # source: published - inventory 2.5: only the given limit is assigned
    assert ax.get_ylim() == (10., 70.)
    assert ax.get_xlim() == pytest.approx((-0.5, 10.5), rel=1e-9)


def test_set_plot_lims_before_or_after_the_glyphs_is_the_same():
    a = fig_()
    kp.set_plot_lims(a, xlim=[2., 6.], ylim=[10., 70.])
    a.line(X, Y)
    b = fig_()
    b.line(X, Y)
    kp.set_plot_lims(b, xlim=[2., 6.], ylim=[10., 70.])
    # source: published - plan "Core design decision": order independence
    assert describe(rendered(a)) == describe(rendered(b))


def test_oplot_slices_draws_23_black_dashdot_lines():
    p = fig_()
    p.line(X, Y)
    kp.oplot_slices(p, [0., 100.])
    ax = rendered(p).axes[0]
    ls = glyph_lines(ax)[1:]
    # source: published - inventory 2.5: oplot_slices draws 23 vertical
    # black dash-dot lines at x = ix*5 - 0.5 for ix = 1..23 (4.5 ... 114.5)
    assert len(ls) == 23
    xs = [line_info(l)['x'] for l in ls]
    assert [x[0] for x in xs] == [ix * 5 - 0.5 for ix in range(1, 24)]
    assert all(x[0] == x[1] for x in xs)
    assert all(np.array_equal(line_info(l)['y'], [0., 100.]) for l in ls)
    assert all(mcolors.to_hex(line_info(l)['colour']) == '#000000'
               for l in ls)
    # Tolerance: ratios of float64 lengths, relative 1e-9.
    assert all(dash_shape(l) == pytest.approx(BOKEH_DASH['dashdot'],
                                              rel=1e-9) for l in ls)
    assert ax.get_legend() is None


def test_log_y_axis():
    p = fig_(y_axis_type='log')
    p.line(X + 1, 10.0 ** X)
    ax = rendered(p).axes[0]
    # source: published - plan 1.1 table: y_axis_type='log' maps to a log y
    # scale
    assert ax.get_yscale() == 'log'
    assert ax.get_xscale() == 'linear'


def test_log_y_axis_with_non_positive_values_does_not_raise(in_tmp):
    p = fig_(y_axis_type='log')
    p.line(X, np.array([-1., 0., 1., 10., 100., np.nan, 5., 7., 3., 2., 1.]))
    p.scatter(X, X - 3, marker='x')
    f = rendered(p)
    kp.save_plot(p, filename='log.png')
    # source: published - plan 1.1 table: "Non-positive values must not
    # raise at render"
    assert f.axes[0].get_yscale() == 'log'
    assert (in_tmp / 'plots' / 'log.png').exists()
    assert np.isfinite(f.axes[0].get_ylim()).all()


# ===========================================================================
# 4. Legend, grid, twin axis
# ===========================================================================
def test_legend_entries_are_merged_by_label():
    p = fig_()
    p.line(X, Y, legend_label='a')
    p.line(X, Y + 1, legend_label='a', line_dash='dashed')
    p.circle(X, Y + 2, legend_label='b')
    ax = rendered(p).axes[0]
    # source: published - plan 1.1 appearance defaults: "legend top-right
    # with entries merged by label"
    assert sorted(legend_labels(ax)) == ['a', 'b']
    assert len(legend_labels(ax)) == 2
    assert len(all_legends(ax.figure)) == 1


def test_no_legend_without_a_legend_label():
    p = fig_()
    p.line(X, Y)
    p.circle(X, Y)
    # source: published - Bokeh draws a legend only for glyphs with
    # legend_label
    assert rendered(p).axes[0].get_legend() is None


def test_legend_is_top_right_by_default():
    p = fig_()
    p.line(X, Y, legend_label='a')
    # source: published - plan 1.1 appearance defaults: "legend top-right"
    assert legend_loc(rendered(p).axes[0]) == 'upper right'


@pytest.mark.parametrize("bokeh_name, mpl_name", [
    # source: published - plan 1.1 table: Bokeh names (top_left, top_center,
    # bottom_center, ...) mapped to Matplotlib's; these are the ones the
    # primitives use (inventory 2.3) plus the default
    ('top_left', 'upper left'), ('top_center', 'upper center'),
    ('bottom_center', 'lower center'), ('top_right', 'upper right')])
@pytest.mark.parametrize("before", [False, True])
def test_legend_location(bokeh_name, mpl_name, before):
    p = fig_()
    if before:
        p.legend.location = bokeh_name
    p.line(X, Y, legend_label='a')
    if not before:
        p.legend.location = bokeh_name
    assert legend_loc(rendered(p).axes[0]) == mpl_name


def test_xgrid_grid_line_color_none_turns_off_vertical_grid_only():
    p = fig_()
    p.line(X, Y)
    p.xgrid.grid_line_color = None
    ax = rendered(p).axes[0]
    # source: published - plan 1.1 table: "Turns off the vertical grid lines"
    assert not grid_visible(ax, 'x')
    assert grid_visible(ax, 'y')


def test_trap_twin_axis_structure():
    p = solvearcs_like()
    f = rendered(p)
    main, twin = f.axes
    # source: published - plan 0.3 render(): "One Axes per recorded figure
    # ... plus one twiny Axes when extra_x_ranges is used"; plan 1.1 table
    # (extra_x_ranges + add_layout(LinearAxis, 'above') -> twiny())
    assert len(f.axes) == 2
    assert main.get_shared_y_axes().joined(main, twin)
    assert twin.get_xlim() == (10., 60.)
    assert main.get_xlim() == (4000., 4100.)
    assert twin.get_xlabel() == 'CCD Y pix'
    assert twin.xaxis.get_label_position() == 'top'
    assert main.get_xlabel() == 'Wavelength (A)'
    assert title_of(main) == 'ATLAS/ARC LINE FITS'
    # glyphs go to the axes of their x_range_name
    n_main = len(glyph_lines(main)) + len(markers(main))
    n_twin = len(glyph_lines(twin)) + len(markers(twin))
    assert (n_main, n_twin) == (3, 5)


def test_trap_twin_axis_legend_is_one_merged_legend_on_the_top_axes():
    f = rendered(solvearcs_like())
    legends = all_legends(f)
    top = max(enumerate(f.axes), key=lambda ia: (ia[1].get_zorder(), ia[0]))[1]
    # source: published - plan 1.1 trap "Legend across the twin axis": one
    # legend of six entries, 'Arc' and 'Atlas' each used by two glyphs,
    # handles collected from both axes, merged by label, drawn on the axes
    # that is on top
    assert len(legends) == 1
    assert top.get_legend() is legends[0]
    labels = legend_labels(top)
    assert len(labels) == 6
    assert sorted(labels) == sorted(
        ['Arc', 'Atlas', 'Fit', 'Peaks', 'Cat lines', 'Resid'])


def test_trap_twin_axis_draw_order_follows_the_description_not_call_order():
    a = rendered(solvearcs_like('main_first'))
    b = rendered(solvearcs_like('twin_first'))
    # source: published - plan "Core design decision" (order independence)
    # and plan 1.1 trap "Draw order with a twin axis": glyphs are placed on
    # the axes of their range name whatever the call order
    assert len(a.axes) == len(b.axes) == 2
    for axa, axb in zip(a.axes, b.axes):
        assert len(glyph_lines(axa)) == len(glyph_lines(axb))
        assert len(markers(axa)) == len(markers(axb))
        assert axa.get_xlim() == axb.get_xlim()
    # the twin axes is drawn after the base axes: its z-order is not lower
    assert a.axes[1].get_zorder() >= a.axes[0].get_zorder()


def test_trap_twin_axis_draw_order_later_range_is_on_top_in_the_png():
    # Two thick lines crossing the same pixels: a red line on the main range
    # recorded first, a blue one on the 'pix' range recorded second.
    p = fig_(plot_width=400, plot_height=300)
    p.x_range = kp.Range1d(0., 10.)
    p.y_range = kp.Range1d(0., 10.)
    p.extra_x_ranges = {'pix': kp.Range1d(0., 10.)}
    p.add_layout(kp.LinearAxis(x_range_name='pix', axis_label='pix'), 'above')
    p.line([0., 10.], [0., 10.], color='red', line_width=30)
    p.line([0., 10.], [0., 10.], color='blue', line_width=30,
           x_range_name='pix')
    f = kp.render(p)
    buf = io.BytesIO()
    f.savefig(buf, format='png')
    buf.seek(0)
    px = np.asarray(Image.open(buf).convert('RGB'))
    x, y = f.axes[0].transData.transform((5., 5.))
    r, g, b = px[int(round(px.shape[0] - y)), int(round(x))]
    # source: published - plan 1.1 trap "Draw order with a twin axis":
    # glyphs on the 'pix' range, recorded last, are on top of the main
    # range's (Bokeh draws in call order; Matplotlib's twin axes are drawn
    # after the base axes, which agrees here)
    assert b > 200 and r < 60


def test_trap_twin_axis_label_fits_inside_the_fixed_size(in_tmp):
    p = solvearcs_like()
    kp.save_plot(p, filename='twin.png')
    f = kp.render(p)
    draw(f)
    canvas = f.canvas if isinstance(f.canvas, FigureCanvasAgg) \
        else FigureCanvasAgg(f)
    renderer = canvas.get_renderer()
    W, H = f.get_size_inches() * f.dpi
    twin = f.axes[1]
    bb = twin.xaxis.label.get_window_extent(renderer)
    # source: published - plan 1.1 trap "Fixed pixel size": the twin-axis
    # label must fit inside the fixed size
    assert png_size(in_tmp / 'plots' / 'twin.png') == (1000, 600)
    # Tolerance: 1 px for text-extent rounding.
    assert bb.y1 <= H + 1 and bb.y0 >= -1 and bb.x0 >= -1 and bb.x1 <= W + 1


def test_add_layout_only_supports_above():
    p = fig_()
    p.extra_x_ranges = {'pix': kp.Range1d(0., 1.)}
    with pytest.raises(TypeError):
        # source: published - plan 1.1 table: "Only 'above' is supported"
        p.add_layout(kp.LinearAxis(x_range_name='pix', axis_label='a'),
                     'below')


# ===========================================================================
# 5. gridplot and render
# ===========================================================================
def _four(titles=('A', 'B', 'C', 'D'), **kw):
    figs = []
    for i, t in enumerate(titles):
        p = kp.figure(title=t, x_axis_label='x', y_axis_label='y',
                      plot_width=999, plot_height=999)
        p.line(X, Y + i, legend_label='L%d' % i)
        figs.append(p)
    return figs


def test_gridplot_layout_is_row_major():
    figs = _four()
    g = kp.gridplot(figs, ncols=2, plot_width=500, plot_height=300,
                    toolbar_location=None)
    f = rendered(g)
    # source: published - plan 0.3 render(): "One Axes per recorded figure
    # (row-major for a grid)"
    assert [title_of(a) for a in f.axes] == ['A', 'B', 'C', 'D']
    pos = [a.get_position() for a in f.axes]
    cx = [(p.x0 + p.x1) / 2 for p in pos]
    cy = [(p.y0 + p.y1) / 2 for p in pos]
    # figure coordinates have y increasing upwards: row 0 has the larger y
    assert cx[0] < cx[1] and cx[2] < cx[3]
    assert cy[0] > cy[2] and cy[1] > cy[3]
    # Tolerance: row/column centres agree within a tenth of a cell width;
    # constrained layout shifts cells by label widths, not by cell sizes.
    assert abs(cy[0] - cy[1]) < 0.05 and abs(cx[0] - cx[2]) < 0.1


def test_gridplot_each_subplot_has_its_own_legend():
    f = rendered(kp.gridplot(_four(), ncols=2, plot_width=500,
                             plot_height=300, toolbar_location=None))
    # source: published - plan 0.3 render(): a legend is drawn per figure
    assert [legend_labels(a) for a in f.axes] == [['L0'], ['L1'], ['L2'],
                                                  ['L3']]


@pytest.mark.parametrize("n, ncols, w, h, size", [
    # source: published - plan 1.1 verification: grid is ncols*plot_width
    # by nrows*plot_height; plan 1.2: "2 amps -> 500x600, 4 amps ncols=2 ->
    # 1000x600" (SubtractOverscan passes plot_width=500, plot_height=300)
    (1, 1, 500, 300, (500, 300)),
    (2, 1, 500, 300, (500, 600)),
    (4, 2, 500, 300, (1000, 600)),
    # source: analytic - 3 figures, 2 columns -> 2 rows
    (3, 2, 500, 300, (1000, 600)),
    (4, 2, 250, 200, (500, 400))])
def test_exact_png_size_of_grids(in_tmp, n, ncols, w, h, size):
    figs = _four(titles=list('ABCD')[:n])
    g = kp.gridplot(figs, ncols=ncols, plot_width=w, plot_height=h,
                    toolbar_location=None)
    assert kp.pixel_size(g) == size
    kp.save_plot(g, filename='grid.png')
    assert png_size(in_tmp / 'plots' / 'grid.png') == size
    assert len(rendered(g).axes) == n


def test_gridplot_size_overrides_each_figures_own_size():
    figs = _four()
    g = kp.gridplot(figs, ncols=2, plot_width=300, plot_height=200)
    # source: published - plan 1.1 table: "plot_width and plot_height are
    # per cell and override each figure's own size" (figures are 999x999)
    assert kp.pixel_size(g) == (600, 400)
    f = kp.render(g)
    assert tuple(f.get_size_inches() * f.dpi) == pytest.approx((600, 400))


def test_same_figures_in_two_gridplots():
    figs = _four()
    g1 = kp.gridplot(figs, ncols=2, plot_width=500, plot_height=300,
                     toolbar_location=None)
    g2 = kp.gridplot(figs, ncols=1, plot_width=200, plot_height=100,
                     toolbar_location=None)
    f1, f2 = rendered(g1), rendered(g2)
    # source: published - plan 1.1 table: "The same figure list may be passed
    # to several gridplot calls" (SubtractOverscan)
    assert kp.pixel_size(g1) == (1000, 600)
    assert kp.pixel_size(g2) == (200, 400)
    assert len(f1.axes) == len(f2.axes) == 4
    assert describe(rendered(g1)) == describe(f1)   # re-rendering is stable


def test_a_figure_can_be_rendered_many_times():
    p = solvearcs_like()
    # source: published - plan "Core design decision", "Reuse": recorded
    # descriptions can be rendered any number of times
    assert describe(rendered(p)) == describe(rendered(p))


def test_render_into_an_existing_figure_clears_and_resizes_it():
    from matplotlib.figure import Figure
    target = Figure()
    first = rendered(solvearcs_like(), fig=target)
    # source: published - plan 0.3: render(obj, fig=None): "with fig given,
    # that figure is cleared, resized to the object's pixel size and drawn
    # into"
    assert first is target
    assert len(target.axes) == 2
    assert tuple(target.get_size_inches() * target.dpi) == pytest.approx(
        (1000, 600))
    g = kp.gridplot(_four(), ncols=2, plot_width=500, plot_height=300)
    again = rendered(g, fig=target)
    assert again is target
    assert len(target.axes) == 4
    assert tuple(target.get_size_inches() * target.dpi) == pytest.approx(
        (1000, 600))
    assert not target.legends
    p = fig_(plot_width=250, plot_height=125)
    p.line(X, Y)
    rendered(p, fig=target)
    assert len(target.axes) == 1 and target.axes[0].get_legend() is None
    assert tuple(target.get_size_inches() * target.dpi) == pytest.approx(
        (250, 125))


def test_pixel_size_and_plot_label():
    p = fig_(title='Label', plot_width=640, plot_height=480)
    q = fig_(title='Other')
    g = kp.gridplot([p, q], ncols=2, plot_width=100, plot_height=50)
    # source: published - plan 0.3 interfaces: pixel_size -> (width, height);
    # plot_label -> title text (a grid: its first figure's title)
    assert kp.pixel_size(p) == (640, 480)
    assert kp.plot_label(p) == 'Label'
    assert kp.plot_label(g) == 'Label'


# ===========================================================================
# 6. NaN handling
# ===========================================================================
def _runs(lines):
    out = []
    for l in lines:
        li = line_info(l)
        ok = np.isfinite(li['x']) & np.isfinite(li['y'])
        run = []
        for xv, good in zip(li['x'], ok):
            if good:
                run.append(xv)
            elif run:
                out.append(run)
                run = []
        if run:
            out.append(run)
    return sorted(out)


def test_trap_nan_leaves_a_gap_in_a_line():
    x = np.arange(7.)
    y = np.array([0., 1., 2., np.nan, 4., 5., 6.])
    p = fig_()
    p.line(x, y)
    # source: published - plan 1.1 trap "NaN handling": Bokeh leaves a gap
    # in a line at NaN; the finite points form two separate runs
    assert _runs(glyph_lines(rendered(p).axes[0])) == [[0., 1., 2.],
                                                       [4., 5., 6.]]


def test_trap_nan_markers_are_dropped():
    x = np.arange(5.)
    y = np.array([1., np.nan, 3., np.nan, 5.])
    p = fig_()
    p.circle(x, y)
    p.scatter(x, y, marker='x')
    ms = markers(rendered(p).axes[0])
    # source: published - plan 1.1 trap "NaN handling": NaN markers dropped
    for m in ms:
        ok = np.isfinite(m['x']) & np.isfinite(m['y'])
        assert list(m['x'][ok]) == [0., 2., 4.]
        assert list(m['y'][ok]) == [1., 3., 5.]


def test_trap_all_nan_series_renders_an_empty_plot(in_tmp):
    x = np.arange(5.)
    y = np.full(5, np.nan)
    p = fig_()
    p.line(x, y, legend_label='nothing')
    p.circle(x, y)
    p.diamond(x, y)
    f = rendered(p)
    kp.save_plot(p, filename='nan.png')
    # source: published - plan 1.1 trap "NaN handling": "an all-NaN series
    # renders an empty plot and does not raise"
    assert len(f.axes) == 1
    assert np.isfinite(f.axes[0].get_xlim()).all()
    assert np.isfinite(f.axes[0].get_ylim()).all()
    assert (in_tmp / 'plots' / 'nan.png').exists()
    assert png_size(in_tmp / 'plots' / 'nan.png') == (400, 300)


# ===========================================================================
# 7. Fixed pixel size
# ===========================================================================
@pytest.mark.parametrize("kind", ['plain', 'log', 'image', 'twin', 'legend'])
def test_exact_png_size_of_a_single_figure(in_tmp, kind):
    if kind == 'twin':
        p = solvearcs_like()
    else:
        p = kp.figure(title='T', x_axis_label='x', y_axis_label='y',
                      plot_width=1000, plot_height=600,
                      **({'y_axis_type': 'log'} if kind == 'log' else {}))
        if kind == 'image':
            p.image([np.arange(20.).reshape(4, 5)], x=0, y=0, dw=5, dh=4,
                    palette="Spectral11")
        else:
            p.line(X + 1, X + 1, legend_label='a')
        if kind == 'legend':
            p.legend.location = 'bottom_center'
    kp.save_plot(p, filename='one.png')
    # source: published - plan 1.1 verification and plan 1.2: a saved PNG of
    # a single figure is exactly plot_width x plot_height pixels (1000x600)
    assert png_size(in_tmp / 'plots' / 'one.png') == (1000, 600)


def test_trap_fixed_pixel_size_with_a_very_long_title(in_tmp):
    title = ('Img # 12345 (BH2) Slicer: Medium, Grating: BM, CWave: 4500, '
             'Filter: KBlue  ATLAS/ARC LINE FITS for a very long title that '
             'would overflow a 1000 pixel wide figure if it were not wrapped '
             'or shrunk, xxxxxxxx yyyyyyyy zzzzzzzz 0123456789 0123456789')
    p = kp.figure(title=title, x_axis_label='x', y_axis_label='y',
                  plot_width=1000, plot_height=600)
    p.line(X, Y)
    kp.save_plot(p, filename='long.png')
    f = kp.render(p)
    draw(f)
    canvas = f.canvas if isinstance(f.canvas, FigureCanvasAgg) \
        else FigureCanvasAgg(f)
    bb = f.axes[0]._left_title.get_window_extent(canvas.get_renderer())
    # source: published - plan 1.1 trap "Fixed pixel size": save with the
    # figure's own size and never a tight bounding box; long titles must fit
    # inside the fixed size (constrained layout)
    assert png_size(in_tmp / 'plots' / 'long.png') == (1000, 600)
    # Tolerance: 1 px for text-extent rounding.
    assert bb.x0 >= -1 and bb.x1 <= 1000 + 1
    assert bb.y0 >= -1 and bb.y1 <= 600 + 1


def test_long_title_in_a_narrow_grid_cell_keeps_the_grid_size(in_tmp):
    figs = []
    for i in range(4):
        p = kp.figure(title='Img # 00123 OSCAN [10:20, 30:40] amp %d, noise: '
                      '1.234 e-/px  ' % i + 'x' * 40, x_axis_label='x',
                      y_axis_label='y', plot_width=500, plot_height=300)
        p.line(X, Y)
        figs.append(p)
    g = kp.gridplot(figs, ncols=2, plot_width=500, plot_height=300,
                    toolbar_location=None)
    kp.save_plot(g, filename='cell.png')
    # source: published - plan 1.1 verification: exact grid size
    assert png_size(in_tmp / 'plots' / 'cell.png') == (1000, 600)


# ===========================================================================
# 8. save_plot
# ===========================================================================
def test_save_plot_writes_plots_dir_and_logs(in_tmp, kcwi_log):
    p = fig_(plot_width=300, plot_height=200)
    p.line(X, Y)
    kp.save_plot(p, filename='a.png')
    # source: published - plan 1.1: save_plot writes plots/<filename>
    # relative to the current directory, creates plots/ if missing, logs the
    # same ">>> Saving to" line
    assert (in_tmp / 'plots' / 'a.png').is_file()
    assert any('>>> Saving to' in r.getMessage() and 'a.png' in r.getMessage()
               for r in kcwi_log.records if r.name == 'KCWI')
    assert png_size(in_tmp / 'plots' / 'a.png') == (300, 200)


def test_save_plot_default_file_name(in_tmp):
    p = fig_()
    p.line(X, Y)
    kp.save_plot(p)
    # source: published - inventory 2.5: default file plots/kcwi_drp_plot.png
    assert (in_tmp / 'plots' / 'kcwi_drp_plot.png').is_file()


def test_save_plot_uses_an_existing_plots_dir(in_tmp):
    (in_tmp / 'plots').mkdir()
    p = fig_()
    p.line(X, Y)
    kp.save_plot(p, filename='b.png')
    assert (in_tmp / 'plots' / 'b.png').is_file()


def test_plot_level_0_save_without_any_viewer(in_tmp):
    p = fig_(plot_width=500, plot_height=300)
    p.line(X, Y, legend_label='a')
    p.xgrid.grid_line_color = None
    kp.oplot_slices(p, [0., 100.])
    kp.set_plot_lims(p, xlim=[0., 115.], ylim=[0., 100.])
    # no StartPlotViewer event, no session: SolveArcs saves at every level
    kp.save_plot(p, filename='arc_00001_nlines.png')
    # source: published - plan 1.2 "Plot level 0": with no viewer started,
    # save_plot on a figure writes its PNG
    assert png_size(in_tmp / 'plots' / 'arc_00001_nlines.png') == (500, 300)


def test_mutating_inputs_after_the_call_does_not_change_the_png(in_tmp):
    x = np.linspace(0., 10., 50)
    y = np.sin(x)
    p = fig_()
    p.line(x, y)
    p.circle(x, y)
    p.quad(top=np.array([1.]), bottom=0, left=np.array([0.]),
           right=np.array([1.]))
    x[:] = 0.
    y[:] = 99.
    kp.save_plot(p, filename='mutated.png')
    ref_x = np.linspace(0., 10., 50)
    q = fig_()
    q.line(ref_x, np.sin(ref_x))
    q.circle(ref_x, np.sin(ref_x))
    q.quad(top=np.array([1.]), bottom=0, left=np.array([0.]),
           right=np.array([1.]))
    kp.save_plot(q, filename='ref.png')
    # source: published - plan 1.1 verification: "Mutating an input array in
    # place after the glyph call does not change the saved PNG"
    assert np.array_equal(png_pixels(in_tmp / 'plots' / 'mutated.png'),
                          png_pixels(in_tmp / 'plots' / 'ref.png'))


def test_save_plot_runs_without_display_or_pyplot_in_a_fresh_process(
        tmp_path):
    out = _run_probe(tmp_path, strict_env='1')
    # source: published - plan 1.1 verification: rendering runs with
    # DISPLAY unset; wrapper rule 4: no pyplot
    assert out['saved_size'] == [300, 200]
    assert out['after_render'] == {'matplotlib.pyplot': False,
                                   'bokeh': False, 'selenium': False}


# ===========================================================================
# 9. Snapshot
# ===========================================================================
def test_snapshot_rendered_data_is_the_original_after_in_place_overwrite():
    x = np.arange(6.)
    y = np.arange(6.) * 2.
    xs = [0., 1., 2.]
    p = fig_()
    p.line(x, y)
    p.circle(xs, xs)
    x[:] = -1.
    y[:] = -1.
    xs[0] = 99.     # a plain list mutated in place
    ax = rendered(p).axes[0]
    # source: published - plan 1.1 rule 3 and plan 1.2 "Snapshot test":
    # record a line from an array, overwrite the array in place, render: the
    # rendered data is the original
    li = line_info(glyph_lines(ax)[0])
    assert list(li['x']) == [0., 1., 2., 3., 4., 5.]
    assert list(li['y']) == [0., 2., 4., 6., 8., 10.]
    assert list(markers(ax)[0]['x']) == [0., 1., 2.]


def test_snapshot_of_an_image_quad_and_y_range_list():
    a = np.arange(6.).reshape(2, 3)
    top = np.array([1., 2.])
    yr = [2., 8.]
    p = fig_(y_range=yr)
    p.image([a], x=0, y=0, dw=3, dh=2, palette="Spectral11")
    p.quad(top=top, bottom=0, left=[0., 1.], right=[1., 2.])
    a[:] = 0.
    top[:] = 50.
    yr[:] = [0., 1.]
    ax = rendered(p).axes[0]
    # source: published - plan 1.1 rule 3: arrays, scalars and lists are
    # copied when recorded
    assert np.array_equal(np.asarray(ax.images[0].get_array()),
                          np.arange(6.).reshape(2, 3))
    assert [r['box'][3] for r in rects(ax)] == [1., 2.]
    assert ax.get_ylim() == (2., 8.)


class StubViewer:
    """The interface of PlotViewer that bokeh_plot/bokeh_clear use."""

    def __init__(self):
        self.shown = []
        self.snaps = []
        self.cleared = 0

    def show(self, obj):
        self.shown.append(obj)
        self.snaps.append(kp.snapshot(obj))

    def clear(self):
        self.cleared += 1


def test_snapshot_handed_to_a_stub_viewer_is_not_changed_by_later_calls():
    from kcwidrp.core.bokeh_plotting import bokeh_plot
    y = np.arange(5.)
    p = fig_(title='first')
    p.line(np.arange(5.), y, legend_label='one')
    viewer = StubViewer()
    bokeh_plot(p, viewer)
    # the primitive carries on: more glyphs, a new title, array overwritten
    p.line(np.arange(5.), y + 10, color='red', legend_label='two')
    p.title.text = 'second'
    p.y_range.start = -100
    y[:] = 0.
    # source: published - plan 1.2 "Snapshot test": "the same for a figure
    # submitted to a stub viewer and then modified by a further glyph call on
    # the primitive's side"; plan 0.3: bokeh_plot(plot, session) calls
    # session.show(plot)
    assert viewer.shown == [p] and viewer.shown[0] is p
    f = rendered(viewer.snaps[0])
    ax = f.axes[0]
    assert title_of(ax) == 'first'
    assert len(glyph_lines(ax)) == 1
    assert list(line_info(glyph_lines(ax)[0])['y']) == [0., 1., 2., 3., 4.]
    assert legend_labels(ax) == ['one']
    assert ax.get_ylim()[0] > -1
    # and the primitive's own object did get the changes
    g = rendered(p).axes[0]
    assert title_of(g) == 'second' and len(glyph_lines(g)) == 2


def test_snapshot_is_independent_for_a_grid():
    figs = _four()
    g = kp.gridplot(figs, ncols=2, plot_width=500, plot_height=300)
    snap = kp.snapshot(g)
    figs[0].line(X, Y + 50)
    figs[0].title.text = 'changed'
    # source: published - plan 0.3: snapshot(obj) -> an independent deep copy
    assert snap is not g
    f = rendered(snap)
    assert title_of(f.axes[0]) == 'A'
    assert len(glyph_lines(f.axes[0])) == 1


# ===========================================================================
# 10. Strict mode (STRICT true: TypeError at the call)
# ===========================================================================
def test_strict_unsupported_method_raises_type_error():
    p = fig_()
    with pytest.raises(TypeError):
        # source: published - plan 1.1 rule 1 / verification: with STRICT
        # true an unsupported method raises TypeError at the call
        p.wedge(X, Y)


@pytest.mark.parametrize("glyph", ['line', 'circle', 'scatter', 'diamond'])
def test_strict_unsupported_keyword_raises_type_error(glyph):
    p = fig_()
    with pytest.raises(TypeError):
        # source: published - plan 1.1 rule 1: unsupported keyword
        getattr(p, glyph)(X, Y, bogus=1)
    if glyph != 'scatter':
        with pytest.raises(TypeError):
            # source: published - plan 1.1 table: `marker` belongs to scatter
            getattr(p, glyph)(X, Y, marker='x')


def test_strict_unsupported_figure_keyword_raises():
    with pytest.raises(TypeError):
        # source: published - plan 1.2 verification mutations:
        # figure(x_axis_type='log') must be rejected
        kp.figure(x_axis_type='log')


@pytest.mark.parametrize("stmt", [
    "p.ygrid.grid_line_color = None",
    "p.xgrid.grid_line_color = 'red'",
    "p.nonexistent_attribute = 3",
    "p.x_range.range_padding = 0.1",
    "p.legend.click_policy = 'hide'"])
def test_strict_unsupported_attribute_raises_type_error(stmt):
    p = fig_()
    p.line(X, Y, legend_label='a')
    with pytest.raises(TypeError):
        # source: published - plan 1.1 rule 1 and table (xgrid supports only
        # None; ygrid, other attributes and range_padding != 0 are not in
        # the table)
        exec(stmt, {'p': p})


@pytest.mark.parametrize("call", [
    "p.line(X, Y, line_dash='longdash')",
    "p.scatter(X, Y, marker='triangle')",
    "p.image([a], x=0, y=0, dw=1, dh=1, palette='Greys9')",
    "p.add_layout(LinearAxis(x_range_name='pix', axis_label='a'), 'below')",
    "setattr(p.legend, 'location', 'center_left')"])
def test_strict_unsupported_literal_value_raises_type_error(call):
    p = fig_()
    p.extra_x_ranges = {'pix': kp.Range1d(0., 1.)}
    env = {'p': p, 'X': X, 'Y': Y, 'a': np.ones((2, 2)),
           'LinearAxis': kp.LinearAxis}
    with pytest.raises(TypeError):
        # source: published - plan 1.1 table: marker 'x' only, dash names
        # 'dashed'/'dashdot'/'dotted' only, palette "Spectral11" only,
        # add_layout 'above' only, legend locations are Bokeh's names
        exec(call, env)


@pytest.mark.parametrize("value", ['linear', 'datetime', 'mercator'])
def test_strict_y_axis_type_other_than_log_raises(value):
    with pytest.raises(TypeError):
        # source: published - plan 1.1 table: "Any other value is unsupported"
        kp.figure(y_axis_type=value)


def test_strict_column_name_form_is_unsupported():
    p = fig_()
    with pytest.raises(TypeError):
        # source: published - plan 1.1 table: p.line('x', 'y', source=...) is
        # not supported by the wrapper
        p.line('x', 'y', line_color='black', source={'x': X, 'y': Y})


def test_strict_other_constructors_reject_unknown_keywords():
    with pytest.raises(TypeError):
        kp.Range1d(start=1, stop=2)
    with pytest.raises(TypeError):
        kp.LinearAxis(x_range_name='a', bogus=1)
    with pytest.raises(TypeError):
        kp.gridplot([fig_()], ncols=1, plot_width=1, plot_height=1,
                    bogus=1)


# ===========================================================================
# 11. Tolerant mode (STRICT false: one warning, nothing raised, usable)
# ===========================================================================
def _one_warning_naming(kcwi_log, lineno, name):
    ws = warnings_of(kcwi_log)
    # source: published - plan 1.1 rule 1 / plan 1.2: "one warning is logged
    # with the caller's file and line"
    assert len(ws) == 1
    msg = ws[0].getMessage()
    assert name in msg
    assert re.search(re.escape(THIS_FILE) + r'\D{1,8}%d\b' % lineno, msg), msg


def _still_usable(p, in_tmp):
    p.line(X, Y + 1, color='red', legend_label='more')
    f = rendered(p)
    kp.save_plot(p, filename='after.png')
    assert (in_tmp / 'plots' / 'after.png').is_file()
    assert any(mcolors.to_hex(line_info(l)['colour']) == '#ff0000'
               for l in glyph_lines(f.axes[0]))


def test_tolerant_unsupported_method(tolerant, kcwi_log, in_tmp):
    p = fig_()
    p.line(X, Y)
    lineno = sys._getframe().f_lineno + 1
    p.wedge(X, Y)
    _one_warning_naming(kcwi_log, lineno, 'wedge')
    _still_usable(p, in_tmp)


def test_tolerant_unsupported_keyword(tolerant, kcwi_log, in_tmp):
    p = fig_()
    lineno = sys._getframe().f_lineno + 1
    p.line(X, Y, bogus=1, color='green')
    _one_warning_naming(kcwi_log, lineno, 'bogus')
    # source: published - plan 1.1 rule 1: "skip that one call or keyword" --
    # the call itself is kept, without the keyword
    lines = glyph_lines(rendered(p).axes[0])
    assert len(lines) == 1
    assert mcolors.to_hex(line_info(lines[0])['colour']) == '#008000'
    kcwi_log.clear()
    _still_usable(p, in_tmp)


def test_tolerant_unsupported_attribute(tolerant, kcwi_log, in_tmp):
    p = fig_()
    p.line(X, Y)
    lineno = sys._getframe().f_lineno + 1
    p.ygrid.grid_line_color = None
    _one_warning_naming(kcwi_log, lineno, 'ygrid')
    # the grid was not turned off: the assignment was skipped
    assert grid_visible(rendered(p).axes[0], 'y')
    kcwi_log.clear()
    _still_usable(p, in_tmp)


def test_tolerant_unsupported_literal_value(tolerant, kcwi_log, in_tmp):
    p = fig_()
    lineno = sys._getframe().f_lineno + 1
    p.line(X, Y, line_dash='longdash')
    _one_warning_naming(kcwi_log, lineno, 'longdash')
    kcwi_log.clear()
    _still_usable(p, in_tmp)


def test_tolerant_unsupported_figure_keyword(tolerant, kcwi_log, in_tmp):
    lineno = sys._getframe().f_lineno + 1
    p = kp.figure(title='t', x_axis_type='log', plot_width=200,
                  plot_height=100)
    _one_warning_naming(kcwi_log, lineno, 'x_axis_type')
    # source: published - plan 1.1 rule 2: figure() still returns a usable
    # recorder object on failure
    assert isinstance(p, kp.PlotFigure)
    kcwi_log.clear()
    p.line(X, Y)
    kp.save_plot(p, filename='f.png')
    assert png_size(in_tmp / 'plots' / 'f.png') == (200, 100)
    assert warnings_of(kcwi_log) == []


def test_tolerant_unsupported_literal_in_legend_location(
        tolerant, kcwi_log, in_tmp):
    p = fig_()
    p.line(X, Y, legend_label='a')
    lineno = sys._getframe().f_lineno + 1
    p.legend.location = 'center_left'
    _one_warning_naming(kcwi_log, lineno, 'center_left')
    _still_usable(p, in_tmp)


def test_tolerant_column_name_form(tolerant, kcwi_log, in_tmp):
    p = fig_()
    lineno = sys._getframe().f_lineno + 1
    p.line('x', 'y', line_color='black', source={'x': X, 'y': Y})
    ws = warnings_of(kcwi_log)
    # source: published - plan 1.1 table / rule 1: unsupported usage is
    # logged, never raised, and the object stays usable.  One call can break
    # more than one rule (string x/y and source=), so only "at least one
    # warning, all naming this line" is asserted.
    assert ws
    assert all(re.search(re.escape(THIS_FILE) + r'\D{1,8}%d\b' % lineno,
                         w.getMessage()) for w in ws)
    kcwi_log.clear()
    _still_usable(p, in_tmp)


def test_tolerant_gridplot_with_garbage_returns_a_grid(
        tolerant, kcwi_log, in_tmp):
    p = fig_()
    p.line(X, Y)
    g = kp.gridplot([p, 'not a figure'], ncols=2, plot_width=100,
                    plot_height=100)
    # source: published - plan 1.1 rule 2: gridplot catches its own
    # exceptions and still returns a usable recorder object
    assert isinstance(g, kp.PlotGrid)
    kp.save_plot(g, filename='g.png')    # does not raise


def test_tolerant_render_failure_inside_matplotlib_is_logged(
        tolerant, kcwi_log, in_tmp, monkeypatch):
    p = fig_()
    p.line(X, Y)

    def boom(*a, **k):
        raise RuntimeError('canvas boom')

    from matplotlib.figure import Figure
    monkeypatch.setattr(FigureCanvasAgg, 'draw', boom)
    monkeypatch.setattr(FigureCanvasAgg, 'print_png', boom)
    monkeypatch.setattr(Figure, 'draw', boom)
    monkeypatch.setattr(Figure, 'savefig', boom)
    # source: published - plan 1.1 verification: "a render that fails inside
    # Matplotlib (forced in a test by patching the canvas to raise) logs a
    # warning and save_plot returns normally"
    assert kp.save_plot(p, filename='boom.png') is None
    assert warnings_of(kcwi_log)


# ===========================================================================
# 12. Failure isolation
# ===========================================================================
def test_save_plot_to_an_unwritable_directory_does_not_raise(
        tolerant, kcwi_log, in_tmp):
    (in_tmp / 'plots').write_text('a file where the directory should be')
    p = fig_()
    p.line(X, Y)
    # source: published - plan 1.2 "Failure isolation": save_plot with an
    # unwritable directory logs a warning and returns; it does not raise
    assert kp.save_plot(p, filename='x.png') is None
    assert warnings_of(kcwi_log)
    assert (in_tmp / 'plots').is_file()


def test_bokeh_plot_with_a_viewer_whose_show_raises(
        tolerant, kcwi_log):
    from kcwidrp.core.bokeh_plotting import bokeh_plot

    class Broken:
        def show(self, obj):
            raise RuntimeError('viewer is broken')

    p = fig_()
    p.line(X, Y)
    # source: published - plan 1.2 "Failure isolation": the same for
    # bokeh_plot with a viewer whose show raises (plan 0.3: bokeh_plot calls
    # session.show(plot); wrapper rule 2: never raises into a primitive)
    assert bokeh_plot(p, Broken()) is None
    assert warnings_of(kcwi_log)


def test_bokeh_plot_calls_show_and_ignores_none():
    from kcwidrp.core.bokeh_plotting import bokeh_plot
    viewer = StubViewer()
    p = fig_()
    bokeh_plot(p, viewer)
    # source: published - plan 0.3: bokeh_plot(plot, session) ->
    # session.show(plot); no-op for None
    assert viewer.shown == [p]
    assert bokeh_plot(p, None) is None


def test_bokeh_clear_with_none_does_not_raise():
    from kcwidrp.core.bokeh_plotting import bokeh_clear
    # source: published - plan 1.2: bokeh_clear(None) must not raise (today's
    # bokeh_clear would raise on None; the new one must not).  STRICT true.
    assert bokeh_clear(None) is None


def test_bokeh_clear_calls_clear_and_survives_a_broken_viewer(
        tolerant, kcwi_log):
    from kcwidrp.core.bokeh_plotting import bokeh_clear
    viewer = StubViewer()
    bokeh_clear(viewer)
    # source: published - plan 0.3: bokeh_clear(session) -> session.clear()
    assert viewer.cleared == 1

    class Broken:
        def clear(self):
            raise RuntimeError('boom')

    assert bokeh_clear(Broken()) is None
    assert warnings_of(kcwi_log)


# ===========================================================================
# 13. Threads
# ===========================================================================
def _thread_plot(i, j):
    x = np.linspace(0., 10., 120)
    p = kp.figure(title='thread %d plot %d' % (i, j), x_axis_label='x',
                  y_axis_label='y', plot_width=300, plot_height=200)
    colour = ['red', 'green', 'blue', 'black'][i]
    p.line(x, np.sin(x * (j + 1)) + i, color=colour, legend_label='s')
    p.circle(x[::10], np.cos(x[::10]), size=3 + j % 3, fill_color=colour)
    return p


def test_four_threads_render_and_save_ten_plots_each_concurrently(in_tmp):
    errors = []
    barrier = threading.Barrier(4)

    def work(i):
        try:
            barrier.wait(timeout=30)
            for j in range(10):
                p = _thread_plot(i, j)
                f = kp.render(p)
                if len(f.axes) != 1:
                    raise AssertionError('axes %d' % len(f.axes))
                kp.save_plot(p, filename='t%d_%d.png' % (i, j))
        except BaseException as exc:    # collected, re-raised below
            errors.append((i, repr(exc)))

    threads = [threading.Thread(target=work, args=(i,)) for i in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=120)
    # source: published - plan 1.2 "Thread test": four threads each render
    # and save ten plots concurrently; all files exist with the right pixel
    # size
    assert errors == []
    for i in range(4):
        for j in range(10):
            path = in_tmp / 'plots' / ('t%d_%d.png' % (i, j))
            assert path.is_file(), path
            assert png_size(path) == (300, 200)
    # and nothing was corrupted by the concurrency: one plot per thread is
    # identical, pixel for pixel, to the same plot saved on its own
    for i in range(4):
        kp.save_plot(_thread_plot(i, 7), filename='solo%d.png' % i)
        assert np.array_equal(
            png_pixels(in_tmp / 'plots' / ('t%d_7.png' % i)),
            png_pixels(in_tmp / 'plots' / ('solo%d.png' % i)))


# ===========================================================================
# 14. Process-level: imports, STRICT from the environment, headless render
# ===========================================================================
PROBE = r'''
import json, os, sys
MODS = ("matplotlib.pyplot", "bokeh", "selenium")

def seen():
    return {m: any(k == m or k.startswith(m + ".") for k in sys.modules)
            for m in MODS}

import kcwidrp.core.kcwi_plotting as kp
after_import = seen()
p = kp.figure(title="t", x_axis_label="x", y_axis_label="y",
              plot_width=300, plot_height=200)
p.line([0, 1, 2], [0, 1, 4], legend_label="a")
kp.save_plot(p, filename="probe.png")
from PIL import Image
with Image.open(os.path.join("plots", "probe.png")) as im:
    size = list(im.size)
print("PROBE" + json.dumps(dict(strict=kp.STRICT, after_import=after_import,
                                after_render=seen(), saved_size=size,
                                display=os.environ.get("DISPLAY"))))
'''


def _run_probe(tmp_path, strict_env):
    env = {k: v for k, v in os.environ.items()
           if k not in ('DISPLAY', 'MPLBACKEND', 'KCWI_PLOT_STRICT')}
    env['PYTHONPATH'] = str(REPO)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    if strict_env is not None:
        env['KCWI_PLOT_STRICT'] = strict_env
    res = subprocess.run([sys.executable, '-c', PROBE], cwd=str(tmp_path),
                         env=env, capture_output=True, text=True,
                         timeout=120)
    assert res.returncode == 0, res.stderr
    line = [l for l in res.stdout.splitlines() if l.startswith('PROBE')][-1]
    out = json.loads(line[len('PROBE'):])
    assert out['display'] is None
    return out


def test_import_leaves_pyplot_bokeh_and_selenium_out_of_sys_modules(tmp_path):
    out = _run_probe(tmp_path, strict_env=None)
    # source: published - plan 1.1 verification: importing
    # kcwidrp.core.kcwi_plotting in a fresh interpreter leaves
    # matplotlib.pyplot, bokeh and selenium out of sys.modules
    assert out['after_import'] == {'matplotlib.pyplot': False,
                                   'bokeh': False, 'selenium': False}


def test_strict_flag_comes_from_the_environment(tmp_path):
    on = _run_probe(tmp_path, strict_env='1')
    off = _run_probe(tmp_path, strict_env=None)
    # source: published - plan 0.3: STRICT is "true when KCWI_PLOT_STRICT=1
    # at import"; plan 1.1 rule 1: the default in a reduction is False
    assert on['strict'] is True
    assert off['strict'] is False
