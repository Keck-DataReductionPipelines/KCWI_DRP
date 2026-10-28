"""Bokeh-shaped plotting wrapper drawn with Matplotlib (bokeh-plan 1.1).

Primitives build plots with the same calls they used with Bokeh
(``figure``, ``p.line``, ``p.x_range.start = ...``, ``gridplot`` ...).  The
objects returned here only *record* those calls; a Matplotlib ``Figure`` is
built from the record when :func:`render` is called, by :func:`save_plot` or
by the live viewer.  Every array is copied when it is recorded, so a primitive
may go on changing its arrays, and a recorded plot can be rendered any number
of times, on any thread.

Exactly the Bokeh features listed in ``docs/plot_inventory.md`` section 2 are
supported.  Anything else is *unsupported*:

* ``STRICT`` true (the test suite, or ``KCWI_PLOT_STRICT=1`` in the
  environment at import): ``TypeError`` (:class:`UnsupportedPlotFeature`) is
  raised at the call;
* ``STRICT`` false (a reduction): one warning naming the calling file and
  line is logged, that call, keyword or assignment is skipped, and the plot
  carries on.

With ``STRICT`` false no plotting call raises into a primitive.

Only ``matplotlib.figure.Figure`` and ``FigureCanvasAgg`` are used; importing
this module imports neither ``matplotlib.pyplot`` nor Bokeh nor Selenium.
"""
import copy
import logging
import math
import numbers
import os
import sys

import numpy as np

logger = logging.getLogger('KCWI')

STRICT = os.environ.get('KCWI_PLOT_STRICT', '') == '1'

__all__ = ['STRICT', 'UnsupportedPlotFeature', 'figure', 'PlotFigure',
           'Range1d', 'LinearAxis', 'gridplot', 'PlotGrid', 'render',
           'pixel_size', 'snapshot', 'plot_label', 'save_plot',
           'get_plot_lims', 'set_plot_lims', 'oplot_slices',
           'configure_plot_driver']


class UnsupportedPlotFeature(TypeError):
    """Raised, when ``STRICT`` is true, for a Bokeh feature the wrapper
    does not implement."""


# ---------------------------------------------------------------------------
# Strict / tolerant handling
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
_OWN_FILES = {os.path.join(_HERE, 'kcwi_plotting.py'),
              os.path.join(_HERE, '_plot_render.py')}


def _caller():
    """'file:line' of the first frame outside this wrapper."""
    try:
        f = sys._getframe(1)
        while f is not None:
            fn = os.path.abspath(f.f_code.co_filename)
            if fn not in _OWN_FILES:
                return "%s:%d" % (f.f_code.co_filename, f.f_lineno)
            f = f.f_back
    except Exception:
        pass
    return "<unknown>"


def _unsupported(what):
    """Report an unsupported feature: raise when STRICT, else warn once."""
    where = _caller()
    if STRICT:
        raise UnsupportedPlotFeature("%s: unsupported plotting feature: %s"
                                     % (where, what))
    logger.warning("%s: unsupported plotting feature: %s (ignored)"
                   % (where, what))


def _failed(what, exc):
    """An unexpected error inside the wrapper: re-raise when STRICT,
    otherwise log a warning (traceback at debug level) and carry on."""
    if STRICT or isinstance(exc, UnsupportedPlotFeature):
        raise exc
    logger.warning("%s: plotting failed in %s: %s: %s (ignored)"
                   % (_caller(), what, type(exc).__name__, exc))
    logger.debug("plotting traceback", exc_info=exc)


class _Null:
    """Returned for an unsupported attribute when STRICT is false: accepts
    any further call, attribute read or assignment and does nothing, so the
    primitive's following lines do not fail."""

    def __call__(self, *args, **kwargs):
        return self

    def __getattr__(self, name):
        if name.startswith('__'):
            raise AttributeError(name)
        return self

    def __setattr__(self, name, value):
        pass

    def __repr__(self):
        return '<unsupported plotting object>'


# ---------------------------------------------------------------------------
# Value checks
# ---------------------------------------------------------------------------

def _is_number(v):
    return isinstance(v, numbers.Real) and not isinstance(v, (bool, np.bool_))


def _is_color(v):
    if v is None:
        return True
    if not isinstance(v, str):
        return False
    from matplotlib.colors import is_color_like
    return is_color_like(v)


LEGEND_LOCATIONS = {
    'top_left': 'upper left', 'top_center': 'upper center',
    'top_right': 'upper right', 'center_left': 'center left',
    'center': 'center', 'center_right': 'center right',
    'bottom_left': 'lower left', 'bottom_center': 'lower center',
    'bottom_right': 'lower right',
    # Bokeh's short aliases
    'top': 'upper center', 'bottom': 'lower center',
    'left': 'center left', 'right': 'center right',
}
LINE_DASHES = ('solid', 'dashed', 'dashdot', 'dotted')
MARKERS = ('circle', 'x')
PALETTES = ('Spectral11',)


# ---------------------------------------------------------------------------
# Small recorded models
# ---------------------------------------------------------------------------

class _Model:
    """Attribute-guarded record.  ``_ATTRS`` maps each supported public
    attribute to a validator returning the value to store (or raising
    ``ValueError`` for an unsupported value)."""

    _ATTRS = {}
    _KIND = 'object'

    def __getattr__(self, name):
        # Only called for names not found normally.
        if name.startswith('_'):
            raise AttributeError(name)
        _unsupported("attribute %s.%s" % (self._KIND, name))
        return _Null()

    def __setattr__(self, name, value):
        if name.startswith('_'):
            object.__setattr__(self, name, value)
            return
        try:
            check = self._ATTRS.get(name)
            if check is None:
                _unsupported("attribute %s.%s" % (self._KIND, name))
                return
            try:
                value = check(self, value)
            except (ValueError, TypeError) as e:
                if isinstance(e, UnsupportedPlotFeature):
                    raise
                _unsupported("%s.%s = %r (%s)" % (self._KIND, name, value, e))
                return
            object.__setattr__(self, name, value)
        except Exception as e:
            _failed("%s.%s assignment" % (self._KIND, name), e)


def _check_end(self, v):
    if v is None:
        if isinstance(self, Range1d):
            raise ValueError("a Range1d end cannot be None")
        return None
    if not _is_number(v):
        raise ValueError("not a number")
    return float(v)


def _check_padding(self, v):
    if not _is_number(v) or v < 0:
        raise ValueError("must be a number >= 0")
    return float(v)


class _DataRange(_Model):
    """The automatic range a figure starts with (Bokeh's ``DataRange1d``):
    data bounds plus ``range_padding`` (a fraction of the span, split
    between the two ends); ``start``/``end`` override either end."""

    _KIND = 'DataRange1d'
    _ATTRS = {'start': _check_end, 'end': _check_end,
              'range_padding': _check_padding}

    def __init__(self):
        object.__setattr__(self, 'start', None)
        object.__setattr__(self, 'end', None)
        object.__setattr__(self, 'range_padding', 0.1)

    def __repr__(self):
        return 'DataRange1d(start=%r, end=%r, range_padding=%r)' % (
            self.start, self.end, self.range_padding)


class Range1d(_Model):
    """A fixed range: ``Range1d(start, end)`` or ``Range1d(start=, end=)``."""

    _KIND = 'Range1d'
    _ATTRS = {'start': _check_end, 'end': _check_end}

    def __init__(self, *args, **kwargs):
        object.__setattr__(self, 'start', 0.0)
        object.__setattr__(self, 'end', 1.0)
        try:
            bad = sorted(set(kwargs) - {'start', 'end'})
            if bad:
                _unsupported("Range1d keyword(s) %s" % ', '.join(bad))
            if len(args) > 2 or (args and ('start' in kwargs)) or \
                    (len(args) > 1 and 'end' in kwargs):
                _unsupported("Range1d arguments %r %r" % (args, kwargs))
                return
            vals = dict(zip(('start', 'end'), args))
            vals.update({k: kwargs[k] for k in ('start', 'end')
                         if k in kwargs})
            for k, v in vals.items():
                setattr(self, k, v)
        except Exception as e:
            _failed("Range1d()", e)

    def __repr__(self):
        return 'Range1d(start=%r, end=%r)' % (self.start, self.end)


class LinearAxis(_Model):
    """A second x axis: ``LinearAxis(x_range_name=, axis_label=)``, added
    with ``p.add_layout(axis, 'above')``."""

    _KIND = 'LinearAxis'
    _ATTRS = {}

    def __init__(self, **kwargs):
        object.__setattr__(self, 'x_range_name', 'default')
        object.__setattr__(self, 'axis_label', '')
        try:
            bad = sorted(set(kwargs) - {'x_range_name', 'axis_label'})
            if bad:
                _unsupported("LinearAxis keyword(s) %s" % ', '.join(bad))
            for k in ('x_range_name', 'axis_label'):
                if k in kwargs:
                    v = kwargs[k]
                    if not isinstance(v, str):
                        _unsupported("LinearAxis(%s=%r)" % (k, v))
                        continue
                    object.__setattr__(self, k, v)
        except Exception as e:
            _failed("LinearAxis()", e)

    def __repr__(self):
        return 'LinearAxis(x_range_name=%r, axis_label=%r)' % (
            self.x_range_name, self.axis_label)


def _check_location(self, v):
    if v not in LEGEND_LOCATIONS:
        raise ValueError("supported: %s" % ', '.join(sorted(LEGEND_LOCATIONS)))
    return v


class _Legend(_Model):
    _KIND = 'legend'
    _ATTRS = {'location': _check_location}

    def __init__(self):
        object.__setattr__(self, 'location', 'top_right')


def _check_grid_color(self, v):
    if v is not None:
        raise ValueError("only None is supported")
    return v


class _Grid(_Model):
    _KIND = 'xgrid'
    _ATTRS = {'grid_line_color': _check_grid_color}

    def __init__(self):
        object.__setattr__(self, 'grid_line_color', '#e5e5e5')


def _check_text(self, v):
    if v is None:
        return ''
    if not isinstance(v, str):
        raise ValueError("not a string")
    return v


class _Title(_Model):
    _KIND = 'title'
    _ATTRS = {'text': _check_text}

    def __init__(self, text=''):
        object.__setattr__(self, 'text', text)

    def __repr__(self):
        return 'Title(text=%r)' % self.text


# ---------------------------------------------------------------------------
# Figure recorder
# ---------------------------------------------------------------------------

_FIGURE_KW = {'title', 'x_axis_label', 'y_axis_label', 'plot_width',
              'plot_height', 'tooltips', 'y_range', 'y_axis_type'}

_LINE_KW = {'color', 'line_color', 'line_width', 'line_dash', 'alpha',
            'line_alpha', 'legend_label', 'x_range_name'}
_MARKER_KW = {'color', 'line_color', 'fill_color', 'line_width', 'alpha',
              'line_alpha', 'size', 'legend_label', 'x_range_name'}
_GLYPH_KW = {
    'line': _LINE_KW,
    'circle': _MARKER_KW,
    'diamond': _MARKER_KW,
    'scatter': _MARKER_KW | {'marker'},
    'quad': {'top', 'bottom', 'left', 'right', 'color', 'line_color',
             'fill_color', 'line_width', 'alpha', 'line_alpha',
             'legend_label', 'x_range_name'},
    'image': {'x', 'y', 'dw', 'dh', 'palette', 'level'},
}


def _check_kw(kind, name, v):
    """Validate one glyph keyword; return the value to record (copied) or
    raise ValueError."""
    if name in ('color', 'line_color', 'fill_color'):
        if not _is_color(v):
            raise ValueError("not a colour name")
        return v
    if name in ('line_width', 'size'):
        if not _is_number(v) or v < 0:
            raise ValueError("must be a number >= 0")
        return float(v)
    if name in ('alpha', 'line_alpha'):
        if not _is_number(v) or not 0 <= v <= 1:
            raise ValueError("must be a number in [0, 1]")
        return float(v)
    if name == 'line_dash':
        if v not in LINE_DASHES:
            raise ValueError("supported: %s" % ', '.join(LINE_DASHES))
        return v
    if name == 'marker':
        if v not in MARKERS:
            raise ValueError("supported: %s" % ', '.join(MARKERS))
        return v
    if name in ('legend_label', 'x_range_name'):
        if not isinstance(v, str):
            raise ValueError("not a string")
        return v
    if name == 'palette':
        if v not in PALETTES:
            raise ValueError("supported: %s" % ', '.join(PALETTES))
        return v
    if name == 'level':
        return v
    if name in ('top', 'bottom', 'left', 'right'):
        return _data(v)
    if name in ('x', 'y', 'dw', 'dh'):
        a = _data(v)
        if a.size != 1:
            raise ValueError("must be a single number")
        return float(a[0])
    raise ValueError("unknown")


def _data(v):
    """Copy one data argument to a 1-D float array."""
    if isinstance(v, (str, bytes)):
        raise ValueError("a string (ColumnDataSource column names are not "
                         "supported; pass arrays)")
    if np.ma.isMaskedArray(v):
        a = np.ma.filled(np.ma.array(v, dtype=float, copy=True), np.nan)
    else:
        a = np.array(v, dtype=float, copy=True)
    return np.atleast_1d(a).ravel()


class PlotFigure(_Model):
    """Records what a primitive does to a Bokeh ``figure``."""

    _KIND = 'figure'

    def __init__(self, **kwargs):
        object.__setattr__(self, '_title', _Title(''))
        object.__setattr__(self, '_x_axis_label', '')
        object.__setattr__(self, '_y_axis_label', '')
        object.__setattr__(self, '_plot_width', 600)
        object.__setattr__(self, '_plot_height', 600)
        object.__setattr__(self, '_y_log', False)
        object.__setattr__(self, '_x_range', _DataRange())
        object.__setattr__(self, '_y_range', _DataRange())
        object.__setattr__(self, '_extra_x_ranges', {})
        object.__setattr__(self, '_layouts', [])
        object.__setattr__(self, '_legend', _Legend())
        object.__setattr__(self, '_xgrid', _Grid())
        object.__setattr__(self, '_glyphs', [])
        self._configure(kwargs)

    def _configure(self, kw):
        bad = sorted(set(kw) - _FIGURE_KW)
        if bad:
            _unsupported("figure() keyword(s) %s" % ', '.join(bad))
        for name in ('title', 'x_axis_label', 'y_axis_label'):
            if name in kw:
                v = kw[name]
                if v is None:
                    v = ''
                if not isinstance(v, str):
                    _unsupported("figure(%s=%r)" % (name, v))
                    continue
                if name == 'title':
                    self._title.text = v
                else:
                    object.__setattr__(self, '_' + name, v)
        for name in ('plot_width', 'plot_height'):
            if name in kw:
                v = kw[name]
                if not _is_number(v) or v < 1 or int(v) != v:
                    _unsupported("figure(%s=%r)" % (name, v))
                    continue
                object.__setattr__(self, '_' + name, int(v))
        if 'y_axis_type' in kw:
            v = kw['y_axis_type']
            if v in ('log', 'linear'):
                object.__setattr__(self, '_y_log', v == 'log')
            else:
                _unsupported("figure(y_axis_type=%r)" % (v,))
        if 'y_range' in kw:
            v = kw['y_range']
            if isinstance(v, Range1d):
                object.__setattr__(self, '_y_range', v)
            elif isinstance(v, (list, tuple)) and len(v) == 2 and \
                    all(_is_number(e) for e in v):
                object.__setattr__(self, '_y_range',
                                   Range1d(float(v[0]), float(v[1])))
            else:
                _unsupported("figure(y_range=%r)" % (v,))
        # 'tooltips' is accepted and ignored: a PNG has no hover.

    def __repr__(self):
        return '<PlotFigure %r, %d glyphs>' % (self._title.text,
                                               len(self._glyphs))

    # -- properties ---------------------------------------------------------

    def _set_range(self, name, v):
        try:
            if not isinstance(v, Range1d):
                _unsupported("%s = %r (only Range1d)" % (name, v))
                return
            object.__setattr__(self, '_' + name, v)
        except Exception as e:
            _failed("%s assignment" % name, e)

    title = property(lambda self: self._title)
    legend = property(lambda self: self._legend)
    xgrid = property(lambda self: self._xgrid)
    x_range = property(lambda self: self._x_range,
                       lambda self, v: self._set_range('x_range', v))
    y_range = property(lambda self: self._y_range,
                       lambda self, v: self._set_range('y_range', v))

    @property
    def extra_x_ranges(self):
        return self._extra_x_ranges

    @extra_x_ranges.setter
    def extra_x_ranges(self, v):
        try:
            if not isinstance(v, dict) or len(v) > 1 or not all(
                    isinstance(k, str) and isinstance(r, Range1d)
                    for k, r in v.items()):
                _unsupported("extra_x_ranges = %r (one {name: Range1d} only)"
                             % (v,))
                return
            object.__setattr__(self, '_extra_x_ranges', dict(v))
        except Exception as e:
            _failed("extra_x_ranges assignment", e)

    _SETTABLE = ('x_range', 'y_range', 'extra_x_ranges')

    def __setattr__(self, name, value):
        if name.startswith('_') or name in self._SETTABLE:
            object.__setattr__(self, name, value)
            return
        _unsupported("attribute figure.%s" % name)

    # -- glyph methods ------------------------------------------------------

    def _glyph(self, kind, args, kwargs):
        try:
            rec = {'kind': kind}
            bad = sorted(set(kwargs) - _GLYPH_KW[kind])
            problems = ["keyword(s) %s" % ', '.join(bad)] if bad else []
            for k, v in kwargs.items():
                if k in bad:
                    continue
                try:
                    rec[k] = _check_kw(kind, k, v)
                except (ValueError, TypeError) as e:
                    problems.append("%s=%r (%s)" % (k, _short(v), e))
            fatal = None
            if kind == 'image':
                if len(args) != 1:
                    fatal = "image() takes one positional argument [array]"
                else:
                    try:
                        rec['image'] = _image_data(args[0])
                    except (ValueError, TypeError) as e:
                        fatal = "image data (%s)" % e
                for k in ('x', 'y', 'dw', 'dh', 'palette'):
                    if k not in rec and fatal is None:
                        fatal = "image() needs %s=" % k
            elif kind == 'quad':
                if args:
                    fatal = "quad() takes keywords only"
                elif not all(k in rec for k in ('top', 'bottom', 'left',
                                                 'right')):
                    fatal = "quad() needs top=, bottom=, left=, right="
                else:
                    try:
                        arrs = np.broadcast_arrays(
                            *(rec[k] for k in ('left', 'right', 'bottom',
                                               'top')))
                    except ValueError:
                        fatal = "quad() edges of different lengths"
                    else:
                        for k, a in zip(('left', 'right', 'bottom', 'top'),
                                        arrs):
                            rec[k] = np.array(a, copy=True)
            else:
                if len(args) != 2:
                    fatal = "%s() takes two positional arguments (x, y)" % kind
                else:
                    try:
                        x, y = _data(args[0]), _data(args[1])
                        x, y = np.broadcast_arrays(x, y)
                        rec['x'], rec['y'] = np.array(x), np.array(y)
                    except (ValueError, TypeError) as e:
                        fatal = "%s() data (%s)" % (kind, e)
            if fatal:
                problems.append(fatal)
            if problems:
                _unsupported("%s(): %s" % (kind, '; '.join(problems)) +
                             (" - glyph skipped" if fatal else ""))
            if fatal:
                return None
            self._glyphs.append(rec)
        except Exception as e:
            _failed("%s()" % kind, e)
        return None

    def line(self, *args, **kwargs):
        return self._glyph('line', args, kwargs)

    def circle(self, *args, **kwargs):
        return self._glyph('circle', args, kwargs)

    def diamond(self, *args, **kwargs):
        return self._glyph('diamond', args, kwargs)

    def scatter(self, *args, **kwargs):
        return self._glyph('scatter', args, kwargs)

    def quad(self, *args, **kwargs):
        return self._glyph('quad', args, kwargs)

    def image(self, *args, **kwargs):
        return self._glyph('image', args, kwargs)

    def add_layout(self, obj, place='center'):
        try:
            if not isinstance(obj, LinearAxis) or place != 'above':
                _unsupported("add_layout(%r, %r) (only a LinearAxis, "
                             "'above')" % (obj, place))
                return None
            self._layouts.append(copy.copy(obj))
        except Exception as e:
            _failed("add_layout()", e)
        return None


def _short(v):
    r = repr(v)
    return r if len(r) < 60 else r[:57] + '...'


def _image_data(v):
    if not isinstance(v, (list, tuple)) or len(v) != 1:
        raise ValueError("expected a list holding one 2-D array")
    a = v[0]
    if np.ma.isMaskedArray(a):
        a = np.ma.filled(np.ma.array(a, dtype=float, copy=True), np.nan)
    else:
        a = np.array(a, dtype=float, copy=True)
    if a.ndim != 2:
        raise ValueError("expected a 2-D array")
    return a


def figure(**kwargs):
    """Start recording a plot, as Bokeh's ``figure``.

    Supported keywords: ``title``, ``x_axis_label``, ``y_axis_label``,
    ``plot_width``, ``plot_height`` (pixels, default 600), ``y_range`` (a
    two-element list/tuple or a ``Range1d``), ``y_axis_type`` (``'log'``;
    ``'linear'`` is the default) and ``tooltips`` (ignored).
    """
    try:
        return PlotFigure(**kwargs)
    except UnsupportedPlotFeature:
        raise
    except Exception as e:
        _failed("figure()", e)
        return PlotFigure()


# ---------------------------------------------------------------------------
# Grid
# ---------------------------------------------------------------------------

class PlotGrid(_Model):
    """Records a ``gridplot``: the figures (by reference, so the same list
    can go into several grids), the number of columns and the cell size."""

    _KIND = 'gridplot'

    def __init__(self, children=(), ncols=1, cell_width=None,
                 cell_height=None):
        object.__setattr__(self, '_children', list(children))
        object.__setattr__(self, '_ncols', max(1, int(ncols)))
        object.__setattr__(self, '_cell_width', cell_width)
        object.__setattr__(self, '_cell_height', cell_height)

    def __setattr__(self, name, value):
        if name.startswith('_'):
            object.__setattr__(self, name, value)
            return
        _unsupported("attribute gridplot.%s" % name)

    def _cell_size(self):
        w, h = self._cell_width, self._cell_height
        if w is None:
            w = max([c._plot_width for c in self._children] or [600])
        if h is None:
            h = max([c._plot_height for c in self._children] or [600])
        return w, h

    def _shape(self):
        n = len(self._children)
        ncols = min(self._ncols, n) if n else 1
        nrows = max(1, int(math.ceil(n / ncols))) if n else 1
        return nrows, ncols

    def __repr__(self):
        return '<PlotGrid %d figures, ncols=%d>' % (len(self._children),
                                                    self._ncols)


def gridplot(children, ncols=None, plot_width=None, plot_height=None,
             toolbar_location=None, **kwargs):
    """Lay recorded figures out in a grid, as Bokeh's ``gridplot``.

    ``plot_width``/``plot_height`` are per cell and override each figure's
    own size; ``toolbar_location`` is accepted and ignored.  The saved
    image is ``ncols * plot_width`` by ``nrows * plot_height`` pixels.
    """
    grid = PlotGrid()
    try:
        if kwargs:
            _unsupported("gridplot() keyword(s) %s"
                         % ', '.join(sorted(kwargs)))
        kids = []
        for c in list(children):
            if isinstance(c, PlotFigure):
                kids.append(c)
            else:
                _unsupported("gridplot() child %r (only figures in a flat "
                             "list)" % (c,))
        if ncols is None or not _is_number(ncols) or ncols < 1:
            _unsupported("gridplot(ncols=%r)" % (ncols,))
            ncols = 1
        size = []
        for name, v in (('plot_width', plot_width),
                        ('plot_height', plot_height)):
            if v is not None and (not _is_number(v) or v < 1 or int(v) != v):
                _unsupported("gridplot(%s=%r)" % (name, v))
                v = None
            size.append(None if v is None else int(v))
        grid = PlotGrid(kids, int(ncols), size[0], size[1])
    except Exception as e:
        _failed("gridplot()", e)
    return grid


# ---------------------------------------------------------------------------
# Rendering and hand-off
# ---------------------------------------------------------------------------

def pixel_size(obj):
    """(width, height) in pixels of the image ``render`` will draw."""
    if isinstance(obj, PlotGrid):
        w, h = obj._cell_size()
        nrows, ncols = obj._shape()
        return ncols * w, nrows * h
    if isinstance(obj, PlotFigure):
        return obj._plot_width, obj._plot_height
    raise TypeError("not a recorded figure or grid: %r" % (obj,))


def plot_label(obj):
    """The title text (for a grid, its first figure's title)."""
    if isinstance(obj, PlotGrid):
        return obj._children[0]._title.text if obj._children else ''
    if isinstance(obj, PlotFigure):
        return obj._title.text
    return ''


def snapshot(obj):
    """An independent deep copy of a recorded figure or grid, for hand-off
    to another thread."""
    return copy.deepcopy(obj)


def render(obj, fig=None):
    """Draw a recorded figure or grid into a Matplotlib ``Figure``.

    With ``fig`` None a new ``Figure`` (100 dpi, Agg canvas) is made; with
    ``fig`` given it is cleared, resized to ``pixel_size(obj)`` and drawn
    into.  One ``Axes`` per recorded figure (row-major for a grid), plus one
    ``twiny`` ``Axes`` for a figure with ``extra_x_ranges``.
    """
    from kcwidrp.core import _plot_render
    if not isinstance(obj, (PlotFigure, PlotGrid)):
        raise TypeError("not a recorded figure or grid: %r" % (obj,))
    return _plot_render.render(obj, fig)


def save_plot(fig, filename=None):
    """Render a recorded figure or grid to ``plots/<filename>`` (relative to
    the current directory; default ``plots/kcwi_drp_plot.png``).

    Needs no viewer.  Never raises when ``STRICT`` is false: a failure is
    logged as a warning and the reduction carries on.  (Errors writing the
    file are only logged even when ``STRICT`` is true.)
    """
    try:
        if not isinstance(fig, (PlotFigure, PlotGrid)):
            _unsupported("save_plot(%r) (not a figure or gridplot)" % (fig,))
            return
        if filename is None:
            fnam = os.path.join('plots', 'kcwi_drp_plot.png')
        else:
            fnam = os.path.join('plots', filename)
        mfig = render(fig)
        try:
            os.makedirs(os.path.dirname(fnam) or '.', exist_ok=True)
            mfig.savefig(fnam, format='png', dpi=mfig.dpi)
        except OSError as e:
            logger.warning("Could not save plot to %s: %s" % (fnam, e))
            return
        logger.info(">>> Saving to %s" % fnam)
    except Exception as e:
        _failed("save_plot()", e)


def configure_plot_driver(firefox_compat=False, prewarm=False):
    """No-op.  PNGs are now written by Matplotlib, with no browser.

    Kept only because ``StartBokeh`` still imports and calls it; remove it
    when ``StartBokeh`` is replaced by ``StartPlotViewer`` (bokeh-plan 3.1).
    """
    return None


# ---------------------------------------------------------------------------
# Helpers used by the primitives (unchanged behaviour)
# ---------------------------------------------------------------------------

def get_plot_lims(data, padding=0.05, clip=True):
    """Get plot limits using data range plus padding fraction"""
    dmin = np.nanmin(data)
    dmax = np.nanmax(data)
    ddel = np.abs(dmax - dmin)
    if ddel < 0.0001 and clip:
        ddel = 1.0
    return dmin - padding * ddel, dmax + padding * ddel


def oplot_slices(fig, yrange):
    """Overplot slices vertically on plot"""
    for ix in range(1, 24):
        sx = ix * 5 - 0.5
        fig.line([sx, sx], yrange, color='black', line_dash='dashdot')


def set_plot_lims(fig, xlim=None, ylim=None):
    """Set figure plot ranges"""
    if xlim:
        fig.x_range.start = xlim[0]
        fig.x_range.end = xlim[1]
    if ylim:
        fig.y_range.start = ylim[0]
        fig.y_range.end = ylim[1]
