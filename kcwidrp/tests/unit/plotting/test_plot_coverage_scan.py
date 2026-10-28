"""Coverage scan for the plotting wrapper (bokeh plan 1.2; upgrade plan s7
amendment 5: a permanent gate).

An ``ast`` scan of every file in ``kcwidrp/primitives`` finds each method,
keyword, attribute chain and literal value that the code uses on a figure
object, and fails the test for anything the wrapper is specified not to
support, or that the scan cannot classify.  Its result is also compared, in
both directions, with the feature summary (section 2) of
``docs/plot_inventory.md``.

The supported set below is hard-coded from the table in section 1.1 of
``docs/bokeh_port_plan.md``.  It is deliberately NOT read from the wrapper,
so that the scan and the wrapper are two independent statements of the same
specification (the structural tests in ``test_kcwi_plotting.py`` check that
the wrapper really supports every item listed here).

The scan needs only ``ast``: it runs without the wrapper and without Bokeh.

Rules of the scan (from the plan):

* A name is a figure only from the point where it is assigned the result of a
  call to ``figure(...)`` (or another figure name), and stops being one when
  it is assigned anything else.  Tracking is per function, in source order,
  so ``p = Pool(); p.map(...); p.close(); p = figure(...)`` (MakeCube,
  FitCenter) does not flag ``map``/``close``.  Figure names are discovered,
  never listed.
* A figure passed to ``set_plot_lims``/``oplot_slices`` (any function defined
  in the scanned wrapper file) is followed into that function; ``bokeh_plot``,
  ``save_plot`` and ``gridplot`` (and the other wrapper entry points) are
  accepted sinks; a list that a figure is appended to may be given to
  ``gridplot``.
* Anything else that touches a figure (passed to an unknown function, stored
  in an attribute, returned, ``setattr``, ``**kwargs`` ...) is reported as
  unclassifiable, with file and line.  Unknown never passes silently.
"""
import ast
import re
import warnings
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[4]
PRIMITIVES = REPO / "kcwidrp" / "primitives"
WRAPPER = REPO / "kcwidrp" / "core" / "kcwi_plotting.py"
INVENTORY = REPO / "docs" / "plot_inventory.md"

# ---------------------------------------------------------------------------
# The supported set.
# source: published - docs/bokeh_port_plan.md section 1.1 (table "Bokeh usage"
# and the list of known mapping traps).
# ---------------------------------------------------------------------------
FIGURE_KW = {
    "title", "x_axis_label", "y_axis_label", "plot_width", "plot_height",
    "tooltips", "y_range", "y_axis_type",
}
_COMMON_GLYPH_KW = {
    "color", "line_color", "fill_color", "line_width", "line_dash", "alpha",
    "line_alpha", "size", "legend_label", "x_range_name",
}
GLYPH_KW = {
    "line": _COMMON_GLYPH_KW,
    "circle": _COMMON_GLYPH_KW,
    "diamond": _COMMON_GLYPH_KW,
    "scatter": _COMMON_GLYPH_KW | {"marker"},
    "quad": _COMMON_GLYPH_KW | {"top", "bottom", "left", "right"},
    "image": {"x", "y", "dw", "dh", "palette", "level"},
}
POSITIONAL_DATA = {"line": 2, "circle": 2, "diamond": 2, "scatter": 2,
                   "quad": 0, "image": 1}
LITERAL_KW = {  # keyword -> supported literal values
    "marker": {"x"},
    "line_dash": {"dashed", "dashdot", "dotted"},
    "palette": {"Spectral11"},
    "y_axis_type": {"log"},
}
LEGEND_LOCATIONS = {"top_left", "top_center", "top_right", "bottom_center"}
ATTR_ASSIGN = {  # attribute chain -> kind of value that is supported
    "x_range": "range1d", "y_range": "range1d",
    "x_range.start": "any", "x_range.end": "any",
    "y_range.start": "any", "y_range.end": "any",
    "x_range.range_padding": "zero", "y_range.range_padding": "zero",
    "legend.location": "legend_location",
    "xgrid.grid_line_color": "none",
    "extra_x_ranges": "range1d_dict",
}
ATTR_READ = {"title.text"}  # source: published - plan "Interfaces fixed ..."
RANGE1D_KW = {"start", "end"}
LINEAR_AXIS_KW = {"x_range_name", "axis_label"}
GRIDPLOT_KW = {"ncols", "plot_width", "plot_height", "toolbar_location"}
ADD_LAYOUT_POSITIONS = {"above"}
UNSUPPORTED_CALLS = {"ColumnDataSource"}  # rewritten away in plan task 2.2
# Wrapper entry points that take a figure (or grid) and only consume it.
SINKS = {"bokeh_plot", "save_plot", "render", "snapshot", "pixel_size",
         "plot_label"}
# Features in the inventory that plan task 2.2 removes from the primitives.
REMOVED_BY_2_2 = {"glyph.kw:line:source"}


class Scan:
    def __init__(self):
        self.features = set()
        self.problems = []  # (file, line, message)

    def feature(self, name):
        self.features.add(name)

    def problem(self, file, line, message):
        self.problems.append((str(file), line, message))

    def describe(self):
        return "\n".join("%s:%s: %s" % (Path(f).name, ln, m)
                         for f, ln, m in sorted(self.problems))


def _call_name(call):
    f = call.func
    if isinstance(f, ast.Name):
        return f.id
    if isinstance(f, ast.Attribute):
        return f.attr
    return None


def _flatten(node):
    """(root, [attr, ...], clean) of an attribute chain.

    root is the innermost expression; clean is False if a call or subscript
    sits between the root and the outermost attribute."""
    attrs = []
    clean = True
    while True:
        if isinstance(node, ast.Attribute):
            attrs.append(node.attr)
            node = node.value
        elif isinstance(node, (ast.Subscript, ast.Call)):
            clean = False
            node = node.value if isinstance(node, ast.Subscript) else node.func
        else:
            break
    return node, list(reversed(attrs)), clean


class _Analyzer:
    def __init__(self, scan, file, funcs, followed, figs=(), tag=None):
        self.scan = scan
        self.file = file
        self.funcs = funcs
        self.followed = followed
        self.figs = set(figs)
        self.lists = set()
        self.seen_calls = set()

    # -- helpers ----------------------------------------------------------
    def bad(self, node, msg):
        self.scan.problem(self.file, node.lineno, "unclassifiable: " + msg)

    def unsup(self, node, msg):
        self.scan.problem(self.file, node.lineno, "unsupported: " + msg)

    def is_figure_call(self, node):
        return isinstance(node, ast.Call) and _call_name(node) == "figure"

    def is_fig(self, node):
        return ((isinstance(node, ast.Name) and node.id in self.figs)
                or self.is_figure_call(node))

    # -- statements -------------------------------------------------------
    def run(self, body):
        for stmt in body:
            self.stmt(stmt)

    def stmt(self, node):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef,
                             ast.ClassDef)):
            return  # analysed as their own scope
        if isinstance(node, ast.Assign):
            if isinstance(node.value, (ast.List, ast.Tuple)) and \
                    node.value.elts and all(
                        self.is_fig(e) for e in node.value.elts):
                for e in node.value.elts:
                    if self.is_figure_call(e):
                        self.expr(e)
            else:
                self.expr(node.value)
            for t in node.targets:
                self.bind(t, node.value, node)
        elif isinstance(node, ast.AnnAssign):
            if node.value is not None:
                self.expr(node.value)
                self.bind(node.target, node.value, node)
        elif isinstance(node, ast.AugAssign):
            self.expr(node.value)
            root, attrs, clean = _flatten(node.target)
            if self.is_fig(root) and attrs:
                self.read_chain(node.target, root, attrs, clean)
                self.store_chain(node.target, root, attrs, clean, node.value)
            elif isinstance(node.target, ast.Name):
                self.figs.discard(node.target.id)
        elif isinstance(node, ast.For):
            lst = (isinstance(node.iter, ast.Name)
                   and node.iter.id in self.lists)
            if not lst:
                self.expr(node.iter)
            if isinstance(node.target, ast.Name):
                if lst:
                    self.figs.add(node.target.id)
                else:
                    self.figs.discard(node.target.id)
            else:
                self.unbind_all(node.target)
            self.run(node.body)
            self.run(node.orelse)
        elif isinstance(node, ast.With) or isinstance(node, ast.AsyncWith):
            for item in node.items:
                self.expr(item.context_expr)
                if item.optional_vars is not None:
                    self.unbind_all(item.optional_vars)
            self.run(node.body)
        elif isinstance(node, ast.Delete):
            for t in node.targets:
                if isinstance(t, ast.Name):
                    self.figs.discard(t.id)
        elif isinstance(node, (ast.Import, ast.ImportFrom, ast.Pass,
                               ast.Break, ast.Continue, ast.Global,
                               ast.Nonlocal)):
            return
        else:
            for field, value in ast.iter_fields(node):
                self.generic(value)

    def generic(self, value):
        if isinstance(value, list):
            for v in value:
                self.generic(v)
        elif isinstance(value, ast.stmt):
            self.stmt(value)
        elif isinstance(value, ast.expr):
            self.expr(value)
        elif isinstance(value, ast.ExceptHandler):
            if value.name:
                self.figs.discard(value.name)
            self.run(value.body)
        elif isinstance(value, ast.AST):
            for _f, v in ast.iter_fields(value):
                self.generic(v)

    def unbind_all(self, target):
        for n in ast.walk(target):
            if isinstance(n, ast.Name):
                self.figs.discard(n.id)
                self.lists.discard(n.id)

    def bind(self, target, value, stmt):
        if isinstance(target, ast.Name):
            name = target.id
            if self.is_fig(value):
                self.figs.add(name)
                self.lists.discard(name)
            elif isinstance(value, (ast.List, ast.Tuple)) and \
                    value.elts and all(self.is_fig(e) for e in value.elts):
                # the elements were handled as bare names by expr(); undo
                self.figs.discard(name)
                self.lists.add(name)
            else:
                self.figs.discard(name)
                self.lists.discard(name)
            return
        root, attrs, clean = _flatten(target)
        if self.is_fig(root) and attrs:
            self.store_chain(target, root, attrs, clean, value)
            return
        if isinstance(target, (ast.Tuple, ast.List)):
            self.unbind_all(target)
            return
        # Attribute/Subscript target on something that is not a figure.
        if self.is_fig(value):
            self.bad(stmt, "figure stored in %s"
                     % ast.unparse(target))
        self.expr_target(target)

    def expr_target(self, target):
        # the object being stored into is read
        for child in ast.iter_child_nodes(target):
            if isinstance(child, ast.expr):
                self.expr(child)

    # -- chains on a figure -------------------------------------------------
    def store_chain(self, target, root, attrs, clean, value):
        if not clean:
            self.bad(target, "attribute chain through call/subscript: "
                     + ast.unparse(target))
            self.expr(value)
            return
        if self.is_figure_call(root):
            self.expr(root)
        chain = ".".join(attrs)
        self.scan.feature("attr:" + chain)
        kind = ATTR_ASSIGN.get(chain)
        if kind is None:
            self.unsup(target, "attribute %s assigned on a figure" % chain)
            self.expr(value)
            return
        if kind == "range1d":
            if not (isinstance(value, ast.Call)
                    and _call_name(value) == "Range1d"):
                self.unsup(target, "%s assigned something other than "
                           "Range1d(...)" % chain)
        elif kind == "range1d_dict":
            if not (isinstance(value, ast.Dict) and all(
                    isinstance(v, ast.Call) and _call_name(v) == "Range1d"
                    for v in value.values)):
                self.unsup(target, "%s must be a dict of Range1d(...)"
                           % chain)
        elif kind == "zero":
            if not (isinstance(value, ast.Constant) and value.value == 0
                    and not isinstance(value.value, bool)):
                self.unsup(target, "%s assigned %s (only the literal 0 is "
                           "supported)" % (chain, ast.unparse(value)))
        elif kind == "none":
            if not (isinstance(value, ast.Constant) and value.value is None):
                self.unsup(target, "%s assigned %s (only None is supported)"
                           % (chain, ast.unparse(value)))
        elif kind == "legend_location":
            if not isinstance(value, ast.Constant):
                self.bad(target, "legend.location is not a literal: %s"
                         % ast.unparse(value))
            else:
                self.scan.feature("literal:legend.location=%s" % value.value)
                if value.value not in LEGEND_LOCATIONS:
                    self.unsup(target, "legend.location = %r"
                               % (value.value,))
        self.expr(value)

    def read_chain(self, node, root, attrs, clean):
        if not clean:
            self.bad(node, "attribute chain through call/subscript: "
                     + ast.unparse(node))
            return
        chain = ".".join(attrs)
        self.scan.feature("read:" + chain)
        if chain not in ATTR_READ:
            self.unsup(node, "attribute %s read on a figure" % chain)

    # -- expressions ----------------------------------------------------------
    def expr(self, node):
        if node is None:
            return
        if isinstance(node, ast.Name):
            if node.id in self.figs:
                self.bad(node, "figure %r used outside a recognised form"
                         % node.id)
            return
        if isinstance(node, ast.Call):
            return self.call(node)
        if isinstance(node, ast.Attribute):
            root, attrs, clean = _flatten(node)
            if self.is_fig(root):
                if self.is_figure_call(root):
                    self.expr(root)
                self.read_chain(node, root, attrs, clean)
                return
            self.expr(root)
            return
        if isinstance(node, ast.NamedExpr):
            self.expr(node.value)
            if isinstance(node.target, ast.Name):
                self.figs.discard(node.target.id)
            return
        if isinstance(node, (ast.Lambda,)):
            return
        for child in ast.iter_child_nodes(node):
            if isinstance(child, ast.expr):
                self.expr(child)
            else:
                self.generic(child)

    def call(self, node):
        func = node.func
        name = _call_name(node)
        root, attrs, clean = _flatten(func)
        # method on a figure: <fig>.<chain>(...)
        if isinstance(func, ast.Attribute) and self.is_fig(root):
            if self.is_figure_call(root):
                self.expr(root)
            self.fig_method(node, attrs, clean)
            return
        # append to a figure list
        if (isinstance(func, ast.Attribute) and func.attr == "append"
                and isinstance(func.value, ast.Name)
                and len(node.args) == 1 and self.is_fig(node.args[0])):
            self.lists.add(func.value.id)
            self.figs.discard(func.value.id)
            if self.is_figure_call(node.args[0]):
                self.expr(node.args[0])
            return
        if name == "figure" and isinstance(func, (ast.Name, ast.Attribute)):
            return self.figure_call(node)
        if name == "Range1d":
            return self.range1d_call(node)
        if name == "LinearAxis":
            return self.linearaxis_call(node)
        if name == "gridplot":
            return self.gridplot_call(node)
        if name in UNSUPPORTED_CALLS:
            self.scan.feature("call:" + name)
            self.unsup(node, "%s(...) is not supported" % name)
        # arguments that are figures or figure lists
        sink = name in SINKS and isinstance(func, (ast.Name, ast.Attribute))
        helper = (isinstance(func, ast.Name) and name in self.funcs
                  and not sink)
        self.expr_func(func)
        for i, a in enumerate(node.args):
            self.call_arg(node, a, i, None, sink, helper, name)
        for kw in node.keywords:
            if kw.arg is None:
                if self.is_fig(kw.value) or (isinstance(kw.value, ast.Name)
                                             and kw.value.id in self.lists):
                    self.bad(kw.value, "figure passed through **kwargs")
                self.expr(kw.value)
            else:
                self.call_arg(node, kw.value, None, kw.arg, sink, helper,
                              name)

    def expr_func(self, func):
        if isinstance(func, ast.Attribute):
            root, attrs, clean = _flatten(func)
            if not self.is_fig(root):
                self.expr(root)
        elif not isinstance(func, ast.Name):
            self.expr(func)

    def call_arg(self, call, arg, pos, kwname, sink, helper, fname):
        is_list = isinstance(arg, ast.Name) and arg.id in self.lists
        figlist_literal = (isinstance(arg, (ast.List, ast.Tuple)) and arg.elts
                           and all(self.is_fig(e) for e in arg.elts))
        if self.is_fig(arg):
            if self.is_figure_call(arg):
                self.expr(arg)
            if sink:
                return
            if helper:
                self.follow(call, fname, pos, kwname, arg)
                return
            self.bad(arg, "figure passed to unknown function %s()" % fname)
            return
        if is_list or figlist_literal:
            if figlist_literal:
                for e in arg.elts:
                    if self.is_figure_call(e):
                        self.expr(e)
            if fname in ("len", "gridplot") or sink:
                return
            self.bad(arg, "figure list passed to unknown function %s()"
                     % fname)
            return
        if isinstance(arg, ast.Starred):
            if self.is_fig(arg.value) or (isinstance(arg.value, ast.Name)
                                          and arg.value.id in self.lists):
                self.bad(arg, "figure passed through *args")
            self.expr(arg.value)
            return
        self.expr(arg)

    def follow(self, call, fname, pos, kwname, arg):
        for file, fdef in self.funcs[fname]:
            params = [a.arg for a in fdef.args.posonlyargs + fdef.args.args]
            if kwname is not None:
                param = kwname if kwname in params + [
                    a.arg for a in fdef.args.kwonlyargs] else None
            else:
                param = params[pos] if pos is not None and pos < len(params) \
                    else None
            if param is None:
                self.bad(call, "cannot map figure argument to a parameter "
                         "of %s()" % fname)
                continue
            key = (file, fdef.lineno, param)
            if key in self.followed:
                continue
            self.followed.add(key)
            sub = _Analyzer(self.scan, file, self.funcs, self.followed,
                            figs={param})
            sub.run(fdef.body)

    def fig_method(self, node, attrs, clean):
        method = ".".join(attrs)
        if not clean:
            self.bad(node, "call through attribute chain/subscript: "
                     + ast.unparse(node.func))
            return
        if len(attrs) != 1 or (attrs[0] not in GLYPH_KW
                               and attrs[0] != "add_layout"):
            self.scan.feature("method:" + method)
            self.unsup(node, "method %s() called on a figure" % method)
            for a in node.args:
                self.expr(a)
            for kw in node.keywords:
                self.expr(kw.value)
            return
        method = attrs[0]
        if method == "add_layout":
            self.scan.feature("method:add_layout")
            for kw in node.keywords:
                self.unsup(node, "add_layout keyword %s=" % kw.arg)
            if len(node.args) != 2:
                self.unsup(node, "add_layout needs (LinearAxis(...), 'above')")
            else:
                axis, place = node.args
                if not (isinstance(axis, ast.Call)
                        and _call_name(axis) == "LinearAxis"):
                    self.unsup(node, "add_layout of something other than "
                               "LinearAxis(...)")
                if not (isinstance(place, ast.Constant)
                        and isinstance(place.value, str)):
                    self.bad(node, "add_layout position is not a literal")
                else:
                    self.scan.feature("literal:add_layout=%s" % place.value)
                    if place.value not in ADD_LAYOUT_POSITIONS:
                        self.unsup(node, "add_layout position %r"
                                   % place.value)
            for a in node.args:
                self.expr(a)
            return
        self.scan.feature("glyph:" + method)
        supported = GLYPH_KW[method]
        for kw in node.keywords:
            if kw.arg is None:
                self.bad(kw.value, "**kwargs in %s()" % method)
                self.expr(kw.value)
                continue
            self.scan.feature("glyph.kw:%s:%s" % (method, kw.arg))
            if kw.arg not in supported:
                self.unsup(kw.value, "keyword %s= on %s()" % (kw.arg, method))
            self.literal(kw, kw.arg)
            self.expr(kw.value)
        npos = 0
        for i, a in enumerate(node.args):
            if isinstance(a, ast.Starred):
                self.bad(a, "*args in %s()" % method)
                self.expr(a.value)
                continue
            npos += 1
            if isinstance(a, ast.Constant) and isinstance(a.value, str):
                self.unsup(a, "string %r as a positional argument of %s() "
                           "(column-name form)" % (a.value, method))
            self.expr(a)
        if npos > POSITIONAL_DATA[method]:
            self.unsup(node, "%d positional arguments to %s() (supported: "
                       "%d)" % (npos, method, POSITIONAL_DATA[method]))

    def literal(self, kw, key):
        if key not in LITERAL_KW:
            return
        v = kw.value
        if not isinstance(v, ast.Constant):
            self.bad(kw.value, "%s= is not a literal: %s"
                     % (key, ast.unparse(v)))
            return
        self.scan.feature("literal:%s=%s" % (key, v.value))
        if v.value not in LITERAL_KW[key]:
            self.unsup(kw.value, "%s=%r" % (key, v.value))

    def figure_call(self, node):
        self.scan.feature("call:figure")
        for a in node.args:
            self.unsup(a, "positional argument to figure()")
            self.expr(a)
        for kw in node.keywords:
            if kw.arg is None:
                self.bad(kw.value, "**kwargs in figure()")
                self.expr(kw.value)
                continue
            self.scan.feature("figure.kw:" + kw.arg)
            if kw.arg not in FIGURE_KW:
                self.unsup(kw.value, "figure() keyword %s=" % kw.arg)
            self.literal(kw, kw.arg)
            self.expr(kw.value)

    def range1d_call(self, node):
        if node.args:
            self.scan.feature("Range1d:positional")
        if len(node.args) > 2:
            self.unsup(node, "Range1d with more than two positional args")
        for kw in node.keywords:
            if kw.arg is None:
                self.bad(kw.value, "**kwargs in Range1d()")
            else:
                self.scan.feature("Range1d.kw:" + kw.arg)
                if kw.arg not in RANGE1D_KW:
                    self.unsup(kw.value, "Range1d keyword %s=" % kw.arg)
        for a in node.args:
            self.expr(a)
        for kw in node.keywords:
            self.expr(kw.value)

    def linearaxis_call(self, node):
        self.scan.feature("call:LinearAxis")
        for a in node.args:
            self.unsup(a, "positional argument to LinearAxis()")
        for kw in node.keywords:
            if kw.arg is None:
                self.bad(kw.value, "**kwargs in LinearAxis()")
            else:
                self.scan.feature("LinearAxis.kw:" + kw.arg)
                if kw.arg not in LINEAR_AXIS_KW:
                    self.unsup(kw.value, "LinearAxis keyword %s=" % kw.arg)
            self.expr(kw.value)

    def gridplot_call(self, node):
        self.scan.feature("call:gridplot")
        if len(node.args) != 1:
            self.unsup(node, "gridplot needs exactly one positional "
                       "argument (the list of figures)")
        for a in node.args:
            if isinstance(a, ast.Name) and a.id in self.lists:
                continue
            if isinstance(a, (ast.List, ast.Tuple)) and a.elts and all(
                    self.is_fig(e) for e in a.elts):
                for e in a.elts:
                    if self.is_figure_call(e):
                        self.expr(e)
                continue
            self.bad(a, "gridplot argument is not a known list of figures: "
                     + ast.unparse(a))
            self.expr(a)
        for kw in node.keywords:
            if kw.arg is None:
                self.bad(kw.value, "**kwargs in gridplot()")
            else:
                self.scan.feature("gridplot.kw:" + kw.arg)
                if kw.arg not in GRIDPLOT_KW:
                    self.unsup(kw.value, "gridplot keyword %s=" % kw.arg)
            self.expr(kw.value)


def scan_sources(sources, wrapper_files=()):
    """Scan ``{filename: source text}``.

    Every function of every file is analysed on its own (module level code
    too).  Functions defined in a file named in ``wrapper_files`` are not
    analysed on their own (the wrapper's internals are its business); they
    are scanned only when a figure is passed to them.
    """
    scan = Scan()
    trees = {}
    for fname, text in sources.items():
        try:
            with warnings.catch_warnings():
                # docstrings in the primitives have "\*" escapes
                warnings.simplefilter("ignore", SyntaxWarning)
                trees[fname] = ast.parse(text, filename=str(fname))
        except SyntaxError as exc:
            scan.problem(fname, exc.lineno or 0, "unclassifiable: cannot "
                         "parse (%s)" % exc.msg)
    funcs = {}
    for fname, tree in trees.items():
        if fname not in wrapper_files:
            continue
        for n in ast.walk(tree):
            if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)):
                funcs.setdefault(n.name, []).append((fname, n))
    followed = set()
    for fname, tree in trees.items():
        if fname in wrapper_files:
            continue
        scopes = [tree.body]
        scopes += [n.body for n in ast.walk(tree)
                   if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]
        for body in scopes:
            _Analyzer(scan, fname, funcs, followed).run(body)
    return scan


def scan_paths(primitive_dir=PRIMITIVES, wrapper=WRAPPER):
    sources = {}
    for path in sorted(Path(primitive_dir).glob("*.py")):
        sources[str(path)] = path.read_text()
    wrapper_files = ()
    if Path(wrapper).exists():
        sources[str(wrapper)] = Path(wrapper).read_text()
        wrapper_files = (str(wrapper),)
    return scan_sources(sources, wrapper_files)


# ---------------------------------------------------------------------------
# The inventory's feature summary (docs/plot_inventory.md section 2), parsed
# into the same feature names the scan produces.
# ---------------------------------------------------------------------------
def _section(text, start, end):
    i = text.index(start)
    j = text.index(end, i)
    return text[i:j]


def _rows(section):
    for line in section.splitlines():
        if line.startswith("|") and not re.match(r"^\|[-| ]+\|$", line):
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            yield cells


def _ticks(cell):
    return re.findall(r"`([^`]*)`", cell)


def inventory_features(path=INVENTORY):
    text = Path(path).read_text()
    feats = set()
    sec21 = _section(text, "### 2.1 ", "### 2.2 ")
    for cells in list(_rows(sec21))[1:]:
        for t in _ticks(cells[0]):
            feats.add("figure.kw:" + t)
    sec22 = _section(text, "### 2.2 ", "### 2.3 ")
    for cells in list(_rows(sec22))[1:]:
        method = _ticks(cells[0])[0]
        feats.add("glyph:" + method)
        for kw in _ticks(cells[1]):
            feats.add("glyph.kw:%s:%s" % (method, kw))
    sec23 = _section(text, "### 2.3 ", "### 2.4 ")
    for cells in list(_rows(sec23))[1:]:
        for t in _ticks(cells[0]):
            if t.startswith("add_layout"):
                feats.add("method:add_layout")
            else:
                feats.add("attr:" + t)
    sec24 = _section(text, "### 2.4 ", "Session and document calls")
    for cells in list(_rows(sec24))[1:]:
        name = _ticks(cells[0])[0]
        usage = cells[2]
        if name == "Range1d":
            if "positional" in usage:
                feats.add("Range1d:positional")
            for kw in ("start", "end"):
                if kw + "=" in usage:
                    feats.add("Range1d.kw:" + kw)
        elif name == "LinearAxis":
            feats.add("call:LinearAxis")
            for kw in re.findall(r"(\w+)=", _ticks(usage)[0]):
                feats.add("LinearAxis.kw:" + kw)
        elif name == "gridplot":
            feats.add("call:gridplot")
            for kw in re.findall(r"(\w+)=", _ticks(usage)[0]):
                feats.add("gridplot.kw:" + kw)
        elif name == "figure":
            feats.add("call:figure")
        elif name == "ColumnDataSource":
            feats.add("call:ColumnDataSource")
    sec26 = _section(text, "### 2.6 ", "### 2.7 ")
    key_map = {"marker=": "marker", "line_dash=": "line_dash",
               "palette=": "palette", "y_axis_type=": "y_axis_type",
               "legend.location": "legend.location",
               "add_layout": "add_layout"}
    for cells in list(_rows(sec26))[1:]:
        label = _ticks(cells[0])
        key = key_map.get(label[0]) if label else None
        if key is None:
            key = key_map.get(cells[0])
        if key is None:
            continue
        for lit in _ticks(cells[1]):
            if lit[:1] not in "'\"":
                continue  # prose such as `oplot_slices`, not a literal value
            feats.add("literal:%s=%s" % (key, lit.strip("'\"")))
    return feats


# ---------------------------------------------------------------------------
# Tests on the real code base
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def real_scan():
    return scan_paths()


def test_primitives_use_only_supported_plot_features(real_scan):
    """The gate: every figure feature used in kcwidrp/primitives is supported."""
    # source: published - plan 1.2 "Coverage scan" and upgrade plan s7
    # amendment 5. Before plan phase 2 the only failures are the items task
    # 2.2 removes (MakeInvsens ColumnDataSource and the column-name line).
    assert real_scan.problems == [], "\n" + real_scan.describe()


def test_scan_finds_the_known_figures(real_scan):
    """Guards the scan itself: a scan that sees nothing must not pass."""
    # source: published - docs/plot_inventory.md section 1 (48 figures in 16
    # primitives); the scan must have seen all six glyph methods, the
    # twin-axis idiom and the grid.
    for feat in ("glyph:line", "glyph:circle", "glyph:scatter",
                 "glyph:diamond", "glyph:quad", "glyph:image",
                 "method:add_layout", "attr:extra_x_ranges", "call:gridplot",
                 "attr:y_range.start", "attr:legend.location"):
        assert feat in real_scan.features


def test_scan_matches_inventory_feature_summary(real_scan):
    """Both directions: scan result == docs/plot_inventory.md section 2."""
    # source: published - docs/plot_inventory.md section 2 (parsed, not copied)
    inventory = inventory_features()
    # the scan sees the wrapper's calls in the code; ignore feature kinds the
    # inventory does not list (reads of attributes, unlisted calls)
    scanned = {f for f in real_scan.features if not f.startswith("read:")}
    only_scan = scanned - inventory
    only_inventory = inventory - scanned - REMOVED_BY_2_2
    # Features the inventory lists but do not belong to figure usage in the
    # scanned files (ColumnDataSource is removed by task 2.2).
    only_inventory -= {"call:ColumnDataSource"}
    only_scan -= {"call:ColumnDataSource"}
    assert (sorted(only_scan), sorted(only_inventory)) == ([], [])


# ---------------------------------------------------------------------------
# Tests of the scan itself: mutations of a primitive, run on source text
# ---------------------------------------------------------------------------
BASE = '''
from bokeh.plotting import figure

def f(self):
    p = figure(title="t", x_axis_label="x", y_axis_label="y",
               plot_width=1, plot_height=1)
    p.line(x, y, color="red", legend_label="a")
    bokeh_plot(p, self.context.bokeh_session)
'''


def _scan(src, wrapper_src=None):
    sources = {"prim.py": src}
    wrapper_files = ()
    if wrapper_src is not None:
        sources["wrap.py"] = wrapper_src
        wrapper_files = ("wrap.py",)
    return scan_sources(sources, wrapper_files)


def _add(line):
    return BASE + "    " + line + "\n"


def test_scan_accepts_the_base_example():
    assert _scan(BASE).problems == []


@pytest.mark.parametrize("line", [
    "p.wedge(x, y)",
    "p.line(x, y, foo=1)",
    "p.circle(x, y, source=s)",
    "p.ygrid.grid_line_color = None",
    "p.xgrid.grid_line_color = 'red'",
    "p.title.text = 'new'",
    "p.legend.click_policy = 'hide'",
    "p.x_range = [0, 1]",
    "p.x_range.range_padding = 0.1",
    "p.legend.location = 'center_left'",
    "p.line(x, y, line_dash='longdash')",
    "p.scatter(x, y, marker='triangle')",
    "p.image([a], x=0, y=0, dw=1, dh=1, palette='Greys9')",
    "p.add_layout(LinearAxis(x_range_name='a'), 'below')",
    "p.add_layout(LinearAxis(bogus=1), 'above')",
    "p.quad(0, 1, 2, 3)",
    "p.line('x', 'y')",
    "p.line(x, y, z)",
    "p.circle(x, y, **kw)",
    "p.line(*args)",
    "p.extra_x_ranges = {'a': 5}",
    "r = Range1d(start=1, stop=2)",
    "g = gridplot([p], ncols=1, bogus=2)",
    "p.legend[0].location = 'top_left'",
    "p.legend.items.append(1)",
    "setattr(p, 'x_range', 1)",
    "helper(p)",
    "self.fig = p",
    "return p",
    "column(children=[p])",
])
def test_scan_flags_unsupported_use(line):
    # source: published - plan 1.2 Verification: each of these must fail
    scan = _scan(_add(line))
    assert scan.problems, "not flagged: " + line
    # the report names the line of the offending statement
    last = BASE.count("\n") + 1
    assert any(ln == last for _f, ln, _m in scan.problems)


def test_scan_flags_unsupported_figure_keyword():
    scan = _scan(BASE.replace('plot_height=1)', "plot_height=1,\n"
                              "               x_axis_type='log')"))
    assert scan.problems


def test_scan_flags_y_axis_type_other_than_log():
    assert not _scan(BASE.replace("plot_height=1)", "plot_height=1, "
                                  "y_axis_type='log')")).problems
    assert _scan(BASE.replace("plot_height=1)", "plot_height=1, "
                              "y_axis_type='datetime')")).problems


def test_scan_flags_figure_bound_to_a_new_name():
    # source: published - plan 1.2: "a figure bound to a new name"
    scan = _scan(_add("q = figure(title='q')\n    q.wedge(x, y)"))
    assert [m for _f, _l, m in scan.problems if "wedge" in m]


def test_scan_flags_wedge_on_a_discovered_figure_with_other_name():
    # figure names are discovered, not listed
    src = BASE.replace("p = figure", "my_fig = figure").replace(
        "p.line", "my_fig.line").replace("bokeh_plot(p", "bokeh_plot(my_fig")
    assert _scan(src).problems == []
    assert _scan(src + "    my_fig.wedge(x, y)\n").problems


def test_scan_flags_column_name_form_and_columndatasource():
    src = BASE + ("    s = ColumnDataSource(data=dict(x=w, y=v))\n"
                  "    p.line('x', 'y', line_color='black', source=s)\n")
    msgs = [m for _f, _l, m in _scan(src).problems]
    assert any("ColumnDataSource" in m for m in msgs)
    assert any("column-name" in m for m in msgs)
    assert any("source=" in m for m in msgs)


def test_scan_flags_chained_assignment_target_of_unsupported_attribute():
    # both targets of a chained assignment are collected
    src = BASE + "    p.x_range.range_padding = p.y_range.padding_bad = 0\n"
    msgs = [m for _f, _l, m in _scan(src).problems]
    assert any("y_range.padding_bad" in m for m in msgs)
    ok = BASE + "    p.x_range.range_padding = p.y_range.range_padding = 0\n"
    sc = _scan(ok)
    assert sc.problems == []
    assert {"attr:x_range.range_padding",
            "attr:y_range.range_padding"} <= sc.features


def test_scan_does_not_flag_pool_map_and_close():
    # source: published - plan 1.2: p.map / p.close on the pools in MakeCube
    # and FitCenter are not figure methods
    src = '''
from multiprocessing import Pool
from bokeh.plotting import figure

def f(self):
    p = Pool()
    out = p.map(work, items)
    p.close()
    p = figure(title="t")
    p.line(x, y)
    bokeh_plot(p, s)
'''
    assert _scan(src).problems == []
    # the same calls on a figure are flagged
    assert _scan(src + "    p.map(work, items)\n    p.close()\n").problems


def test_scan_figure_name_stops_being_a_figure_when_rebound():
    src = '''
def f(self):
    p = figure(title="t")
    p.line(x, y)
    p = Pool()
    p.map(work, items)
'''
    assert _scan(src).problems == []


def test_scan_alias_of_a_figure_is_a_figure():
    src = BASE + "    q = p\n    q.wedge(x, y)\n"
    assert _scan(src).problems


def test_scan_flags_figure_in_loop_over_figure_list():
    src = '''
def f(self):
    plts = []
    for i in range(3):
        p = figure(title="t")
        p.line(x, y)
        plts.append(p)
    for q in plts:
        q.wedge(x, y)
    bokeh_plot(gridplot(plts, ncols=2, plot_width=1, plot_height=1,
                        toolbar_location=None), s)
'''
    sc = _scan(src)
    assert len(sc.problems) == 1 and "wedge" in sc.problems[0][2]


def test_scan_gridplot_of_figure_list_and_len_are_accepted():
    src = '''
def f(self):
    plts = []
    p = figure(title="t")
    p.line(x, y)
    plts.append(p)
    if len(plts) > 0:
        save_plot(gridplot(plts, ncols=2, plot_width=1, plot_height=1,
                           toolbar_location=None), filename="a.png")
'''
    assert _scan(src).problems == []


def test_scan_follows_a_figure_into_wrapper_helpers():
    wrapper = '''
def set_plot_lims(fig, xlim=None, ylim=None):
    fig.x_range.start = xlim[0]
    fig.x_range.end = xlim[1]

def oplot_slices(fig, yrange):
    fig.line([1, 1], yrange, color="black", line_dash="dashdot")

def bad_helper(fig):
    fig.wedge(1, 2)

def indirect(fig):
    set_plot_lims(fig, xlim=[0, 1])
    fig.ygrid.grid_line_color = None
'''
    ok = BASE + "    set_plot_lims(p, xlim=xl)\n    oplot_slices(p, yr)\n"
    sc = _scan(ok, wrapper)
    assert sc.problems == []
    assert {"attr:x_range.start", "attr:x_range.end",
            "literal:line_dash=dashdot"} <= sc.features
    # the helper's own use of the figure is judged, with the helper's line
    sc = _scan(BASE + "    bad_helper(p)\n", wrapper)
    assert [(Path(f).name, m) for f, _l, m in sc.problems] == [
        ("wrap.py", "unsupported: method wedge() called on a figure")]
    # helpers calling helpers are followed
    sc = _scan(BASE + "    indirect(p)\n", wrapper)
    assert len(sc.problems) == 1 and "ygrid" in sc.problems[0][2]
    # keyword form of the figure argument is followed too
    sc = _scan(BASE + "    bad_helper(fig=p)\n", wrapper)
    assert len(sc.problems) == 1


def test_scan_does_not_analyse_wrapper_functions_on_their_own():
    wrapper = '''
def internals():
    f = figure(**kw)
    f.anything()
'''
    assert _scan(BASE, wrapper).problems == []


def test_scan_records_literal_values_and_forms():
    src = '''
def f(self):
    p = figure(y_axis_type='log', y_range=[0, 1])
    p.scatter(x, y, marker='x', line_dash='dotted')
    p.x_range = Range1d(1, 2)
    p.y_range = Range1d(start=1, end=2)
    p.legend.location = "top_left"
    p.add_layout(LinearAxis(x_range_name="pix", axis_label="a"), 'above')
    p.image([a], x=0, y=0, dw=1, dh=1, palette="Spectral11", level="image")
'''
    sc = _scan(src)
    assert sc.problems == []
    # source: analytic - the literals written in the snippet above
    assert {"literal:y_axis_type=log", "literal:marker=x",
            "literal:line_dash=dotted", "literal:legend.location=top_left",
            "literal:add_layout=above", "literal:palette=Spectral11",
            "Range1d:positional", "Range1d.kw:start", "Range1d.kw:end",
            "LinearAxis.kw:axis_label", "figure.kw:y_range"} <= sc.features


def test_scan_non_literal_where_a_literal_is_required_is_unclassifiable():
    msgs = [m for _f, _l, m in _scan(_add("p.line(x, y, line_dash=style)"))
            .problems]
    assert msgs and msgs[0].startswith("unclassifiable")


def test_scan_inventory_parser_reads_the_expected_sections():
    # source: published - docs/plot_inventory.md section 2 (spot checks of
    # the parser, so a parser that returns nothing cannot make the
    # comparison pass)
    feats = inventory_features()
    for f in ("figure.kw:y_axis_type", "glyph.kw:quad:fill_color",
              "glyph.kw:image:palette", "attr:xgrid.grid_line_color",
              "attr:x_range.range_padding", "method:add_layout",
              "LinearAxis.kw:axis_label", "gridplot.kw:toolbar_location",
              "Range1d:positional", "Range1d.kw:start",
              "literal:line_dash=dotted", "literal:legend.location=top_center",
              "literal:marker=x", "literal:y_axis_type=log",
              "literal:palette=Spectral11", "literal:add_layout=above"):
        assert f in feats, f
