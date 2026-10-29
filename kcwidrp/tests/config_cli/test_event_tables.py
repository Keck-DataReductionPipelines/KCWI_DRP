"""Event tables (plan T.6): the ``event_table`` dict literals of
``Kcwi_pipeline`` (kcwi_pipeline.py) and ``Keck_RTI_Pipeline``
(keck_rti_pipeline.py) are consistent.

Each entry is ``event_name: (action_name, new_state, next_event)``
(keckdrpframework ``BasePipeline`` docstring). The tests parse the dict
literal with ``ast``, because Python silently keeps the last of two equal
keys in a dict display.

How the framework resolves ``action_name`` (``BasePipeline._get_action``,
``_find_import_action``), which the resolver below follows:

1. a method of the pipeline class (``action_planner``,
   ``add_to_dataframe_only``): called directly, it needs no ``_perform``;
2. a name in the pipeline's module namespace (both modules do
   ``from kcwidrp.primitives.kcwi_file_primitives import *``): a class there
   is wrapped (``klass(action, context).apply()``);
3. a module ``<prefix>.<action_name>`` for each prefix in the framework
   config's ``primitive_path`` (``kcwidrp/configs/framework.cfg``:
   ``"kcwidrp.primitives", ""``) holding a class with the CamelCased module
   name, or a class named ``action_name`` in the package ``<prefix>``.
   keckdrpframework's own primitives are NOT on that path.

A class that is run through ``apply()`` must have ``_perform``
(``BasePrimitive.apply`` calls it).

"Next event": ``None`` ends the chain; a string must be a key of the same
table (the framework pushes it on the high-priority queue, and
``event_to_action`` looks it up in ``event_table``, falling back to a logged
"not found" no-op, which would silently end the chain).

Reachability: an event is reachable if it is an entry point or the next event
of a reachable event. Entry points are the string literals passed to
``push_event`` / ``append_event`` / ``push_hi_event`` anywhere in kcwidrp,
the values assigned to ``default_ingestion_event`` in the scripts, and the
framework default ``default_ingestion_event`` from framework.cfg. Unreachable
events are listed in the failure message of
``test_unreachable_events_are_the_known_ones`` and in the warning that test
issues.
"""
import ast
import importlib
import inspect
import pathlib
import warnings

import pytest
from keckdrpframework.config.framework_config import ConfigClass
from keckdrpframework.pipelines.base_pipeline import BasePipeline

import kcwidrp
from kcwidrp.core.kcwi_pkg_resources import get_resource_path

PKG_ROOT = pathlib.Path(kcwidrp.__file__).resolve().parent

PIPELINES = {
    "kcwi_pipeline": ("kcwidrp.pipelines.kcwi_pipeline", "Kcwi_pipeline"),
    "keck_rti_pipeline": ("kcwidrp.pipelines.keck_rti_pipeline",
                          "Keck_RTI_Pipeline"),
}

# Events the framework itself defines (BasePipeline.event_table0); they may be
# pushed or used as next events without being in the pipeline's table.
FRAMEWORK_EVENTS = {"noop", "echo", "no_event", "info"}

# Actions that are methods of the pipeline class, not primitives.
PIPELINE_METHOD_ACTIONS = {"action_planner", "add_to_dataframe_only"}

# Events that are in a table but cannot be reached from any entry point.
# Each needs a reason. Keyed by pipeline module.
KNOWN_UNREACHABLE = {
    "kcwi_pipeline": {
        "next_file_stop": "last entry of the table; nothing pushes it, no "
                          "event has it as next event, no script uses it",
        "start_bokeh": "bokeh plan 3.1: alias of start_plot_viewer kept for "
                       "old callers; the scripts now push start_plot_viewer",
    },
    "keck_rti_pipeline": {
        "next_file_stop": "same as in kcwi_pipeline",
        "start_bokeh": "same as in kcwi_pipeline",
    },
}


# --- parsing ---------------------------------------------------------------------

def _source_tree(module_name):
    path = PKG_ROOT.parent.joinpath(*module_name.split(".")).with_suffix(".py")
    import warnings as _w
    with _w.catch_warnings():
        _w.simplefilter("ignore", SyntaxWarning)
        return path, ast.parse(path.read_text(), filename=str(path))


def _event_table_node(module_name, class_name):
    path, tree = _source_tree(module_name)
    for cls in ast.walk(tree):
        if isinstance(cls, ast.ClassDef) and cls.name == class_name:
            for stmt in cls.body:
                if (isinstance(stmt, ast.Assign)
                        and any(isinstance(t, ast.Name)
                                and t.id == "event_table"
                                for t in stmt.targets)):
                    assert isinstance(stmt.value, ast.Dict), (
                        f"{path}:{stmt.lineno}: event_table is not a dict "
                        f"literal")
                    return path, stmt.value
    raise AssertionError(f"{path}: no event_table in class {class_name}")


def _entries(module_name, class_name):
    """Return [(key, lineno, (action, state, next))] in file order."""
    path, node = _event_table_node(module_name, class_name)
    out = []
    for k, v in zip(node.keys, node.values):
        assert isinstance(k, ast.Constant) and isinstance(k.value, str), (
            f"{path}:{k.lineno}: event key is not a string literal")
        assert isinstance(v, ast.Tuple) and len(v.elts) == 3, (
            f"{path}:{v.lineno}: event {k.value!r} is not a 3-tuple")
        parts = []
        for e in v.elts:
            assert isinstance(e, ast.Constant), (
                f"{path}:{e.lineno}: event {k.value!r} has a non-literal "
                f"field")
            assert e.value is None or isinstance(e.value, str)
            parts.append(e.value)
        out.append((k.value, k.lineno, tuple(parts)))
    return out


def _table(module_name, class_name):
    return {k: v for k, _, v in _entries(module_name, class_name)}


# --- resolution (follows keckdrpframework BasePipeline) ------------------------

def _camel(name):
    # BasePipeline._find_import_action.to_camel_case
    out, flag = [], True
    for c in name:
        if flag:
            out.append(c.upper())
            flag = False
            continue
        if c == "_":
            flag = True
            continue
        out.append(c)
    return "".join(out)


def _primitive_path():
    cfg = ConfigClass(get_resource_path("kcwidrp", "configs/framework.cfg"))
    return cfg.primitive_path


def _resolve(module_name, class_name, action):
    """Return ('method'|'class'|'function', object) or None."""
    module = importlib.import_module(module_name)
    pipeline = getattr(module, class_name)
    if callable(getattr(pipeline, action, None)):
        return "method", getattr(pipeline, action)
    if hasattr(module, action):
        obj = getattr(module, action)
        return ("class" if isinstance(obj, type) else "function"), obj
    for prefix in _primitive_path():
        full = f"{prefix}.{action}" if prefix else action
        package = ".".join(full.split(".")[:-1])
        last = full.split(".")[-1]
        for mod_name, attr in ((full, _camel(last)), (package, last)):
            if not mod_name:
                continue
            try:
                mod = importlib.import_module(mod_name)
            except Exception:
                continue
            if hasattr(mod, attr):
                return "class", getattr(mod, attr)
    return None


def _framework_finds_action(module_name, class_name, action):
    """The framework's own lookup, as a cross-check of _resolve."""
    from types import SimpleNamespace
    import logging
    module = importlib.import_module(module_name)
    pipeline_cls = getattr(module, class_name)
    config = SimpleNamespace(primitive_path=_primitive_path(),
                             print_trace=False)
    context = SimpleNamespace(logger=logging.getLogger("t6"), config=config,
                              pipeline_logger=logging.getLogger("t6"))
    pipeline = pipeline_cls(context)
    return pipeline._get_action("", action) is not None


# --- entry points ------------------------------------------------------------------

def _pushed_event_literals():
    """{event_name: [file:line]} for literal first args of push/append."""
    found = {}
    for path in sorted(PKG_ROOT.rglob("*.py")):
        if "tests" in path.relative_to(PKG_ROOT).parts:
            continue
        _, tree = _source_tree(
            ".".join(path.relative_to(PKG_ROOT.parent).with_suffix("").parts))
        for node in ast.walk(tree):
            if (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute)
                    and node.func.attr in ("push_event", "append_event",
                                           "push_hi_event")
                    and node.args):
                arg = node.args[0]
                if isinstance(arg, ast.Constant) and isinstance(arg.value,
                                                                 str):
                    found.setdefault(arg.value, []).append(
                        f"{path.relative_to(PKG_ROOT.parent)}:{node.lineno}")
    return found


def _default_ingestion_events():
    names = {}
    cfg = ConfigClass(get_resource_path("kcwidrp", "configs/framework.cfg"))
    names[cfg.default_ingestion_event] = ["configs/framework.cfg"]
    for path in sorted((PKG_ROOT / "scripts").glob("*.py")):
        _, tree = _source_tree(f"kcwidrp.scripts.{path.stem}")
        for node in ast.walk(tree):
            if (isinstance(node, ast.Assign)
                    and any(isinstance(t, ast.Attribute)
                            and t.attr == "default_ingestion_event"
                            for t in node.targets)
                    and isinstance(node.value, ast.Constant)):
                names.setdefault(node.value.value, []).append(
                    f"{path.relative_to(PKG_ROOT.parent)}:{node.lineno}")
    return names


def _entry_points():
    roots = {}
    for source in (_pushed_event_literals(), _default_ingestion_events()):
        for name, where in source.items():
            roots.setdefault(name, []).extend(where)
    return roots


def _reachable(table, roots):
    seen = set()
    todo = [r for r in roots if r in table]
    while todo:
        ev = todo.pop()
        if ev in seen:
            continue
        seen.add(ev)
        nxt = table[ev][2]
        if nxt in table:
            todo.append(nxt)
    return seen


# --- tests -------------------------------------------------------------------------

@pytest.fixture(params=sorted(PIPELINES), ids=sorted(PIPELINES))
def pipeline(request):
    return request.param, *PIPELINES[request.param]


@pytest.mark.parametrize("key", [
    pytest.param("kcwi_pipeline", id="kcwi_pipeline"),
    pytest.param("keck_rti_pipeline", id="keck_rti_pipeline"),
])
def test_no_duplicate_event_keys(key):
    module_name, class_name = PIPELINES[key]
    keys = [k for k, _, _ in _entries(module_name, class_name)]
    dups = sorted({k for k in keys if keys.count(k) > 1})
    # source: specification (plan T.6): a duplicate key is silently dropped
    # by Python, keeping the last entry
    assert not dups, f"duplicate event keys in {module_name}: {dups}"


def test_every_action_resolves_to_a_primitive(pipeline):
    _, module_name, class_name = pipeline
    problems = []
    for key, line, (action, _state, _next) in _entries(module_name,
                                                       class_name):
        found = _resolve(module_name, class_name, action)
        if found is None:
            problems.append(f"{key} (line {line}): action {action!r} "
                            f"resolves to nothing")
            continue
        kind, obj = found
        if kind == "method":
            if action not in PIPELINE_METHOD_ACTIONS:
                problems.append(f"{key} (line {line}): {action!r} is a "
                                f"pipeline method not on the allow-list")
        elif kind == "class":
            if not callable(getattr(obj, "_perform", None)):
                problems.append(f"{key} (line {line}): class {action} has "
                                f"no _perform")
        else:
            problems.append(f"{key} (line {line}): {action!r} is a plain "
                            f"function; expected a primitive class")
        # the framework must find the same thing
        if not _framework_finds_action(module_name, class_name, action):
            problems.append(f"{key} (line {line}): the framework's own "
                            f"lookup does not find {action!r}")
    # source: specification (plan T.6) and BasePipeline._get_action
    assert not problems, "\n".join(problems)


def test_pipeline_method_allow_list_is_exact(pipeline):
    _, module_name, class_name = pipeline
    actions = {a for _, _, (a, _, _) in _entries(module_name, class_name)}
    cls = getattr(importlib.import_module(module_name), class_name)
    for name in PIPELINE_METHOD_ACTIONS:
        assert name in actions, f"stale PIPELINE_METHOD_ACTIONS: {name}"
        assert inspect.isfunction(cls.__dict__.get(name)), (
            f"{name} is not defined in {class_name}")


def test_every_next_event_exists(pipeline):
    _, module_name, class_name = pipeline
    table = _table(module_name, class_name)
    problems = []
    for key, line, (_a, _s, nxt) in _entries(module_name, class_name):
        if nxt is None:
            continue
        if nxt not in table:
            note = (" (a framework event, but it would end the chain)"
                    if nxt in FRAMEWORK_EVENTS else "")
            problems.append(f"{key} (line {line}): next event {nxt!r} is "
                            f"not in the table{note}")
    # source: specification (plan T.6); BasePipeline.event_to_action turns an
    # unknown event into a logged no-op
    assert not problems, "\n".join(problems)


def test_pushed_event_names_exist(pipeline):
    _, module_name, class_name = pipeline
    table = _table(module_name, class_name)
    roots = _entry_points()
    # Events pushed by name by kcwidrp code must exist in a pipeline table or
    # be framework events. Checked against the union of both tables, because
    # the scripts for one pipeline push into that pipeline only, but the
    # tables are meant to agree (see test_tables_agree).
    unknown = {n: w for n, w in roots.items()
               if n not in table and n not in FRAMEWORK_EVENTS}
    # source: specification (plan T.6)
    assert not unknown, f"events pushed by name but not in the table: {unknown}"


def test_unreachable_events_are_the_known_ones(pipeline):
    key, module_name, class_name = pipeline
    table = _table(module_name, class_name)
    roots = _entry_points()
    unreachable = sorted(set(table) - _reachable(table, roots))
    known = KNOWN_UNREACHABLE[key]
    if unreachable:
        warnings.warn(f"{module_name}: unreachable events: {unreachable}",
                      stacklevel=1)
    # source: specification (plan T.6): unreachable events are reported; the
    # test fails only for ones not on the commented allow-list
    assert set(unreachable) == set(known), (
        f"{module_name}: unreachable events {unreachable}, allow-listed "
        f"{sorted(known)}; entry points: {sorted(roots)}")


def test_tables_differ_only_by_the_rti_alert():
    # The RTI pipeline is a copy of the main one plus an `alert_rti` event
    # (SendHTTP) after the two flux calibration events.
    a = _table(*PIPELINES["kcwi_pipeline"])
    b = _table(*PIPELINES["keck_rti_pipeline"])
    only_b = {k for k in b if k not in a}
    changed = {k for k in a if k in b and a[k] != b[k]}
    # source: characterization @ upgrade-start
    assert only_b == {"alert_rti"}
    assert set(a) - set(b) == set()
    assert changed == {"object_flux_calibrate", "nandshuff_flux_calibrate"}
    for k in changed:
        assert a[k][:2] == b[k][:2] and a[k][2] is None
        assert b[k][2] == "alert_rti"


def test_event_counts():
    # source: characterization @ upgrade-start (guards the AST parser itself:
    # a parser that silently found nothing would pass every other test).
    # Entries counted in the dict literal, duplicates included; the RTI
    # table had 99 before its duplicate `alert_rti` was removed (B23).
    # bokeh 3.1 added `start_plot_viewer` to both tables (+1 each).
    counts = {key: len(_entries(*PIPELINES[key])) for key in PIPELINES}
    assert counts == {"kcwi_pipeline": 98, "keck_rti_pipeline": 99}, counts
