"""B9 (plan T.11): each launch script initialises the plot session in time.

The scripts cannot be imported (smart_reduce_kcwi fails to import today and
main() cannot be run in pieces), so the source is parsed with ast.

Rule checked, per script: in each function that creates the framework
(``<x> = Framework(...)``), a call ``init_plot_session(<x>.context)`` (bare
name or attribute, e.g. ``bokeh_plotting.init_plot_session``) appears after
the framework creation and before the earliest ``<x>.append_event(...)`` or
``<x>.ingest_data(...)`` call, by line number. Only the function's own body is
searched: nested function definitions (``process``, ``process_subset``,
``process_list``) are skipped, because they are called later and their
``append_event`` lines precede the framework creation textually.
"""
import ast
import importlib.util

import pytest

SCRIPTS = [
    "kcwidrp.scripts.reduce_kcwi",
    "kcwidrp.scripts.kcwi_rti",
    "kcwidrp.scripts.smart_reduce_kcwi",
]
QUEUEING = ("append_event", "ingest_data")


def _own_nodes(func):
    """Yield nodes in func's body, not descending into nested defs."""
    stack = [n for n in func.body
             if not isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]
    while stack:
        node = stack.pop()
        yield node
        for child in ast.iter_child_nodes(node):
            if not isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef,
                                      ast.Lambda)):
                stack.append(child)


def _is_name(node, name):
    return isinstance(node, ast.Name) and node.id == name


def _called_name(call):
    f = call.func
    if isinstance(f, ast.Name):
        return f.id
    if isinstance(f, ast.Attribute):
        return f.attr
    return None


@pytest.mark.parametrize("module", SCRIPTS)
def test_script_calls_init_plot_session_before_queueing(module):
    spec = importlib.util.find_spec(module)
    assert spec is not None and spec.origin, f"cannot locate {module}"
    path = spec.origin
    with open(path) as f:
        tree = ast.parse(f.read(), filename=path)

    checked = 0
    for func in ast.walk(tree):
        if not isinstance(func, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        nodes = list(_own_nodes(func))
        # framework creation: <x> = Framework(...)
        created = []  # (lineno, name)
        for n in nodes:
            if (isinstance(n, ast.Assign) and isinstance(n.value, ast.Call)
                    and _called_name(n.value) == "Framework"):
                for t in n.targets:
                    if isinstance(t, ast.Name):
                        created.append((n.lineno, t.id))
        for create_line, name in created:
            checked += 1
            where = f"{path}: {func.name}()"
            inits = [
                n.lineno for n in nodes
                if isinstance(n, ast.Call)
                and _called_name(n) == "init_plot_session"
                and len(n.args) == 1
                and isinstance(n.args[0], ast.Attribute)
                and n.args[0].attr == "context"
                and _is_name(n.args[0].value, name)
            ]
            queued = [
                n.lineno for n in nodes
                if isinstance(n, ast.Call)
                and isinstance(n.func, ast.Attribute)
                and n.func.attr in QUEUEING
                and _is_name(n.func.value, name)
            ]
            # source: specification (plan T.11 / §13 B9, B10)
            assert inits, (
                f"{where}: no call init_plot_session({name}.context) found; "
                f"it must follow the framework creation (line {create_line})"
                f" and precede the first {name}.append_event/"
                f"{name}.ingest_data call")
            assert queued, (
                f"{where}: found no {name}.append_event/{name}.ingest_data "
                f"call, so the ordering rule cannot be checked")
            first_queue = min(queued)
            # source: specification (plan T.11 / §13 B9, B10)
            assert any(create_line < i < first_queue for i in inits), (
                f"{where}: init_plot_session({name}.context) at line(s) "
                f"{inits} is not between the framework creation (line "
                f"{create_line}) and the first {name}.append_event/"
                f"{name}.ingest_data call (line {first_queue})")

    # source: specification (plan T.11 / §13 B9, B10)
    assert checked >= 1, f"{path}: no function creating a Framework found"
