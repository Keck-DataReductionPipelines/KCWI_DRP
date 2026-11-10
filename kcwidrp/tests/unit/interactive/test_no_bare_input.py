"""AST scan: no bare ``input(...)`` and no ``time.sleep(...plot_pause...)``
in the primitives and launch scripts (plan P.3 / P.4 validation).

Specification (docs/drp_upgrade_plan.md section 8, P.3 and P.4 validation):
every prompt goes through ``kcwidrp/core/kcwi_interactive.py``.  The scan
covers ``kcwidrp/primitives/*.py``, ``kcwidrp/scripts/reduce_kcwi.py`` and
``kcwidrp/scripts/kcwi_rti.py``; the helper module itself is not scanned.

* ``input(...)`` (or ``builtins.input(...)``) is never allowed: the
  allow-list has no entry for it.
* ``time.sleep(<expression mentioning plot_pause>)`` is a ``pause`` that was
  not converted.  A few call-site sleeps legitimately stay, because the
  helper's ``pause`` sleeps one ``plot_pause`` and these need a different
  duration or sit outside a prompt.  They are listed in ``SLEEP_ALLOWED``
  as ``{(file, function_name): reason}``.

Before P.3/P.4 are done this test fails, naming every bare site: expected.
"""
import ast
import pathlib
import warnings

REPO = pathlib.Path(__file__).resolve().parents[4]
PKG = REPO / "kcwidrp"

SCANNED = sorted((PKG / "primitives").glob("*.py")) + [
    PKG / "scripts" / "reduce_kcwi.py", PKG / "scripts" / "kcwi_rti.py"]

# (file relative to the repository root, enclosing function) -> reason.
# Sleeps only; there must never be an entry for input(.
SLEEP_ALLOWED = {
    ("kcwidrp/primitives/MakeInvsens.py", "_perform"):
        "the plot_level 1 branches after the efficiency, inverse "
        "sensitivity, calibrated and residual figures sleep "
        "2. * plot_pause, which pause() (one plot_pause) does not offer",
    ("kcwidrp/primitives/SubtractOverscan.py", "_perform"):
        "the level-2 sleep: the per-amplifier and gridplot figures sleep "
        "plot_pause in a branch that only prompts at plot_level >= 3",
    ("kcwidrp/primitives/MakeCube.py", "_perform"):
        "after the user quits the slice review the remaining slices are "
        "still shown for plot_pause each",
}


def _owners(tree):
    """Map each node to the name of its innermost enclosing function."""
    owner = {}

    def visit(node, name):
        for child in ast.iter_child_nodes(node):
            owner[child] = name
            visit(child, child.name if isinstance(
                child, (ast.FunctionDef, ast.AsyncFunctionDef)) else name)

    visit(tree, "<module>")
    return owner


def _is_input_call(call):
    f = call.func
    if isinstance(f, ast.Name):
        return f.id == "input"
    return (isinstance(f, ast.Attribute) and f.attr == "input"
            and isinstance(f.value, ast.Name) and f.value.id == "builtins")


def _is_plot_pause_sleep(call):
    f = call.func
    named = (isinstance(f, ast.Name) and f.id == "sleep") or (
        isinstance(f, ast.Attribute) and f.attr == "sleep"
        and isinstance(f.value, ast.Name) and f.value.id == "time")
    args = list(call.args) + [k.value for k in call.keywords]
    return named and any("plot_pause" in ast.unparse(a) for a in args)


def violations(source, rel, allowed=()):
    """List "file:line: ..." for every violation in ``source``."""
    with warnings.catch_warnings():
        # docstrings in the primitives hold "\*" (invalid escape): not ours
        warnings.simplefilter("ignore", SyntaxWarning)
        tree = ast.parse(source, filename=rel)
    owner = _owners(tree)
    found = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = owner.get(node)
        if _is_input_call(node):
            found.append("%s:%d: bare input() in %s" % (rel, node.lineno,
                                                       func))
        elif _is_plot_pause_sleep(node) and (rel, func) not in allowed:
            found.append("%s:%d: time.sleep(...plot_pause...) in %s"
                         % (rel, node.lineno, func))
    return sorted(found, key=lambda s: int(s.split(":")[1]))   # one file


def test_scanner_flags_each_kind_and_honors_the_allow_list():
    src = ("import time, builtins\n"
           "def f(self):\n"
           "    input('a')\n"
           "    builtins.input('b')\n"
           "    time.sleep(2. * self.config.instrument.plot_pause)\n"
           "    time.sleep(1)\n"
           "    x = 'input(' + 'plot_pause'\n")
    # source: analytic - lines 3, 4 and 5 are the only violations
    assert [v.split(":")[1] for v in violations(src, "m.py")] == \
        ["3", "4", "5"]
    assert [v.split(":")[1]
            for v in violations(src, "m.py", {("m.py", "f")})] == ["3", "4"]


def test_allow_list_is_sleeps_only_and_explained():
    assert len(SLEEP_ALLOWED) == 3
    assert all(isinstance(k, tuple) and len(k) == 2 and v.strip()
               for k, v in SLEEP_ALLOWED.items())


def test_scanned_files_exist():
    assert len(SCANNED) > 30
    assert all(p.exists() for p in SCANNED)


def test_no_bare_input_or_plot_pause_sleep():
    found = []
    for path in SCANNED:
        found += violations(path.read_text(),
                            path.relative_to(REPO).as_posix(),
                            SLEEP_ALLOWED)
    assert not found, ("%d bare prompt site(s):\n  " % len(found)
                       + "\n  ".join(found))
