"""reduce_kcwi / kcwi_rti argument handling (plan T.6).

Seams used (no part of the reduction is run):

* Argument parsing: both scripts have a module-level
  ``_parse_arguments(in_args)`` that takes the full argv list (``in_args[0]``
  is the program name) and returns the ``argparse.Namespace``; the tests call
  it directly.
* Bug B8 (channel check): the decision is written inline in ``main()``, as
  ``if args.blue and <condition on frame>:`` (and the same for ``args.red``),
  in the ``-f`` and ``-l`` branches. The tests parse the script source,
  extract the condition expression from each such ``if`` and evaluate it with
  ``eval`` on a sample path bound to ``frame``. The code under test is the
  script's own expression, not a copy.
* Bug B11 (RTI ignores Q): in the ``-f`` branch the tests parse ``main()``
  and check which name is passed to ``framework.ingest_data``. (The branch
  cannot be run without a framework, a config and ``input()``.)
"""
import ast
import importlib
import importlib.util
import warnings

import pytest

SCRIPTS = ["reduce_kcwi", "kcwi_rti"]


def _module(script):
    return importlib.import_module(f"kcwidrp.scripts.{script}")


def _parse(script, *argv):
    extra = ["--rti-ingesttype", "lev1"] if script == "kcwi_rti" else []
    return _module(script)._parse_arguments([script, *extra, *argv])


# --- documented behaviour of -f -l -b -r -------------------------------------

@pytest.mark.parametrize("script", SCRIPTS)
def test_frames_option_takes_a_list(script):
    # -f / --frames: 'input image files (full path, list ok)' (help text)
    # source: specification (option help text in the script)
    assert _parse(script, "-f", "a.fits", "b.fits").frames == [
        "a.fits", "b.fits"]
    assert _parse(script, "--frames", "/d/kb1.fits").frames == ["/d/kb1.fits"]
    assert _parse(script, "-f").frames == []
    assert _parse(script).frames is None


@pytest.mark.parametrize("script", SCRIPTS)
def test_list_option_takes_one_file_name(script):
    # -l / --list: 'File containing a list of files to be processed'
    # source: specification (option help text in the script)
    assert _parse(script, "-l", "bias.txt").file_list == "bias.txt"
    assert _parse(script, "--list", "x.lst").file_list == "x.lst"
    assert _parse(script).file_list is None


@pytest.mark.parametrize("script", SCRIPTS)
def test_channel_flags(script):
    # -b / --blue: 'KCWI Blue processing'; -r / --red: 'KCWI Red processing'
    # source: specification (option help text in the script)
    args = _parse(script)
    assert (args.blue, args.red) == (False, False)
    args = _parse(script, "-b")
    assert (args.blue, args.red) == (True, False)
    args = _parse(script, "--blue")
    assert (args.blue, args.red) == (True, False)
    args = _parse(script, "-r")
    assert (args.blue, args.red) == (False, True)
    args = _parse(script, "--red")
    assert (args.blue, args.red) == (False, True)


@pytest.mark.parametrize("script", SCRIPTS)
def test_options_combine(script):
    args = _parse(script, "-b", "-f", "x.fits", "y.fits", "-k")
    # source: specification (option help text in the script)
    assert args.blue and not args.red and args.skipsky
    assert args.frames == ["x.fits", "y.fits"]
    args = _parse(script, "-r", "-l", "list.txt")
    assert args.red and args.file_list == "list.txt" and args.frames is None


def test_rti_requires_ingest_type():
    # kcwi_rti: --rti-ingesttype is required, choices lev1 / lev2
    mod = _module("kcwi_rti")
    with pytest.raises(SystemExit) as exc:
        mod._parse_arguments(["kcwi_rti", "-b"])
    # source: argparse exit status for a usage error
    assert exc.value.code == 2
    with pytest.raises(SystemExit):
        mod._parse_arguments(["kcwi_rti", "--rti-ingesttype", "lev3"])
    assert mod._parse_arguments(
        ["kcwi_rti", "--rti-ingesttype", "lev2"]).rti_ingesttype == "lev2"


# --- AST seams -------------------------------------------------------------------

def _main_tree(script):
    spec = importlib.util.find_spec(f"kcwidrp.scripts.{script}")
    with open(spec.origin) as f, warnings.catch_warnings():
        warnings.simplefilter("ignore", SyntaxWarning)
        tree = ast.parse(f.read(), filename=spec.origin)
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == "main":
            return spec.origin, node
    raise AssertionError(f"{spec.origin}: no main()")


def _is_args_attr(node, attr):
    return (isinstance(node, ast.Attribute) and node.attr == attr
            and isinstance(node.value, ast.Name) and node.value.id == "args")


def _channel_predicates(script, channel):
    """Condition expressions of `if args.<channel> and <cond on frame>:`."""
    _, main = _main_tree(script)
    found = []
    for node in ast.walk(main):
        if (isinstance(node, ast.If) and isinstance(node.test, ast.BoolOp)
                and isinstance(node.test.op, ast.And)
                and len(node.test.values) == 2
                and _is_args_attr(node.test.values[0], channel)):
            expr = ast.Expression(node.test.values[1])
            ast.fix_missing_locations(expr)
            found.append((node.lineno, compile(expr, f"<{script}:{node.lineno}>",
                                               "eval")))
    return found


def _wrong_channel(script, channel, path):
    """Linenos at which the script's own condition flags `path`."""
    return [line for line, code in _channel_predicates(script, channel)
            if eval(code, {"frame": path})]


@pytest.mark.parametrize("script", SCRIPTS)
@pytest.mark.parametrize("channel, wrong_file", [
    ("blue", "/data/obs/kr260115_00012.fits"),
    ("red", "/data/obs/kb260115_00012.fits"),
])
def test_channel_check_flags_a_wrong_channel_file(script, channel, wrong_file):
    # Positive control for the extraction: a red frame with -b (or a blue
    # frame with -r) triggers the 'wrong channel' prompt in both the -f and
    # the -l branch. KCWI frames are kb<date>_<n>.fits (blue), kr... (red).
    preds = _channel_predicates(script, channel)
    # source: reduce_kcwi.py / kcwi_rti.py main(): one check per branch
    assert len(preds) == 2, preds
    assert len(_wrong_channel(script, channel, wrong_file)) == 2


@pytest.mark.parametrize("script", SCRIPTS)
@pytest.mark.parametrize("channel, right_file", [
    ("blue", "/data/obs/kb260115_00012.fits"),
    ("red", "/data/obs/kr260115_00012.fits"),
])
def test_channel_check_accepts_matching_files(script, channel, right_file):
    # source: specification (plan T.6 / bug B8 description)
    assert _wrong_channel(script, channel, right_file) == []


@pytest.mark.xfail(strict=True, reason=(
    "B8: the channel check looks for 'kr'/'kb' in the whole path, so a "
    "directory name containing it triggers a false wrong-channel prompt "
    "(reduce_kcwi.py:436-469, kcwi_rti.py:412-440); fixed in PR3b"))
@pytest.mark.parametrize("script", SCRIPTS)
@pytest.mark.parametrize("channel, path", [
    # blue frame in a directory whose name contains 'kr'
    ("blue", "/data/kr_nights/2026-01-15/kb260115_00012.fits"),
    # red frame in a directory whose name contains 'kb'
    ("red", "/data/kbnights/2026-01-15/kr260115_00012.fits"),
])
def test_channel_check_ignores_directory_names(script, channel, path):
    # source: specification (plan section 13, B8): only the file name tells
    # the channel; the frame has the requested channel
    assert _wrong_channel(script, channel, path) == []


def _ingest_frames_arg(script):
    """The second argument of framework.ingest_data(...) in the -f branch."""
    origin, main = _main_tree(script)
    for node in ast.walk(main):
        if isinstance(node, ast.If) and _is_args_attr(node.test, "frames"):
            for sub in ast.walk(ast.Module(body=node.body, type_ignores=[])):
                if (isinstance(sub, ast.Call)
                        and isinstance(sub.func, ast.Attribute)
                        and sub.func.attr == "ingest_data"):
                    return origin, sub.lineno, ast.unparse(sub.args[1])
    raise AssertionError(f"{origin}: no ingest_data in the -f branch")


@pytest.mark.parametrize("script", [
    "reduce_kcwi",
    pytest.param("kcwi_rti", marks=pytest.mark.xfail(strict=True, reason=(
        "B11: in the -f branch of kcwi_rti.py, answering Q at the "
        "wrong-channel prompt empties the local 'frames' list, but ingest_data is then "
        "called with args.frames, so the quit is ignored "
        "(kcwi_rti.py:405-423); fixed in PR3b"))),
])
def test_f_branch_ingests_the_filtered_list(script):
    origin, line, arg = _ingest_frames_arg(script)
    # source: specification (plan section 13, B11): reduce_kcwi.py ingests
    # the filtered `frames` list; the RTI script must do the same
    assert arg == "frames", (
        f"{origin}:{line}: ingest_data is passed {arg!r}, not the filtered "
        f"local list 'frames'")
