"""Entry points (plan T.6): every ``[tool.poetry.scripts]`` target of
pyproject.toml imports, the function exists, and ``<name> --help`` exits 0.

``--help`` is run in a subprocess as
``python -c "from <module> import <func>; sys.argv = [<name>, '--help']; <func>()"``
so that it works from the source tree without an installed console script.
Each subprocess has ``DISPLAY`` unset, runs in an empty temporary directory
and must leave it empty: ``--help`` of the launch scripts exits inside
``argparse`` before the script creates ``logs/``, ``plots/``, a config file or
starts a bokeh server (``subprocess.Popen('bokeh serve')``); if one of them
ever did, the empty-directory check fails. The subprocesses run concurrently
so the file takes about as long as the slowest import (bokeh, ~1 s).

``wb`` and ``wr`` take file names only: they have no ``--help`` and treat the
word as a FITS file name (B24).
"""
import importlib
import os
import subprocess
import sys
import tomllib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]

# Entry points that fail today, with the plan ID that fixes them (same IDs as
# tests/unit/test_imports.py). strict=True: remove the entry with the fix.
IMPORT_FAILURES = {
    "kcwi_masksky_ds9": "B14: kcwidrp.scripts.kcwi_masksky_ds9 needs the "
        "undeclared package 'regions'; pyproject still declares pyregion",
}

# Entry points whose --help fails although they import.
HELP_FAILURES = {
    "wb": "B24: wb treats '--help' as a FITS file name; "
          "get_log_string returns 2 values on a bad file where wb_main "
          "unpacks 3 (wb.py:167), so it dies with ValueError",
    "wr": "B24: wr treats '--help' as a FITS file name; "
          "get_log_string returns 2 values on a bad file where wr_main "
          "unpacks 3 (wr.py:165), so it dies with ValueError",
}


def _entry_points():
    with open(REPO_ROOT / "pyproject.toml", "rb") as f:
        scripts = tomllib.load(f)["tool"]["poetry"]["scripts"]
    out = {}
    for name, target in scripts.items():
        module, func = target.split(":")
        out[name] = (module, func)
    return out


ENTRY_POINTS = _entry_points()


def _run_help(name, module, func, cwd):
    env = {k: v for k, v in os.environ.items() if k != "DISPLAY"}
    env["PYTHONPATH"] = str(REPO_ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["PYTHONWARNINGS"] = "ignore"
    code = (f"import sys; from {module} import {func} as f; "
            f"sys.argv = [{name!r}, '--help']; f()")
    try:
        proc = subprocess.run([sys.executable, "-c", code], cwd=cwd, env=env,
                              capture_output=True, text=True, timeout=60,
                              stdin=subprocess.DEVNULL)
    except subprocess.TimeoutExpired as exc:
        return None, f"timeout after {exc.timeout} s", cwd
    return proc.returncode, proc.stdout + proc.stderr, cwd


@pytest.fixture(scope="module")
def help_results(tmp_path_factory):
    jobs = {}
    with ThreadPoolExecutor(max_workers=len(ENTRY_POINTS)) as pool:
        for name, (module, func) in ENTRY_POINTS.items():
            cwd = tmp_path_factory.mktemp(f"help_{name}")
            jobs[name] = pool.submit(_run_help, name, module, func, cwd)
    return {name: job.result() for name, job in jobs.items()}


def _params(failures):
    params = []
    for name in ENTRY_POINTS:
        marks = []
        if name in failures:
            marks.append(pytest.mark.xfail(strict=True,
                                           reason=failures[name]))
        params.append(pytest.param(name, id=name, marks=marks))
    return params


def test_expected_entry_points_are_declared():
    # source: pyproject.toml [tool.poetry.scripts] at upgrade-start; guards
    # the parser and notices a new script that needs a case above
    assert sorted(ENTRY_POINTS) == sorted([
        "reduce_kcwi", "kcwi_masksky_ds9", "smart_reduce_kcwi", "modhead",
        "start_kcwi_rti", "wb", "wr", "check_cals"])


@pytest.mark.parametrize("name", _params(IMPORT_FAILURES))
def test_entry_point_imports_and_function_exists(name):
    module_name, func = ENTRY_POINTS[name]
    module = importlib.import_module(module_name)
    # source: specification (plan T.6)
    assert callable(getattr(module, func, None)), (
        f"{module_name} has no callable {func}")


@pytest.mark.parametrize("name", _params({**IMPORT_FAILURES,
                                          **HELP_FAILURES}))
def test_entry_point_help_exits_zero(name, help_results):
    code, output, cwd = help_results[name]
    # source: specification (plan T.6): --help exits 0
    assert code == 0, f"{name} --help exited {code}:\n{output[-1500:]}"
    assert "usage" in output.lower(), output[-500:]
    # --help must not create files (logs/, plots/, config, bokeh output) nor
    # start a server
    assert sorted(p.name for p in cwd.iterdir()) == []
