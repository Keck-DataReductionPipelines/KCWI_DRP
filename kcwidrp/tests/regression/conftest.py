"""Regression-suite hooks and fixtures (T.1 skip hook, T.7 runner).

* Every test marked ``regression`` is skipped when ``$KCWI_REGRESSION_DATA``
  is unset (reason below).
* ``test_regression`` is parametrized by manifest entry (session-scoped
  fixture ``entry_run``) and by tier, ids like ``ds01_blue_L1-tier2_structure``.
  pytest groups the tests of one entry, so only one entry's run is alive.
* The options ``--regression-profile`` and ``--keep-runs`` are registered in
  the root conftest (options cannot be added from a subdirectory).
"""
import os
import sys
from pathlib import Path

import pytest

from kcwidrp.tests.regression import compare, manifest as mf, runner
from kcwidrp.tests.regression.report import SessionResults
from kcwidrp.tests.regression.session import EntryRun

SKIP_REASON = "KCWI_REGRESSION_DATA not set"
UNCONFIGURED = "unconfigured"
MANIFEST_ERROR = "manifest-error"


def pytest_collection_modifyitems(config, items):
    """Skip every test marked ``regression`` when the data directory is unset."""
    if os.environ.get("KCWI_REGRESSION_DATA"):
        return
    skip = pytest.mark.skip(reason=SKIP_REASON)
    for item in items:
        if item.get_closest_marker("regression") is not None:
            item.add_marker(skip)


def _entry_ids(config):
    """(ids, error message or None) for the selected profile."""
    if not os.environ.get("KCWI_REGRESSION_DATA"):
        return [UNCONFIGURED], None
    try:
        mfst = mf.load_manifest(os.environ["KCWI_REGRESSION_DATA"])
    except mf.ManifestError as exc:
        return [MANIFEST_ERROR], str(exc)
    profile = config.getoption("--regression-profile")
    return [e["id"] for e in mf.select_entries(mfst, profile)], None


def pytest_generate_tests(metafunc):
    if "entry_run" in metafunc.fixturenames:
        ids, _err = _entry_ids(metafunc.config)
        metafunc.parametrize("entry_run", ids, indirect=True, scope="session")
    if "tier" in metafunc.fixturenames:
        tiers = sorted(compare.TIER_NAMES)
        metafunc.parametrize(
            "tier", tiers,
            ids=["tier%d_%s" % (t, compare.TIER_NAMES[t]) for t in tiers])


@pytest.fixture(scope="session")
def regression_session(request):
    """Session directory, results.json and report.md (written as it goes)."""
    config = request.config
    data_dir = Path(os.environ["KCWI_REGRESSION_DATA"])
    session_dir = runner.new_session_dir()
    meta = {"profile": config.getoption("--regression-profile"),
            "data_dir": str(data_dir), "drp": runner.git_describe(),
            "python": sys.version.split()[0]}
    results = SessionResults(session_dir, meta)
    results.write()
    config._kcwi_regression_session = session_dir
    yield results
    results.write()


@pytest.fixture(scope="session")
def entry_run(request, regression_session):
    """Reduce one entry (checksums first), tear it down after its last test."""
    entry_id = request.param
    if entry_id == MANIFEST_ERROR:
        pytest.fail(_entry_ids(request.config)[1], pytrace=False)
    data_dir = Path(os.environ["KCWI_REGRESSION_DATA"])
    mfst = mf.load_manifest(data_dir)
    entry = mf.get_entry(mfst, entry_id)
    run = EntryRun(
        data_dir, mfst, entry, regression_session,
        regression_session.session_dir,
        timeout=int(os.environ.get("KCWI_REGRESSION_TIMEOUT",
                                   runner.DEFAULT_TIMEOUT)),
        keep_runs=request.config.getoption("--keep-runs"))
    run.start()
    yield run
    run.finish()


def pytest_terminal_summary(terminalreporter, config):
    session = getattr(config, "_kcwi_regression_session", None)
    if session is not None:
        terminalreporter.write_line("regression report: %s" %
                                    (Path(session) / "report.md"))
