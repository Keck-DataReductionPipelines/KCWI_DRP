"""The regression test module run inside pytester (T.7): skip reason, ids,
ordering of entry runs and teardown.

The real ``test_regression.py`` and ``regression/conftest.py`` are exercised
with a fake reduction substituted for ``runner.default_command``.
"""
from pathlib import Path

import pytest

from kcwidrp.tests.regression import compare, runner
from kcwidrp.tests.regression import manifest as mf
from kcwidrp.tests.regression.session import EntryRun  # noqa: F401

pytest_plugins = ["pytester"]

REG = Path(__file__).resolve().parents[2] / "regression"

CONFTEST = '''
import pytest
from kcwidrp.tests.regression.options import add_options as pytest_addoption
from kcwidrp.tests.regression.conftest import *  # noqa: F401,F403

def pytest_configure(config):
    config.addinivalue_line("markers", "regression: needs data")
'''


def setup_pytester(pytester):
    pytester.makeconftest(CONFTEST)
    pytester.makepyfile(test_regression=(REG / "test_regression.py").read_text())


def test_all_regression_tests_skip_when_the_data_dir_is_unset(pytester,
                                                              monkeypatch):
    monkeypatch.delenv("KCWI_REGRESSION_DATA", raising=False)
    setup_pytester(pytester)
    result = pytester.runpytest_inprocess("-rs", "-p", "no:cacheprovider")
    result.assert_outcomes(skipped=8)           # one placeholder x 8 tiers
    result.stdout.fnmatch_lines(["*KCWI_REGRESSION_DATA not set*"])


# The session below runs three fake reductions in subprocesses (about 0.4 s
# each), so it cannot meet the 2 s limit.
@pytest.mark.slow
def test_entries_run_one_at_a_time_ids_and_teardown(pytester, make_world,
                                                    command_factory,
                                                    monkeypatch, tmp_path):
    w = make_world()
    # references: bless the L1 entry's first run by hand (copy the outputs)
    import shutil
    mfst = mf.load_manifest(w.data)
    cmd = command_factory(w.script, fail_ids=("ds9_blue_L0",))
    for eid in ("ds9_blue_L0", "ds9_blue_L1"):
        entry = mf.get_entry(mfst, eid)
        d = tmp_path / "ref" / eid
        runner.run_entry(w.data, mfst, entry, d, command=cmd)
        ref = w.data / "reference" / eid
        ref.mkdir(parents=True)
        shutil.copytree(d / "redux", ref / "redux")
        shutil.copytree(d / "plots", ref / "plots")
        shutil.copy(d / "kcwib.proc", ref)
    # make the third entry (ds9_blue_L0) fail tier 1 at the session run
    live_at_launch = {}
    out = tmp_path / "out"

    def fake(entry, run_dir):
        sess = Path(run_dir).parent
        live_at_launch[entry["id"]] = sorted(
            p.name for p in sess.iterdir() if (p / "redux").is_dir())
        return cmd(entry, run_dir)

    monkeypatch.setattr(runner, "default_command", fake)
    monkeypatch.setenv("KCWI_REGRESSION_DATA", str(w.data))
    monkeypatch.setenv("KCWI_REGRESSION_OUT", str(out))
    setup_pytester(pytester)
    result = pytester.runpytest_inprocess(
        "--regression-profile=full", "-p", "no:cacheprovider", "-v", "-m",
        "regression")
    names = [l.split("::")[-1].split()[0] for l in result.outlines
             if l.startswith("test_regression.py::")]
    # ids are <entry>-<tier>_<name>, so -k can select entry or tier
    assert "test_regression[ds9_blue_L1-tier2_structure]" in names
    # all tests of an entry are consecutive: the entries appear in blocks
    entries = [n.split("[")[1].split("-")[0] for n in names]
    blocks = [e for i, e in enumerate(entries) if i == 0 or entries[i - 1] != e]
    assert sorted(blocks) == sorted(set(entries)) and len(blocks) == 3
    session = next(out.iterdir())
    # L1 and L3 entries passed or were skipped; the failing entry kept redux.
    # Order of launch is the manifest order: L0 (fails), L1 (passes), L3.
    assert live_at_launch["ds9_blue_L1"] == ["ds9_blue_L0"]   # failed one kept
    assert live_at_launch["ds9_blue_L3"] == ["ds9_blue_L0"]   # L1 torn down
    assert (session / "ds9_blue_L0" / "redux").is_dir()
    assert not (session / "ds9_blue_L1" / "redux").exists()
    for kept in ("logs", "prompt_record.jsonl", "reduce.out", "kcwi.cfg",
                 "files.txt", "entry_report.md"):
        if kept == "prompt_record.jsonl" and not (
                session / "ds9_blue_L1" / kept).exists():
            continue
        assert (session / "ds9_blue_L1" / kept).exists(), kept
    assert (session / "report.md").is_file()
    assert (session / "results.json").is_file()
    text = (session / "report.md").read_text()
    assert "| ds9_blue_L0 |" in text and "**FAIL**" in text
