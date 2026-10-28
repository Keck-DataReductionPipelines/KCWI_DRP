"""bless needs tiers 0 and 1 to have status exactly "pass" (T.7 follow-up).

A tier that did not run (deselected with ``-k``), was made informational, or
was skipped does not vouch for the run, so the run must not become a
reference. Written from the rule in docs/testing.md section 10.
"""
import pytest

from kcwidrp.tests.regression import bless as blessmod
from kcwidrp.tests.regression import compare
from kcwidrp.tests.regression import manifest as mf
from kcwidrp.tests.regression.report import SessionResults
from kcwidrp.tests.regression.session import EntryRun

EID = "ds9_blue_L1"


def run_session(w, tmp_path, tiers):
    """A kept run of ds9_blue_L1 where only ``tiers`` were evaluated."""
    session_dir = tmp_path / "s"
    results = SessionResults(session_dir)
    entry = mf.get_entry(w.mfst, EID)
    run = EntryRun(w.data, w.mfst, entry, results, session_dir,
                   command=w.command, keep_runs=True).start()
    for t in tiers:
        run.tier(t)
    run.finish()
    return session_dir


def try_bless(w, session_dir):
    return blessmod.bless(w.data, EID, "initial: x", run=session_dir,
                          assume_yes=True)


def assert_refused_and_unchanged(w):
    assert not (w.data / "reference" / EID).exists()
    assert mf.get_entry(mf.load_manifest(w.data), EID)["bless_history"] == []


def test_refuses_when_tiers_0_and_1_never_ran(make_world, tmp_path):
    w = make_world()
    session_dir = run_session(w, tmp_path, range(2, 8))
    with pytest.raises(blessmod.BlessError, match="did not run") as exc:
        try_bless(w, session_dir)
    assert "tier 0" in str(exc.value) or "tier 1" in str(exc.value)
    assert_refused_and_unchanged(w)


def test_refuses_when_tier_1_never_ran(make_world, tmp_path):
    w = make_world()
    session_dir = run_session(w, tmp_path, [0])
    assert SessionResults.load(session_dir).tier_status(EID, 0) == "pass"
    with pytest.raises(blessmod.BlessError, match="did not run") as exc:
        try_bless(w, session_dir)
    assert "tier 1" in str(exc.value)
    assert_refused_and_unchanged(w)


def test_refuses_when_tier_0_never_ran(make_world, tmp_path):
    w = make_world()
    session_dir = run_session(w, tmp_path, range(1, 8))
    with pytest.raises(blessmod.BlessError, match="did not run") as exc:
        try_bless(w, session_dir)
    assert "tier 0" in str(exc.value)
    assert_refused_and_unchanged(w)


@pytest.mark.parametrize("status", ["info", "skipped"])
def test_refuses_when_tier_1_is_not_pass(make_world, tmp_path, status):
    w = make_world()
    session_dir = run_session(w, tmp_path, range(8))
    # Rewrite tier 1 in results.json, as an informational ('against') entry or
    # a skipped tier would have left it.
    results = SessionResults.load(session_dir)
    results.set_tier(EID, compare.TierResult(1, status, "forced"))
    results.write()
    with pytest.raises(blessmod.BlessError, match="not 'pass'") as exc:
        try_bless(w, session_dir)
    assert "tier 1" in str(exc.value)
    assert_refused_and_unchanged(w)


def test_blesses_when_every_tier_ran_and_passed(make_world, tmp_path):
    w = make_world()
    session_dir = run_session(w, tmp_path, range(8))
    res = SessionResults.load(session_dir)
    assert [res.tier_status(EID, t) for t in (0, 1)] == ["pass", "pass"]
    try_bless(w, session_dir)
    assert (w.data / "reference" / EID / "redux" / "ds_icubes.fits").is_file()
    hist = mf.get_entry(mf.load_manifest(w.data), EID)["bless_history"]
    assert [h["reason"] for h in hist] == ["initial: x"]
