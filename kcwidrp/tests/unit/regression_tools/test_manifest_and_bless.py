"""draft-manifest, manifest validation and the bless workflow (T.7)."""
import hashlib
import json

import pytest

from kcwidrp.tests.regression import __main__ as cli
from kcwidrp.tests.regression import bless as blessmod
from kcwidrp.tests.regression import manifest as mf
from kcwidrp.tests.regression.report import SessionResults
from kcwidrp.tests.regression.session import EntryRun


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_draft_manifest_entries_and_checksums(make_world):
    w = make_world()
    manifest, notes = mf.draft_manifest(w.data)
    ds = manifest["datasets"]["ds9_fake"]
    lev0 = w.data / "raw" / "ds9_fake" / "lev0"
    # source: independent -- hashlib over the raw bytes
    assert ds["files"] == {"lev0/KB.1.fits": sha(lev0 / "KB.1.fits"),
                           "lev0/KB.2.fits": sha(lev0 / "KB.2.fits")}
    ids = [e["id"] for e in manifest["entries"]]
    assert ids == ["ds9_blue_L0", "ds9_blue_L1", "ds9_blue_L3"]
    by = {e["id"]: e for e in manifest["entries"]}
    assert [by[i]["profile"] for i in ids] == ["full", "quick", "full"]
    assert by["ds9_blue_L3"]["answers"] == "answers_ds9_fake_blue_3.txt"
    assert by["ds9_blue_L3"]["answers_sha256"] == sha(
        w.data / "answers" / "answers_ds9_fake_blue_3.txt")
    assert by["ds9_blue_L1"]["config_overrides"] == {
        "enable_bokeh": False, "plot_pause": 0, "plot_level": 1}
    assert by["ds9_blue_L1"]["sky_file"] == "kcwi.sky"
    assert any("no reference/ds9_blue_L1" in n for n in notes)
    assert mf.validate_manifest(manifest) == []


def test_draft_reads_provenance_from_reference_history(make_world):
    from astropy.io import fits
    w = make_world()
    redux = w.data / "reference" / "ds9_blue_L1" / "redux"
    redux.mkdir(parents=True)
    h = fits.PrimaryHDU()
    h.header.add_history("kcwidrp version=9.8.7")
    h.header.add_history("git version=v9-1-gabc")
    h.writeto(redux / "x_icubes.fits")
    manifest, _ = mf.draft_manifest(w.data)
    prov = {e["id"]: e["provenance"] for e in manifest["entries"]}
    assert prov["ds9_blue_L1"]["drp_version"] == "9.8.7"
    assert prov["ds9_blue_L1"]["drp_git_version"] == "v9-1-gabc"
    assert prov["ds9_blue_L0"]["drp_version"] is None


def test_validation_requires_a_reason_for_every_tolerance(make_world):
    m = make_world().mfst
    m["entries"][1]["tolerances"] = {"*_icubes.fits": {"rtol": 1e-9}}
    errors = mf.validate_manifest(m)
    assert len(errors) == 1 and "reason" in errors[0]
    m["entries"][1]["tolerances"] = {}
    m["entries"][2]["answers"] = None
    assert "answers" in mf.validate_manifest(m)[0]


def test_draft_cli_refuses_to_overwrite(make_world, capsys):
    w = make_world()
    assert cli.main(["--data-dir", str(w.data), "draft-manifest"]) == 2
    assert cli.main(["--data-dir", str(w.data), "draft-manifest",
                     "--force"]) == 0
    assert "3 entries, 2 checksums" in capsys.readouterr().out


def test_verify_checksums_cli(make_world, capsys):
    w = make_world()
    assert cli.main(["--data-dir", str(w.data), "verify-checksums"]) == 0
    (w.data / "raw" / "ds9_fake" / "lev0" / "KB.1.fits").write_bytes(b"x")
    assert cli.main(["--data-dir", str(w.data), "verify-checksums"]) == 1


# -- bless ----------------------------------------------------------------------

@pytest.fixture
def first_run(make_world, tmp_path):
    """A kept run of an entry that has no reference yet, in a session."""
    w = make_world()
    session_dir = tmp_path / "out" / "20260101-000000"
    results = SessionResults(session_dir)
    entry = mf.get_entry(w.mfst, "ds9_blue_L1")
    run = EntryRun(w.data, w.mfst, entry, results, session_dir,
                   command=w.command).start()
    for t in range(8):
        run.tier(t)
    run.finish()
    w.session, w.results, w.run = session_dir, results, run
    return w


def test_without_a_reference_tiers_are_skipped_and_outputs_kept(first_run):
    st = {t: first_run.results.tier_status("ds9_blue_L1", t) for t in range(8)}
    assert st == {0: "pass", 1: "pass", 2: "skipped", 3: "skipped",
                  4: "skipped", 5: "pass", 6: "skipped", 7: "info"}
    assert (first_run.session / "ds9_blue_L1" / "redux").is_dir()


def bless_args(w, reason, **kw):
    return dict(data_dir=w.data, entry_id="ds9_blue_L1", reason=reason,
                run=w.session, assume_yes=True, **kw)


def test_bless_records_its_reason_and_creates_the_reference(first_run, capsys):
    w = first_run
    blessmod.bless(**bless_args(w, "initial: first reference"), today="2026-02-03")
    ref = w.data / "reference" / "ds9_blue_L1"
    assert sorted(p.name for p in ref.iterdir()) == [
        "kcwib.proc", "plots", "redux"]
    entry = mf.get_entry(mf.load_manifest(w.data), "ds9_blue_L1")
    hist = entry["bless_history"]
    assert len(hist) == 1
    assert hist[0]["reason"] == "initial: first reference"
    assert hist[0]["date"] == "2026-02-03"
    assert hist[0]["old_provenance"]["drp_version"] is None
    assert hist[0]["new_provenance"]["drp_git_version"]
    assert hist[0]["diff_summary"]["2"]["status"] == "skipped"
    assert entry["provenance"] == hist[0]["new_provenance"]
    assert "Bless ds9_blue_L1" in capsys.readouterr().out
    # a second bless replaces the reference and keeps both history entries
    (ref / "redux" / "stale.txt").write_text("old")
    blessmod.bless(**bless_args(w, "environment: new numpy"))
    assert not (ref / "redux" / "stale.txt").exists()
    assert not (w.data / "reference" / ".ds9_blue_L1.old").exists()
    assert len(mf.get_entry(mf.load_manifest(w.data),
                            "ds9_blue_L1")["bless_history"]) == 2


@pytest.mark.parametrize("reason", [
    None, "", "no colon here", "because: I want to", "environment",
    "environment:   ", "P.5: missing id", "M.9: no such stage"])
def test_bless_refuses_missing_or_disallowed_reasons(first_run, reason):
    with pytest.raises(blessmod.BlessError):
        blessmod.bless(**bless_args(first_run, reason))
    assert not (first_run.data / "reference" / "ds9_blue_L1").exists()


@pytest.mark.parametrize("reason", [
    "initial: x", "environment: x", "T.11 viewer-off: x", "bokeh port 4.3: x",
    "prompt text PR3a: x", "P.5 B1: x", "P.5 B33: x", "M.2c header: x",
    "M.2e initial: x", "M.3: x", "M.4: x", "M.5: x", "matrix default: x",
    "release version: x"])
def test_the_plans_bless_reasons_are_accepted(reason):
    # source: plan T.7 bless table (+ 'initial')
    blessmod.parse_reason(reason)


def test_bless_cli_without_reason_is_refused(first_run, capsys):
    rc = cli.main(["--data-dir", str(first_run.data), "bless", "--entry",
                   "ds9_blue_L1", "--run", str(first_run.session), "--yes"])
    assert rc == 2 and "--reason is required" in capsys.readouterr().err


def test_bless_after_outputs_were_deleted_is_refused(first_run):
    import shutil
    shutil.rmtree(first_run.session / "ds9_blue_L1" / "redux")
    with pytest.raises(blessmod.BlessError, match="outputs .* were deleted"):
        blessmod.bless(**bless_args(first_run, "initial: x"))


def test_bless_of_a_passing_torn_down_run_is_refused(world, tmp_path):
    session_dir = tmp_path / "s"
    results = SessionResults(session_dir)
    entry = mf.get_entry(world.mfst, "ds9_blue_L1")
    run = EntryRun(world.data, world.mfst, entry, results, session_dir,
                   command=world.command).start()
    for t in range(8):
        run.tier(t)
    run.finish()
    assert run.redux_deleted
    with pytest.raises(blessmod.BlessError, match="deleted"):
        blessmod.bless(world.data, "ds9_blue_L1", "environment: x",
                       run=session_dir, assume_yes=True)


def test_bless_refuses_a_run_whose_completion_failed(make_world, tmp_path,
                                                     command_factory):
    w = make_world()
    session_dir = tmp_path / "s"
    results = SessionResults(session_dir)
    entry = mf.get_entry(w.mfst, "ds9_blue_L1")
    run = EntryRun(w.data, w.mfst, entry, results, session_dir,
                   command=command_factory(w.script,
                                           fail_ids=("ds9_blue_L1",))).start()
    for t in range(8):
        run.tier(t)
    run.finish()
    with pytest.raises(blessmod.BlessError, match="tier 1 failed"):
        blessmod.bless(w.data, "ds9_blue_L1", "initial: x", run=session_dir,
                       assume_yes=True)


def test_declined_confirmation_changes_nothing(first_run):
    with pytest.raises(blessmod.BlessError, match="not confirmed"):
        blessmod.bless(first_run.data, "ds9_blue_L1", "initial: x",
                       run=first_run.session, input_fn=lambda prompt: "no")
    assert not (first_run.data / "reference" / "ds9_blue_L1").exists()
    assert mf.get_entry(mf.load_manifest(first_run.data),
                        "ds9_blue_L1")["bless_history"] == []


def test_bless_replaces_a_symlinked_reference_without_touching_its_target(
        first_run, tmp_path):
    w = first_run
    target = tmp_path / "s3_run"
    (target / "redux").mkdir(parents=True)
    (target / "redux" / "keep.txt").write_text("S.3 output")
    (w.data / "reference").mkdir()
    (w.data / "reference" / "ds9_blue_L1").symlink_to(target)
    blessmod.bless(**bless_args(w, "environment: x"))
    assert (target / "redux" / "keep.txt").read_text() == "S.3 output"
    ref = w.data / "reference" / "ds9_blue_L1"
    assert not ref.is_symlink() and (ref / "redux" / "ds_icubes.fits").is_file()


def test_blessing_an_against_entry_turns_it_into_self(first_run):
    w = first_run
    m = mf.load_manifest(w.data)
    e = mf.get_entry(m, "ds9_blue_L1")
    e["compare"] = {"mode": "against", "entry": "ds9_blue_L0",
                    "informational_tiers": [3, 4]}
    mf.save_manifest(w.data, m)
    blessmod.bless(**bless_args(w, "M.2e initial: first matrix cube"))
    assert mf.get_entry(mf.load_manifest(w.data),
                        "ds9_blue_L1")["compare"] == {"mode": "self"}
