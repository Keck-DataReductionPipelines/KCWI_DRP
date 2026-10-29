"""Tier 1 (completion): channel-specific allowed ERROR lines and object cubes.

The runs come from the fake reduction in conftest.py (FAKE_SCRIPT), which
writes the two benign ERROR lines only for red (or with --benign-errors) and
leaves out the OBJECT cube with --drop-cube.
"""
import subprocess
import sys

import pytest
from astropy.table import Table

from kcwidrp.tests.regression import compare
from kcwidrp.tests.regression import manifest as mf

CORRECT = ("2026-01-01 10:00:02:KCWI:ERROR:CorrectDefects.py: Defect list "
           "not found for /x/data/defect_L2U2L1U1_2x2.dat")
MAKECUBE = ("2026-01-01 10:00:03:KCWI:ERROR:MakeCube.py: No reference "
            "geometry, cannot make cube!")


def fake_run(world, tmp_path, *args):
    """Run the fake reduction in tmp_path/run; return that directory."""
    run = tmp_path / "run"
    run.mkdir()
    subprocess.run([sys.executable, str(world.script)] + list(args), cwd=str(run),
                   check=True, stdin=subprocess.DEVNULL)
    return run


def entry_for(world, channel):
    entry = dict(mf.get_entry(world.mfst, "ds9_blue_L1"))
    entry["channel"] = channel
    return entry


def hand_tree(tmp_path, proc_name, types, stems, cubes):
    """A run dir with a fixed-width proctab and the listed cubes."""
    run = tmp_path / "hand"
    (run / "redux").mkdir(parents=True)
    Table({"TYPE": types, "filename": ["/d/lev0/%s.fits" % s for s in stems]}
          ).write(str(run / proc_name), format="ascii.fixed_width")
    for stem in cubes:
        (run / "redux" / (stem + "_icube.fits")).write_bytes(b"x")
    return run


def test_red_run_with_benign_errors_passes_tier1(world, tmp_path):
    run = fake_run(world, tmp_path, "-r")
    r = compare.tier1_result(run, entry_for(world, "red"))
    assert r.status == "pass", r.details


def test_blue_run_with_benign_errors_fails_tier1(world, tmp_path):
    run = fake_run(world, tmp_path, "--benign-errors")
    r = compare.tier1_result(run, entry_for(world, "blue"))
    assert r.status == "fail"
    text = "\n".join(r.details)
    # source: the two benign lines written by FAKE_SCRIPT, not allowed on blue
    assert "unexpected ERROR line" in text
    assert "MakeCube" in text and "CorrectDefects" in text


def test_clean_blue_run_passes_tier1(world, tmp_path):
    run = fake_run(world, tmp_path)
    r = compare.tier1_result(run, entry_for(world, "blue"))
    assert r.status == "pass", r.details


def test_scan_logs_channel_none_blue_red(tmp_path):
    logs = tmp_path / "logs"
    logs.mkdir()
    (logs / "KCWI.log").write_text(
        "t:KCWI:INFO: ok\n%s\n%s\n"
        "t:KCWI:ERROR:Failed executing primitive Foo: boom\n" % (
            CORRECT, MAKECUBE))
    (logs / "DRPF.log").write_text("t:DRPF:ERROR: x while invoking Bar\n")
    # source: 2 failure lines (pattern match), 2 benign ERROR lines, 2 logs
    for ch, n_err in ((None, 2), ("blue", 2), ("red", 0)):
        failures, errors, nlogs = compare.scan_logs(tmp_path, ch)
        assert len(failures) == 2 and nlogs == 2
        assert len(errors) == n_err
    assert compare.scan_logs(tmp_path)[1] == compare.scan_logs(tmp_path, None)[1]
    _f, errors, _n = compare.scan_logs(tmp_path, "blue")
    assert any("CorrectDefects" in e for e in errors)
    assert any("MakeCube" in e for e in errors)


def test_unrelated_error_line_fails_red_run(world, tmp_path):
    run = fake_run(world, tmp_path, "-r")
    with open(run / "logs" / "KCWI.log", "a") as f:
        f.write("t:KCWI:ERROR:Something new broke\n")
    assert compare.scan_logs(run, "red")[1] == [
        "KCWI.log: t:KCWI:ERROR:Something new broke"]
    r = compare.tier1_result(run, entry_for(world, "red"))
    assert r.status == "fail" and "Something new broke" in "\n".join(r.details)


def test_dropped_cube_fails_tier1_without_reference(world, tmp_path):
    run = fake_run(world, tmp_path, "--drop-cube")
    r = compare.tier1_result(run, entry_for(world, "blue"))
    assert r.status == "fail"
    assert "KB.3" in "\n".join(r.details)


def test_dropped_cube_fails_tier1_with_reference(world, tmp_path):
    ref = world.data / "reference" / "ds9_blue_L1"
    assert ref.is_dir()  # the world fixture blessed a reference
    run = fake_run(world, tmp_path, "--drop-cube")
    r = compare.tier1_result(run, mf.get_entry(world.mfst, "ds9_blue_L1"))
    assert r.status == "fail" and "KB.3" in "\n".join(r.details)


def test_missing_cubes_all_present(tmp_path):
    run = hand_tree(tmp_path, "kcwib.proc", ["BIAS", "MARC", "OBJECT"],
                    ["KB.1", "KB.2", "KB.3"], ["KB.2", "KB.3"])
    assert compare.missing_cubes(run, "blue") == []


def test_missing_cubes_names_missing_marc(tmp_path):
    run = hand_tree(tmp_path, "kcwir.proc", ["BIAS", "MARC", "OBJECT"],
                    ["KR.1", "KR.2", "KR.3"], ["KR.3"])
    problems = compare.missing_cubes(run, "red")
    assert len(problems) == 1 and "KR.2" in problems[0]


def test_missing_cubes_uses_the_channels_own_proctab(tmp_path):
    # a blue proctab does not satisfy a red check
    run = hand_tree(tmp_path, "kcwib.proc", ["OBJECT"], ["KB.3"], ["KB.3"])
    assert len(compare.missing_cubes(run, "red")) == 1


def test_missing_cubes_no_proctab_is_one_problem(tmp_path):
    (tmp_path / "redux").mkdir()
    assert len(compare.missing_cubes(tmp_path, "blue")) == 1


def test_missing_cubes_proctab_without_object_rows_is_one_problem(tmp_path):
    run = hand_tree(tmp_path, "kcwib.proc", ["BIAS", "BIAS"],
                    ["KB.1", "KB.2"], [])
    assert len(compare.missing_cubes(run, "blue")) == 1


def test_missing_cubes_records_cubes_expected(tmp_path):
    run = hand_tree(tmp_path, "kcwib.proc",
                    ["BIAS", "ARCLAMP", "MARC", "OBJECT", "OBJECT"],
                    ["KB.1", "KB.2", "KB.2", "KB.3", "KB.4"], ["KB.2"])
    data = {}
    problems = compare.missing_cubes(run, "blue", data)
    # source: one MARC + two OBJECT rows expected; KB.3 and KB.4 cubes absent
    assert data["cubes_expected"] == 3 and len(problems) == 2
