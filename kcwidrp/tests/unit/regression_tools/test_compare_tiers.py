"""Comparison tiers 1-7 on synthetic run/reference trees (T.7 validation).

The reference and the run come from the fake reduction in conftest.py, so the
unmodified run equals its reference except for the header content that tier
2/4 must ignore (HISTORY, a path string).
"""
import shutil

import numpy as np
import pytest
from astropy.io import fits

from kcwidrp.tests.regression import compare
from kcwidrp.tests.regression import manifest as mf

PROC = "kcwib.proc"


@pytest.fixture
def run_copy(world, baseline_run, tmp_path):
    dst = tmp_path / "run"
    shutil.copytree(baseline_run, dst)
    return dst


@pytest.fixture
def ref_dir(world):
    return world.data / "reference" / "ds9_blue_L1"


def tiers(run, ref, tol=None, diag=None):
    c = compare.compare_outputs(run, ref, PROC, tol, diag)
    return {2: compare.structure_result(c), 3: compare.pixels_result(c),
            4: compare.headers_result(c), 6: compare.plots_result(c)}


def edit_fits(path, func):
    with fits.open(path, mode="update") as hdul:
        func(hdul)


def test_run_equals_its_reference_in_every_tier(run_copy, ref_dir):
    res = tiers(run_copy, ref_dir)
    # source: analytic -- the fake reduction is deterministic (seeded RNG), so
    # the only differences are the ignored time/path cards.
    assert {t: r.status for t, r in res.items()} == {
        2: "pass", 3: "pass", 4: "pass", 6: "pass"}


def test_scaling_one_extension_by_one_ppm_fails_tier3_naming_it(
        run_copy, ref_dir, tmp_path):
    # source: analytic -- (1 + 1e-6) in float32 changes nearly every value.
    edit_fits(run_copy / "redux" / "ds_icubes.fits",
              lambda h: h["FLUX"].data.__imul__(np.float32(1 + 1e-6)))
    res = tiers(run_copy, ref_dir, diag=tmp_path / "diff")
    assert res[3].status == "fail"
    assert "ds_icubes.fits[FLUX]" in res[3].details[0]
    assert res[2].status == res[4].status == "pass"
    f = res[3].data["findings"][0]
    assert (f["file"], f["ext"]) == ("ds_icubes.fits", "FLUX")
    # the diagnostic image of the worst plane was written (Agg, no display)
    assert compare.png_size(f["diff_png"]) is not None   # a valid PNG
    assert str(tmp_path / "diff") in f["diff_png"]


def test_manifest_tolerance_lets_a_ppm_change_pass(run_copy, ref_dir):
    edit_fits(run_copy / "redux" / "ds_icubes.fits",
              lambda h: h["FLUX"].data.__imul__(np.float32(1 + 1e-6)))
    tol = {"ds_icubes.fits::FLUX": {"rtol": 1e-5, "atol": 0.0,
                                    "reason": "test"}}
    assert tiers(run_copy, ref_dir, tol)[3].status == "pass"
    # a tolerance for another extension does not apply
    other = {"ds_icubes.fits::MASK": tol["ds_icubes.fits::FLUX"]}
    assert tiers(run_copy, ref_dir, other)[3].status == "fail"


def test_integer_extension_is_exact(run_copy, ref_dir):
    def bump(h):
        h["MASK"].data[0, 0, 0] = 1
    edit_fits(run_copy / "redux" / "ds_icubes.fits", bump)
    res = tiers(run_copy, ref_dir,
                {"*": {"rtol": 0.5, "atol": 0.5, "reason": "huge"}})
    assert res[3].status == "fail"
    assert "ds_icubes.fits[MASK]" in res[3].details[0]


def test_deleted_output_fails_tier2(run_copy, ref_dir):
    (run_copy / "redux" / "ds_int.fits").unlink()
    res = tiers(run_copy, ref_dir)
    assert res[2].status == "fail"
    assert "redux/ds_int.fits" in res[2].details[0]


def test_extra_output_and_changed_shape_fail_tier2(run_copy, ref_dir):
    shutil.copy(run_copy / "redux" / "ds_int.fits",
                run_copy / "redux" / "extra.fits")
    res = tiers(run_copy, ref_dir)
    assert res[2].status == "fail" and "extra.fits" in res[2].details[0]
    (run_copy / "redux" / "extra.fits").unlink()
    fits.PrimaryHDU(np.zeros((5, 5))).writeto(
        run_copy / "redux" / "ds_int.fits", overwrite=True)
    res = tiers(run_copy, ref_dir)
    assert res[2].status == "fail"
    assert any("shape/dtype" in d for d in res[2].details)


def test_changed_header_value_fails_tier4_only(run_copy, ref_dir):
    def change(h):
        h[0].header["CWAVE"] = 4500.6
    edit_fits(run_copy / "redux" / "ds_icubes.fits", change)
    res = tiers(run_copy, ref_dir)
    assert res[4].status == "fail"
    assert "CWAVE" in res[4].details[0]
    assert res[2].status == res[3].status == "pass"


def test_changed_string_header_value_fails_tier4(run_copy, ref_dir):
    def change(h):
        h[0].header["OBJECT"] = "galaxy"
    edit_fits(run_copy / "redux" / "ds_icubes.fits", change)
    assert tiers(run_copy, ref_dir)[4].status == "fail"


def test_added_header_keyword_fails_tier2(run_copy, ref_dir):
    def add(h):
        h[0].header["NEWKEY"] = 1
    edit_fits(run_copy / "redux" / "ds_icubes.fits", add)
    res = tiers(run_copy, ref_dir)
    assert res[2].status == "fail" and "NEWKEY" in res[2].details[0]


def test_ignored_cards_do_not_matter(run_copy, ref_dir):
    def change(h):
        h[0].header["INFILE"] = "/somewhere/else/files.txt"
        h[0].header.add_history("another line")
    edit_fits(run_copy / "redux" / "ds_icubes.fits", change)
    res = tiers(run_copy, ref_dir)
    assert res[2].status == res[4].status == "pass"


def test_changed_observation_date_fails_tier4(run_copy, ref_dir):
    # source: analytic -- DATE* cards are copied from the raw header, so a
    # change means the output came from different input; it must not be
    # ignored (compare.py docstring).
    def change(h):
        h[0].header["DATE-OBS"] = "1999-01-01T00:00:00"
    edit_fits(run_copy / "redux" / "ds_icubes.fits", change)
    res = tiers(run_copy, ref_dir)
    assert res[4].status == "fail" and "DATE-OBS" in " ".join(res[4].details)


def test_proctab_column_difference_fails_tier3_but_directory_does_not(
        run_copy, ref_dir):
    p = run_copy / PROC
    text = p.read_text()
    p.write_text(text.replace(str(run_copy), "/elsewhere"))
    assert tiers(run_copy, ref_dir)[3].status == "pass"   # base names equal
    p.write_text(text.replace("ARCLAMP", "DOMEFLAT"))
    res = tiers(run_copy, ref_dir)
    assert res[3].status == "fail" and "TYPE" in res[3].details[0]


def test_non_fits_file_difference_fails_tier3(run_copy, ref_dir):
    (run_copy / "redux" / "ds_atlines.txt").write_text("other\n")
    assert tiers(run_copy, ref_dir)[3].status == "fail"


def test_plot_name_and_size_fail_tier6(run_copy, ref_dir, world):
    (run_copy / "plots" / "p1.png").unlink()
    res = tiers(run_copy, ref_dir)
    assert res[6].status == "fail" and "p1.png" in res[6].details[0]
    # a different size: copy p2 (9x7) over p1 (8x6)
    shutil.copy(run_copy / "plots" / "p2.png", run_copy / "plots" / "p1.png")
    res = tiers(run_copy, ref_dir)
    assert res[6].status == "fail"
    assert "(9, 7)" in res[6].details[0] and "(8, 6)" in res[6].details[0]


def test_diff_arrays_statistics():
    a = np.zeros((2, 3, 4))
    b = a.copy()
    b[1, 2, 3] = 0.5
    b[0, 0, 0] = -0.25
    s = compare.diff_arrays(a, b)
    # source: analytic -- two of 24 values differ; the larger |diff| is 0.5 at
    # index (1, 2, 3); relative difference is |a-b|/|b| = 1 at both.
    assert (s["n_diff"], s["location"]) == (2, (1, 2, 3))
    assert s["max_abs"] == 0.5 and s["fraction"] == 2 / 24
    assert s["max_rel"] == 1.0
    assert compare.diff_arrays(a, a) is None
    nan = a.copy()
    nan[0, 0, 0] = np.nan
    assert compare.diff_arrays(nan, nan) is None       # NaN equals NaN
    assert compare.diff_arrays(nan, a)["n_nan_mismatch"] == 1


# -- tier 1 ------------------------------------------------------------------

def entry_l1(world):
    return mf.get_entry(world.mfst, "ds9_blue_L1")


def test_tier1_passes_on_clean_run_with_benign_errors(world, run_copy):
    # the fake log carries the two documented benign ERROR lines
    r = compare.tier1_result(run_copy, entry_l1(world))
    assert r.status == "pass", r.details


def test_tier1_fails_on_primitive_failure_naming_it_though_exit_is_zero(
        world, run_copy):
    log = run_copy / "logs" / "KCWI.log"
    log.write_text(log.read_text() +
                   "t:KCWI:ERROR:Failed executing primitive MakeCube: boom\n")
    (run_copy / "logs" / "DRPF.log").open("a").write(
        "t:DRPF:ERROR: Exception boom while invoking WavelengthCorrections\n")

    class Res:
        status, exit_code, reason = "ok", 0, ""
    r = compare.tier1_result(run_copy, entry_l1(world), run_result=Res())
    assert r.status == "fail"
    assert r.data["failed_primitives"] == ["MakeCube", "WavelengthCorrections"]
    assert r.data["exit_code"] == 0


def test_tier1_fails_on_unlisted_error_line(world, run_copy):
    log = run_copy / "logs" / "KCWI.log"
    log.write_text(log.read_text() + "t:KCWI:ERROR:Something new broke\n")
    r = compare.tier1_result(run_copy, entry_l1(world))
    assert r.status == "fail" and "Something new broke" in r.details[0]


def test_tier1_fails_without_logs_or_with_prompts_at_level_1(world, run_copy):
    (run_copy / "prompt_record.jsonl").write_text(
        '{"prompt": "go? ", "answer": ""}\n')
    r = compare.tier1_result(run_copy, entry_l1(world))
    assert r.status == "fail" and "unexpected prompt" in r.details[0]
    shutil.rmtree(run_copy / "logs")
    (run_copy / "prompt_record.jsonl").unlink()
    r = compare.tier1_result(run_copy, entry_l1(world))
    assert r.status == "fail" and "no logs" in r.details[0]


# -- tiers 5 and 7 -----------------------------------------------------------

def test_tier5_notes_no_metrics_and_runs_registered_ones(world, run_copy,
                                                        monkeypatch):
    r = compare.tier5_result(run_copy, None, entry_l1(world))
    assert r.status == "pass" and "no metrics registered" in r.summary
    monkeypatch.setitem(compare.METRICS, "m", lambda *a: (False, "too big"))
    r = compare.tier5_result(run_copy, None, entry_l1(world))
    assert r.status == "fail" and r.details == ["m: too big"]


def test_tier7_timings_from_event_gaps(run_copy):
    r = compare.tier7_result(run_copy, 3.5)
    # source: analytic -- fake DRPF.log events: start 10:00:00, make_cube at
    # :05 and :12, other at :14. Gaps: make_cube 5 + 7 = 12 s in 2 calls,
    # other 2 s, start 0 s.
    prim = r.data["primitives"]
    assert prim["make_cube"] == {"seconds": 12.0, "calls": 2}
    assert prim["other"]["seconds"] == 2.0 and prim["start"]["seconds"] == 0.0
    assert r.status == "info" and r.data["wall_time"] == 3.5


def test_informational_tiers_are_reported_not_failed():
    entry = {"compare": {"mode": "against", "entry": "x",
                         "informational_tiers": [3, 4]}}
    r = compare.apply_informational(compare.TierResult(3, "fail", "bad"), entry)
    assert r.status == "info"
    r = compare.apply_informational(compare.TierResult(2, "fail", "bad"), entry)
    assert r.status == "fail"
