"""Value prompts: what each accepts today, and the value the reduction uses.

Plan 8 (D3, D12): for input that is valid today, Track P's conversion of the
value prompts (P.4) must give identical results.  These tests run each
primitive at the plot level that reaches the prompt, type the answers of a
table of VALID inputs (docs/interactive_prompts.md, "Accepted input (D12)"
of each prompt) and check the value the reduction ends up using: the
primitive's action argument, an output header keyword, or the output data.
They pass on the code before P.4 and must keep passing after it.

Only the *number* of prompts is compared (D3), never their text: P.3/P.4
reword the prompts.  Input that is invalid today (B1c, B2, B4, B5, B6) and
the prompts whose valid behaviour is a bug fixed in PR3b (B1, B3, B7, B12)
are not exercised here, except where noted.

Sites covered: ReadAtlas offset (P40), MakeInvsens limits (P14), line mask
(P15), fit order (P19), MakeMasterFlat peak (P25), NandshuffSubtractSky
scale (P37).  The red non-standard branch (P38) is reached, but its value cannot be
read back (see ``test_red_branch_asks_for_the_scale_once``).
"""
from types import SimpleNamespace

import numpy as np
import pytest
from astropy.io import fits

from kcwidrp.core.bspline import Bspline as real_bspline
from kcwidrp.core.kcwi_pkg_resources import get_resource_path
from kcwidrp.primitives import MakeInvsens as make_invsens_module
from kcwidrp.primitives import MakeMasterFlat as make_master_flat_module
from kcwidrp.primitives import ReadAtlas as read_atlas_module
from kcwidrp.primitives.MakeInvsens import MakeInvsens
from kcwidrp.primitives.MakeMasterFlat import MakeMasterFlat, \
    bm_ledge_position
from kcwidrp.primitives.NandshuffSubtractSky import NandshuffSubtractSky
from kcwidrp.primitives.ReadAtlas import ReadAtlas
from kcwidrp.primitives.kcwi_file_primitives import KCCDData, \
    kcwi_fits_writer

HEAD = "e6334c8"


@pytest.fixture(autouse=True)
def no_png_export(monkeypatch):
    # save_plot exports PNGs through selenium and a browser: not wanted here
    for module in (read_atlas_module, make_invsens_module,
                   make_master_flat_module):
        monkeypatch.setattr(module, "save_plot", lambda *a, **k: None)


# ======================================================== ReadAtlas (P40)
# Input: the arc of the reference bar is the FeAr atlas itself, sampled on
# the preliminary wavelength scale (0.2 A / px, as the atlas), so the
# cross-correlation has a single sharp peak.  The calculated offset is not
# asserted; only that <cr> keeps it and that typed integers replace it.
AT_NX, AT_DISP, AT_CWAVE = 2000, 0.2, 5000.0
NEXT_ATLAS = 2          # prompts before the loop: "Next? <cr>" (1), then loop


def atlas_refdisp():
    # source: independent - the atlas file's own header (astropy.io.fits)
    with fits.open(get_resource_path("kcwidrp", "data/fear.fits")) as hdul:
        return hdul[0].header["CDELT1"]


def run_atlas(run_primitive, synthetic_ccddata, answers, plot_level=2,
              **overrides):
    with fits.open(get_resource_path("kcwidrp", "data/fear.fits")) as hdul:
        flux = hdul[0].data.astype(float)
        head = hdul[0].header
        wave = head["CRVAL1"] + head["CDELT1"] * np.arange(len(flux))
    obswav = (np.arange(AT_NX) - AT_NX // 2) * AT_DISP + AT_CWAVE
    arc = np.interp(obswav, wave, flux)
    frame = synthetic_ccddata(IMTYPE="MARC", BGRATNAM="BM")
    config = {"plot_level": plot_level, "REFBAR": 57}
    config.update(overrides)
    return run_primitive(
        ReadAtlas, frame,
        context_extra={"arcs": [arc.copy() for _ in range(60)],
                       "prelim_disp": AT_DISP},
        action_extra={"illum": "FEAR", "atsig": 2.0, "cwave": AT_CWAVE},
        config_overrides=config, answers=answers)


@pytest.mark.parametrize("typed, expected", [
    ("12", 12), ("-3", -3), ("+4", 4), (" 5 ", 5), ("0", 0),
    ("007", 7), ("100000", 100000),
], ids=["12", "-3", "+4", "padded", "0", "leading-zeros", "large"])
def test_atlas_offset_valid_integers_are_used(run_primitive,
                                              synthetic_ccddata, typed,
                                              expected):
    out = run_atlas(run_primitive, synthetic_ccddata, ["", typed, ""])
    assert len(run_primitive.last.prompts) == 3
    # source: analytic - P40: the typed integer is the offset in atlas px
    assert out.offset_pix == expected
    # source: independent - offset_wave = offset_pix * CDELT1 of the atlas
    # file.  Tolerance: one multiplication of float64 values, exact.
    assert out.offset_wave == expected * atlas_refdisp()


def test_atlas_offset_last_answer_wins(run_primitive, synthetic_ccddata):
    out = run_atlas(run_primitive, synthetic_ccddata,
                    ["", "12", "-3", "8", ""])
    assert len(run_primitive.last.prompts) == 5
    assert out.offset_pix == 8      # source: analytic - the last number typed


@pytest.mark.parametrize("bad", ["1.5", "abc", "q", "1e2", "4 5"])
def test_atlas_offset_invalid_text_is_asked_again_and_value_kept(
        run_primitive, synthetic_ccddata, bad):
    # source: analytic - P40: text and decimals are refused with "Try again"
    # and the offset accepted before stays in force
    out = run_atlas(run_primitive, synthetic_ccddata,
                    ["", "12", bad, "", ])
    assert len(run_primitive.last.prompts) == 4
    assert out.offset_pix == 12
    assert out.offset_wave == 12 * atlas_refdisp()


def test_atlas_offset_accepting_keeps_calculated_offset(run_primitive,
                                                        synthetic_ccddata):
    # <cr> at the entry keeps the calculated offset: the same one the
    # reduction uses at plot_level 1, where nothing is asked
    out2 = run_atlas(run_primitive, synthetic_ccddata, ["", ""])
    assert len(run_primitive.last.prompts) == 2
    out1 = run_atlas(run_primitive, synthetic_ccddata, None, plot_level=1)
    assert run_primitive.last.prompts == []
    assert out2.offset_pix == out1.offset_pix
    assert out2.offset_wave == out1.offset_wave
    # source: analytic - the arc is the atlas on a scale that is exact for
    # it, so the calculated offset is small compared with the 12 px typed
    # in the other tests, and it is an integer number of atlas pixels
    assert abs(out2.offset_pix) < 12


def test_atlas_offset_requested_with_atoff_is_the_default(run_primitive,
                                                          synthetic_ccddata):
    # source: analytic - config ATOFF is the requested offset in px
    out = run_atlas(run_primitive, synthetic_ccddata, ["", ""], ATOFF=6)
    assert (out.offset_pix, out.offset_wave) == (6, 6 * atlas_refdisp())


# ===================================================== MakeInvsens (P14,
# P15, P19).  Same synthetic standard-star cube as test_run_MakeInvsens.py
# (copied, not imported): a star in slice 2 at y = 10 in 600 wavelength
# planes, 4000 - 4599 A, WAVGOOD 4050 - 4550.
STD = "bd28d4211"
NW, NY, NSL = 600, 20, 5
WGOOD = (4050.0, 4550.0)
WALL = (4000.0, 4599.0)
INVSENS = "kb260101_00001_invsens.fits"
NEXT = ""           # <cr> for the "Next?" pauses
FIT_DONE = [NEXT, NEXT, NEXT]       # P16 - P18 before each P19


def std_cube(synthetic_ccddata):
    frame = synthetic_ccddata(
        IMTYPE="OBJECT", OBJECT=STD, TARGNAME=STD, DARCOR=True,
        XPOSURE=100.0, CRVAL3=WALL[0], CD3_3=1.0, CRPIX3=1.0,
        WAVGOOD0=WGOOD[0], WAVGOOD1=WGOOD[1], WAVALL0=WALL[0],
        WAVALL1=WALL[1], DARPADY=2, SKYCOR=True, AIRMASS=1.0, DICHFRAC=1.0,
        BGRATNAM="BL")
    w = WALL[0] + np.arange(NW)
    spec = 1.0e4 * (1.0 + 0.3 * np.sin((w - 4000.0) / 150.0))
    prof = np.exp(-0.5 * ((np.arange(NY) - 10.0) / 1.5) ** 2)
    cube = np.zeros((NW, NY, NSL))
    for sl, amp in ((1, 0.3), (2, 1.0), (3, 0.3)):
        cube[:, :, sl] = amp * spec[:, None] * prof[None, :]
    cube += 1.0
    return KCCDData(cube, meta=frame.header, unit="electron",
                    uncertainty=np.sqrt(cube))


def run_invsens(run_primitive, synthetic_ccddata, tmp_redux, answers):
    run_primitive(MakeInvsens, std_cube(synthetic_ccddata),
                  config_overrides={"plot_level": 2}, answers=answers)
    with fits.open(tmp_redux.redux / INVSENS) as hdul:
        return hdul[0].header.copy()


pytestmark_invsens = pytest.mark.filterwarnings(
    "ignore:Polyfit may be poorly conditioned")


@pytestmark_invsens
@pytest.mark.parametrize("typed, lo, hi", [
    ("4100 4500", 4100.0, 4500.0),
    ("4100.5 4500.25", 4100.5, 4500.25),
    ("4.1e3 4.5e3", 4100.0, 4500.0),
    ("  4100   4500  ", 4100.0, 4500.0),
    ("4100 4500 trailing words are ignored", 4100.0, 4500.0),
    ("4000 4599", 4000.0, 4599.0),
], ids=["ints", "decimals", "exponent", "spaces", "extra-words",
        "whole-range"])
def test_fit_limits_valid_pair_sets_invfw(run_primitive, synthetic_ccddata,
                                          tmp_redux, typed, lo, hi):
    hdr = run_invsens(run_primitive, synthetic_ccddata, tmp_redux,
                      [typed, "", "", ] + FIT_DONE + [""])
    # D3: limits (2 asks), masks (1), P16-P18, order (1)
    assert len(run_primitive.last.prompts) == 2 + 1 + 3 + 1
    # source: analytic - P14: the typed pair is the fit range, recorded in
    # INVFW0 / INVFW1 (WAVALL0 <= lo <= hi <= WAVALL1 here)
    assert (hdr["INVFW0"], hdr["INVFW1"]) == (lo, hi)


@pytestmark_invsens
def test_fit_limits_cr_keeps_wavgood_and_last_pair_wins(run_primitive,
                                                         synthetic_ccddata,
                                                         tmp_redux):
    hdr = run_invsens(run_primitive, synthetic_ccddata, tmp_redux,
                      ["", ""] + FIT_DONE + [""])
    assert (hdr["INVFW0"], hdr["INVFW1"]) == WGOOD     # source: P14 default
    hdr = run_invsens(run_primitive, synthetic_ccddata, tmp_redux,
                      ["4100 4500", "4200 4400", "", ""] + FIT_DONE + [""])
    assert len(run_primitive.last.prompts) == 3 + 1 + 3 + 1
    assert (hdr["INVFW0"], hdr["INVFW1"]) == (4200.0, 4400.0)


@pytestmark_invsens
@pytest.mark.parametrize("bad", ["4500 4100", "3900 4500", "4100 4700"],
                         ids=["reversed", "below-range", "above-range"])
def test_fit_limits_refused_pair_leaves_wavgood(run_primitive,
                                                synthetic_ccddata,
                                                tmp_redux, bad):
    # source: analytic - P14: a pair that is reversed or outside
    # WAVALL0 - WAVALL1 is refused and asked again; no limits were accepted
    # before, so the limits are still WAVGOOD (the B6 case, with an accepted
    # pair before the bad one, is not asserted here)
    hdr = run_invsens(run_primitive, synthetic_ccddata, tmp_redux,
                      [bad, "", ""] + FIT_DONE + [""])
    assert len(run_primitive.last.prompts) == 2 + 1 + 3 + 1
    assert (hdr["INVFW0"], hdr["INVFW1"]) == WGOOD


@pytestmark_invsens
def test_line_masks_are_saved_and_asked_until_cr(run_primitive,
                                                  synthetic_ccddata,
                                                  tmp_redux):
    answers = ["",                                     # limits: accept
               "4200 4210 Hbeta wing", "4300 4320",    # two masks
               "4100 4300",                            # inside the limits
               "",                                     # done
               ] + FIT_DONE + [""]
    run_invsens(run_primitive, synthetic_ccddata, tmp_redux, answers)
    assert len(run_primitive.last.prompts) == 1 + 4 + 3 + 1
    # source: characterization @ e6334c8 (P15) - the saved file has one line
    # per mask, "%.2f %.2f" and "# comment" when a comment was typed; the
    # numbers are the typed ones
    text = (tmp_redux.root / (STD + ".lmsk")).read_text()
    assert text == ("4200.00 4210.00 # Hbeta wing\n"
                    "4300.00 4320.00\n"
                    "4100.00 4300.00\n")


@pytestmark_invsens
@pytest.mark.parametrize("typed", ["3000 3100", "4010 4020", "4540 4560"],
                         ids=["outside", "below-limits", "partly-outside"])
def test_line_mask_outside_the_limits_is_not_added(run_primitive,
                                                    synthetic_ccddata,
                                                    tmp_redux, typed):
    # source: analytic - P15: both ends must lie strictly inside the fit
    # limits (here WAVGOOD 4050 - 4550); otherwise nothing is added, and
    # with no mask no .lmsk file is written
    run_invsens(run_primitive, synthetic_ccddata, tmp_redux,
                ["", typed, ""] + FIT_DONE + [""])
    assert len(run_primitive.last.prompts) == 1 + 2 + 3 + 1
    assert not (tmp_redux.root / (STD + ".lmsk")).exists()


@pytestmark_invsens
def test_line_mask_must_lie_inside_new_limits(run_primitive,
                                               synthetic_ccddata, tmp_redux):
    # source: analytic - P15: "inside the current fit limits": 4120 - 4130
    # is inside WAVGOOD but outside the limits typed first
    run_invsens(run_primitive, synthetic_ccddata, tmp_redux,
                ["4200 4500", "", "4120 4130", "4300 4310", ""]
                + FIT_DONE + [""])
    assert (tmp_redux.root / (STD + ".lmsk")).read_text() == \
        "4300.00 4310.00\n"


@pytestmark_invsens
@pytest.mark.parametrize("typed, order", [
    ("5", 5), ("+5", 5), (" 5 ", 5), ("3", 3), ("12", 12), ("05", 5),
], ids=["5", "+5", "padded", "3", "12", "leading-zero"])
def test_fit_order_valid_integer_sets_invford(run_primitive,
                                              synthetic_ccddata, tmp_redux,
                                              typed, order):
    hdr = run_invsens(run_primitive, synthetic_ccddata, tmp_redux,
                      ["", ""] + FIT_DONE + [typed] + FIT_DONE + [""])
    # D3: limits, masks, P16-18, order, P16-18, order
    assert len(run_primitive.last.prompts) == 2 + 3 + 1 + 3 + 1
    assert hdr["INVFORD"] == order       # source: analytic - P19


@pytestmark_invsens
def test_fit_order_default_and_last_answer_wins(run_primitive,
                                                 synthetic_ccddata,
                                                 tmp_redux):
    hdr = run_invsens(run_primitive, synthetic_ccddata, tmp_redux,
                      ["", ""] + FIT_DONE + [""])
    # source: P19 help text - default 7 for the BL grating times DICHFRAC (1)
    assert hdr["INVFORD"] == 7
    hdr = run_invsens(run_primitive, synthetic_ccddata, tmp_redux,
                      ["", ""] + FIT_DONE + ["5"] + FIT_DONE + ["4"]
                      + FIT_DONE + [""])
    assert len(run_primitive.last.prompts) == 2 + 3 + 1 + 3 + 1 + 3 + 1
    assert hdr["INVFORD"] == 4


@pytestmark_invsens
@pytest.mark.parametrize("bad", ["5.5", "abc", "q", "4 5"])
def test_fit_order_invalid_text_is_asked_again(run_primitive,
                                                synthetic_ccddata, tmp_redux,
                                                bad):
    # source: analytic - P19: text and decimals are refused and the order is
    # asked for again; the next valid order is used.  Before Track P a refused
    # order also showed P16-P18 again; P.4 asks again at once, which D12
    # allows for input that is invalid today (one P19 prompt fewer per pass)
    hdr = run_invsens(run_primitive, synthetic_ccddata, tmp_redux,
                      ["", ""] + FIT_DONE + [bad, "6"] + FIT_DONE + [""])
    assert len(run_primitive.last.prompts) == 2 + 4 + 1 + 4
    assert hdr["INVFORD"] == 6


# ================================================== MakeMasterFlat (P25)
# Same synthetic BLUE BM dome flat as test_run_MakeMasterFlat.py (copied):
# a +20 % ledge step at LEDGE + 5 A and a lower +3 % one at LEDGE - 8 A, so
# the derivative has two peaks and the prompt is reached at plot_level 2.
NROW, NCOL, DROW = 800, 40, 0.1
BCWAVE = 4500.0
LEDGE = bm_ledge_position(BCWAVE, True)
EXTRA = {"new_type": "MDOME", "stack_type": "SDOME", "min_files": 1}
BIG = (LEDGE + 5.0, 2000.0)
SMALL = (LEDGE - 8.0, 300.0)


class ReachedMasterFit(Exception):
    pass


@pytest.fixture
def capture(monkeypatch):
    """Capture the data at the master-flat fit and stop the run there."""
    got = SimpleNamespace(x=None, y=None)

    def iterfit(xdata, ydata, *args, **kwargs):
        if len(kwargs.get("fullbkpt", ())) == 100:      # the ledge fit
            return real_bspline.iterfit(xdata, ydata, *args, **kwargs)
        got.x = np.asarray(xdata, dtype=float).copy()
        got.y = np.asarray(ydata, dtype=float).copy()
        raise ReachedMasterFit()
    monkeypatch.setattr(make_master_flat_module, "Bspline",
                        SimpleNamespace(iterfit=iterfit))
    return got


def setup_flat(synthetic_ccddata, tmp_redux, step_list):
    rows = np.arange(NROW)
    wave_row = (LEDGE - 40.0 + DROW * rows).astype(np.float32)
    wave = np.repeat(wave_row[:, None], NCOL, axis=1)
    y = np.full(wave_row.shape, 10000.0)
    for w, h in step_list:
        y += h / (1.0 + np.exp(-(wave_row.astype(float) - w) / 0.3))
    flat = y.astype(np.float32)
    flat_img = np.repeat(flat[:, None], NCOL, axis=1)
    pos = np.repeat(np.arange(NCOL, dtype=np.float32)[None, :], NROW, axis=0)
    sl = np.full((NROW, NCOL), 9.0, dtype=np.float32)
    common = dict(BGRATNAM="BM", BCWAVE=BCWAVE, GROUPID="kb260101_00010")
    marc = synthetic_ccddata(IMTYPE="MARC", FRAMENO=3,
                             OFNAME="kb260101_00003.fits", **common)
    for suffix, data in (("wavemap", wave), ("slicemap", sl),
                         ("posmap", pos)):
        img = KCCDData(data, meta=marc.header.copy(), unit="adu")
        img.header["HISTORY"] = "synthetic " + suffix
        kcwi_fits_writer(img, output_file="kb260101_00003.fits",
                         output_dir=str(tmp_redux.redux), suffix=suffix)
    tmp_redux.proctab.update_proctab(frame=marc, suffix="marc",
                                     newtype="MARC",
                                     filename="kb260101_00003.fits")
    stack = synthetic_ccddata(IMTYPE="SDOME", FRAMENO=10,
                              OFNAME="kb260101_00010.fits", **common)
    img = KCCDData(flat_img, meta=stack.header.copy(), unit="electron")
    img.header["HISTORY"] = "synthetic stacked dome flat"
    kcwi_fits_writer(img, output_file="kb260101_00010.fits",
                     output_dir=str(tmp_redux.redux), suffix="sdome")
    tmp_redux.proctab.update_proctab(frame=stack, suffix="sdome",
                                     filename="kb260101_00010.fits")
    frame = synthetic_ccddata(IMTYPE="DOMEFLAT", FRAMENO=11,
                              OFNAME="kb260101_00011.fits", **common)
    return frame, dict(zip(wave_row.astype(float), flat.astype(float)))


def first_corrected(capture, original):
    """Wavelength of the first row the ledge correction changed."""
    orig = np.array([original[x] for x in capture.x])
    changed = capture.y != orig
    return capture.x[changed].min() if changed.any() else None


@pytest.mark.parametrize("typed, index", [
    ("40", 40), ("14", 14), ("030", 30), ("0", 0),
], ids=["40", "14", "leading-zero", "0"])
def test_ledge_peak_index_places_the_ledge(run_primitive, synthetic_ccddata,
                                           tmp_redux, capture, typed, index):
    frame, original = setup_flat(synthetic_ccddata, tmp_redux, [SMALL, BIG])
    with pytest.raises(ReachedMasterFit):
        run_primitive(MakeMasterFlat, frame, action_extra=EXTRA,
                      config_overrides={"plot_level": 2},
                      answers=["", typed, ""])
    # D3: ledge figure, peak, ledge-region figure
    assert len(run_primitive.last.prompts) == 3
    # source: analytic - P25: index n is at ledge - 15 + 0.5 n A; the
    # correction starts at the first 0.1 A row at or above it
    apk = LEDGE - 15.0 + 0.5 * index
    first = first_corrected(capture, original)
    assert first is not None
    assert apk <= first < apk + DROW + 1e-3


@pytest.mark.parametrize("bad", ["-3", "4.5", "abc", "+40"])
def test_ledge_peak_invalid_text_is_asked_again(run_primitive,
                                                synthetic_ccddata,
                                                tmp_redux, capture, bad):
    # source: analytic - P25: only digits are accepted; anything else is
    # asked again, with no effect on the reduction
    frame, original = setup_flat(synthetic_ccddata, tmp_redux, [SMALL, BIG])
    with pytest.raises(ReachedMasterFit):
        run_primitive(MakeMasterFlat, frame, action_extra=EXTRA,
                      config_overrides={"plot_level": 2},
                      answers=["", bad, "40", ""])
    assert len(run_primitive.last.prompts) == 4
    first = first_corrected(capture, original)
    assert LEDGE + 5.0 <= first < LEDGE + 5.0 + DROW + 1e-3


# ============================================ NandshuffSubtractSky (P37)
# Synthetic BLUE aborted (inverted) N&S frame, as test_run_Nandshuff...py
# (copied): sky rows 20-39, object rows 40-59.
SHUF = 20
NS = {"nasmask": True, "numopen": 2}
SKY, OBJ = 100.0, 350.0


def blue_frame(synthetic_ccddata):
    f = synthetic_ccddata(SHUFROWS=SHUF, NSHFUP=0, NSHFDN=0)
    f.data[:] = 0.0
    f.data[SHUF:2 * SHUF, :] = SKY
    f.data[2 * SHUF:3 * SHUF, :] = OBJ
    return f


def run_ns(run_primitive, frame, answers):
    return run_primitive(NandshuffSubtractSky, frame, action_extra=NS,
                         config_overrides={"plot_level": 2},
                         answers=answers).ccddata


@pytest.mark.parametrize("typed, scale", [
    ("1.5", 1.5), ("0.95", 0.95), ("1e0", 1.0), ("2", 2.0), (" 0.5 ", 0.5),
    ("3.", 3.0), (".25", 0.25),
], ids=["1.5", "0.95", "1e0", "2", "padded", "3.", ".25"])
def test_sky_scale_valid_float_scales_the_sky(run_primitive,
                                              synthetic_ccddata, typed,
                                              scale):
    out = run_ns(run_primitive, blue_frame(synthetic_ccddata), [typed])
    assert len(run_primitive.last.prompts) == 1
    assert out.header["NASSCL"] == scale        # source: analytic - P37
    # source: analytic - object panel minus scale * sky panel, here
    # 350 - scale * 100 (a flat sky of 100 and a flat object of 350).
    # Tolerance: float64 product and difference of two exactly-held values:
    # exact up to one rounding.
    np.testing.assert_allclose(
        out.data[SHUF + 1:2 * SHUF, :], OBJ - scale * SKY, rtol=1e-12)


def test_sky_scale_zero_is_accepted_until_b7_is_fixed(run_primitive,
                                                      synthetic_ccddata):
    # Valid today (B7 makes it invalid in PR3b: delete this test then, the
    # strict xfail in test_run_NandshuffSubtractSky.py takes over): 0 means
    # no sky subtraction.
    out = run_ns(run_primitive, blue_frame(synthetic_ccddata), ["0"])
    assert len(run_primitive.last.prompts) == 1
    assert out.header["NASSCL"] == 0.0
    assert np.all(out.data[SHUF + 1:2 * SHUF, :] == OBJ)   # source: analytic


@pytest.mark.parametrize("bad", ["-1", "abc", "1,5", "1 2"])
def test_sky_scale_invalid_is_asked_again(run_primitive, synthetic_ccddata,
                                          bad):
    # source: analytic - P37: negative numbers and text are refused
    out = run_ns(run_primitive, blue_frame(synthetic_ccddata), [bad, "0.5"])
    assert len(run_primitive.last.prompts) == 2
    assert out.header["NASSCL"] == 0.5


def test_red_branch_asks_for_the_scale_once(run_primitive,
                                            synthetic_ccddata):
    # P38: the red non-standard branch asks once (D3).  Its value cannot be
    # read back: today the branch raises IndexError right after the prompt
    # (B12, PR3b) and the NASSCL card is not yet written, so only the prompt
    # count is checked.  Add a NASSCL assertion in the commit that fixes B12.
    f = synthetic_ccddata(camera="RED", binning=(1, 1), SHUFROWS=SHUF,
                          NSHFUP=0, NSHFDN=1)
    f.data[:] = 0.0
    f.data[SHUF:2 * SHUF, :] = OBJ
    f.data[2 * SHUF:3 * SHUF, :] = SKY
    try:
        run_ns(run_primitive, f, ["1.25"])
    except IndexError:          # B12
        pass
    assert len(run_primitive.last.prompts) == 1
