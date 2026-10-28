"""MakeInvsens prompt sites run through ``run_primitive`` (task T.4).

Prompt sites (docs/plot_inventory.md, section 3.1):

* P13 (plot_level >= 3): per slice, ``Next? <cr>, q - quit: ``.
* P14 (plot_level > 1): wavelength limits, ``New values? <float> <float>
  (or <cr> - done): `` -- two floats inside [WAVALL0, WAVALL1] set the fit
  limits (INVFW0 / INVFW1), <cr> accepts.
* P15 (plot_level >= 2): line masks, ``Mask line: wavelength start stop
  (Ang) comment (str) <float> <float> <str> (<cr> - done): `` -- each
  accepted mask is written to the local ``<std>.lmsk`` file.
* P16-P18 (plot_level >= 2): ``Next? <cr>: `` after each of three figures.
* P19 (plot_level >= 2): ``Current fit order = %d, New fit order? <int>,
  <cr> - done: `` -- an integer refits with that polynomial order (INVFORD).

Known bugs pinned here as strict xfails (D8): B2 (one-number line mask),
B3 (recorded fit order after a failed fit), B5 (an answer containing "q"
quits, at P13), B6 (an invalid limits entry resets accepted limits), B13
(the local .lmsk header is written back without newlines) and one new bug
(B32, the shipped .lmsk files are never found).  The passing tests
characterize today's behaviour for valid answers (D3/D12).

Input: a synthetic DAR-corrected standard-star cube (the frame's OBJECT
names a standard shipped in data/stds), built in memory.  At plot_level >= 1
the primitive builds Bokeh figures (bokeh_plot returns at once, since the
harness's bokeh_session is None) and exports four PNGs with save_plot, which
needs a browser; save_plot is replaced by a no-op in this module's tests.
"""
import logging

import numpy as np
import pytest
from astropy.io import fits

from kcwidrp.primitives import MakeInvsens as make_invsens_module
from kcwidrp.primitives.MakeInvsens import MakeInvsens
from kcwidrp.primitives.kcwi_file_primitives import KCCDData

HEAD = "c0fed5c"
pytestmark = pytest.mark.filterwarnings(
    # order-7 polynomial over 500 A with x in Angstrom from the start
    "ignore:Polyfit may be poorly conditioned")

STD = "bd28d4211"
NW, NY, NSL = 600, 20, 5
WGOOD = (4050.0, 4550.0)       # WAVGOOD0/1 -> default fit limits
WALL = (4000.0, 4599.0)        # WAVALL0/1 -> allowed range for new limits
OUT = "kb260101_00001_invsens.fits"

# Exact prompt strings (docs/plot_inventory.md P13-P19)
P13 = "Next? <cr>, q - quit: "
P14 = "New values? <float> <float> (or <cr> - done): "
P15 = ("Mask line: wavelength start stop (Ang) comment (str) <float> "
       "<float> <str> (<cr> - done): ")
NEXT = "Next? <cr>: "


def P19(order):
    return ("Current fit order = %d, New fit order? <int>, <cr> - done: "
            % order)


@pytest.fixture(autouse=True)
def no_png_export(monkeypatch):
    # save_plot exports PNGs through selenium and a browser: not available
    # (and not wanted) in unit tests.
    monkeypatch.setattr(make_invsens_module, "save_plot",
                        lambda *a, **k: None)


def std_cube(synthetic_ccddata, std=STD):
    """A DAR-corrected cube with a star in slice 2 at y = 10."""
    frame = synthetic_ccddata(
        IMTYPE="OBJECT", OBJECT=std, TARGNAME=std, DARCOR=True,
        XPOSURE=100.0, CRVAL3=WALL[0], CD3_3=1.0, CRPIX3=1.0,
        WAVGOOD0=WGOOD[0], WAVGOOD1=WGOOD[1], WAVALL0=WALL[0],
        WAVALL1=WALL[1], DARPADY=2, SKYCOR=True, AIRMASS=1.0, DICHFRAC=1.0,
        BGRATNAM="BL")
    w = WALL[0] + np.arange(NW)
    spec = 1.0e4 * (1.0 + 0.3 * np.sin((w - 4000.0) / 150.0))
    yy = np.arange(NY)
    prof = np.exp(-0.5 * ((yy - 10.0) / 1.5) ** 2)
    cube = np.zeros((NW, NY, NSL))
    for sl, amp in ((1, 0.3), (2, 1.0), (3, 0.3)):
        cube[:, :, sl] = amp * spec[:, None] * prof[None, :]
    cube += 1.0
    ccd = KCCDData(cube, meta=frame.header, unit="electron",
                   uncertainty=np.sqrt(cube))
    return ccd


def run(run_primitive, synthetic_ccddata, plot_level, answers=None, std=STD):
    return run_primitive(MakeInvsens, std_cube(synthetic_ccddata, std),
                         config_overrides={"plot_level": plot_level},
                         answers=answers)


def out_header(tmp_redux):
    with fits.open(tmp_redux.redux / OUT) as hdul:
        return hdul[0].header.copy()


FIT_DONE = [NEXT, NEXT, NEXT]           # P16-P18 before each P19


# --------------------------------------------------------------- passing
def test_plot_level_0_uses_defaults_without_prompts(run_primitive,
                                                    synthetic_ccddata,
                                                    tmp_redux):
    run(run_primitive, synthetic_ccddata, 0)
    hdr = out_header(tmp_redux)
    assert run_primitive.last.prompts == []
    # source: characterization @ c0fed5c - fit limits WAVGOOD0/1, fit order
    # 7 for BL (times DICHFRAC = 1)
    assert (hdr["INVFW0"], hdr["INVFW1"]) == WGOOD
    assert hdr["INVFORD"] == 7


def test_valid_answers_set_limits_masks_and_order(run_primitive,
                                                  synthetic_ccddata,
                                                  tmp_redux):
    # source: characterization @ c0fed5c (inventory P14, P15, P19): valid
    # answers and the prompts they lead to
    answers = ["4100 4500", "",                  # P14: new limits, accept
               "4200 4210 Hbeta wing", "",       # P15: one mask, done
               "", "", "", "5",                  # P16-18, P19: order 5
               "", "", "", ""]                   # P16-18, P19: done
    run(run_primitive, synthetic_ccddata, 2, answers)
    assert run_primitive.last.prompts == (
        [P14, P14, P15, P15] + FIT_DONE + [P19(7)] + FIT_DONE + [P19(5)])
    hdr = out_header(tmp_redux)
    assert (hdr["INVFW0"], hdr["INVFW1"]) == (4100.0, 4500.0)
    assert hdr["INVFORD"] == 5
    # source: characterization @ c0fed5c - the mask is saved locally
    with open(tmp_redux.root / (STD + ".lmsk")) as f:
        assert f.read() == "4200.00 4210.00 # Hbeta wing\n"


def test_limits_outside_allowed_range_are_reset(run_primitive,
                                                synthetic_ccddata,
                                                tmp_redux):
    # source: characterization @ c0fed5c (inventory P14) - limits in the
    # wrong order are rejected and reset to WAVGOOD (no earlier accepted
    # limits, so this is not B6)
    answers = ["4500 4100", "", "", "", "", "", ""]
    run(run_primitive, synthetic_ccddata, 2, answers)
    hdr = out_header(tmp_redux)
    assert (hdr["INVFW0"], hdr["INVFW1"]) == WGOOD


def test_slice_review_quit_with_q(run_primitive, synthetic_ccddata,
                                  tmp_redux):
    # source: characterization @ c0fed5c (inventory P13) - "q" stops the
    # per-slice review; the remaining prompts follow
    answers = ["q", "", "", "", "", "", ""]
    run(run_primitive, synthetic_ccddata, 3, answers)
    assert run_primitive.last.prompts == [P13, P14, P15] + FIT_DONE + \
        [P19(7)]


def test_slice_review_shows_every_slice(run_primitive, synthetic_ccddata,
                                        tmp_redux):
    # source: characterization @ c0fed5c (inventory P13) - <cr> goes on
    answers = [""] * NSL + ["", "", "", "", "", ""]
    run(run_primitive, synthetic_ccddata, 3, answers)
    assert run_primitive.last.prompts == [P13] * NSL + [P14, P15] + \
        FIT_DONE + [P19(7)]


# ----------------------------------------------------------- known bugs
@pytest.mark.xfail(strict=True, raises=AssertionError, reason=(
    "B6: an invalid limits entry resets both limits to WAVGOOD, discarding "
    "the limits accepted earlier in the loop (MakeInvsens.py:369-379)"))
def test_invalid_limits_keep_accepted_limits(run_primitive,
                                             synthetic_ccddata, tmp_redux):
    # source: analytic @ c0fed5c - plan section 8.2 (ask_value keeps the
    # previously accepted value on invalid input)
    answers = ["4100 4500", "abc", "", "", "", "", "", ""]
    run(run_primitive, synthetic_ccddata, 2, answers)
    hdr = out_header(tmp_redux)
    assert (hdr["INVFW0"], hdr["INVFW1"]) == (4100.0, 4500.0)


@pytest.mark.xfail(strict=True, raises=IndexError, reason=(
    "B2: a line mask with only one number raises IndexError at "
    "qstr.split()[1]; only ValueError is caught (MakeInvsens.py:474-481)"))
def test_one_number_line_mask_is_asked_again(run_primitive,
                                             synthetic_ccddata, tmp_redux):
    # source: analytic @ c0fed5c - plan section 8.2 - invalid input re-prompts
    answers = ["", "4200", "", "", "", "", ""]
    run(run_primitive, synthetic_ccddata, 2, answers)
    assert run_primitive.last.prompts[:3] == [P14, P15, P15]
    assert (tmp_redux.redux / OUT).exists()


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=(
    "B3: after a failed refit (order -1: np.polyfit raises ValueError) the "
    "fit is the old order-7 one, but ford keeps -1, which is shown at the "
    "next prompt, logged and written to INVFORD (MakeInvsens.py:655-683)"))
def test_failed_refit_keeps_recorded_order(run_primitive, synthetic_ccddata,
                                           tmp_redux):
    # source: analytic @ c0fed5c - plan section 13 B3 - the recorded order must
    # be the order of the fit actually used (7, since the order -1 fit failed)
    answers = ["", "", "", "", "", "-1", "", "", "", ""]
    run(run_primitive, synthetic_ccddata, 2, answers)
    hdr = out_header(tmp_redux)
    assert run_primitive.last.prompts[-1] == P19(7)
    assert hdr["INVFORD"] == 7


@pytest.mark.xfail(strict=True, raises=pytest.fail.Exception, reason=(
    "B5: 'Q' in qstr.upper() quits the per-slice review on any answer "
    "containing q (MakeInvsens.py:196-198), so 'eq' ends the review and the "
    "answers meant for the remaining slices are left over"))
def test_answer_containing_q_does_not_quit(run_primitive, synthetic_ccddata,
                                           tmp_redux):
    # source: analytic @ c0fed5c - plan section 8.2 - quit is an exact match on
    # q/quit; other input containing a q re-prompts.  So "eq" is asked again,
    # then every slice is shown.
    answers = ["eq"] + [""] * NSL + ["", "", "", "", "", ""]
    run(run_primitive, synthetic_ccddata, 3, answers)
    assert run_primitive.last.prompts[:NSL + 1] == [P13] * (NSL + 1)


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=(
    "B13: header lines of an existing .lmsk file are read with strip() and "
    "written back without a newline (MakeInvsens.py:404-408, 488-493), so "
    "the header and the first mask become one comment line"))
def test_local_line_mask_file_survives_a_run(run_primitive,
                                             synthetic_ccddata, tmp_redux):
    lmsk = tmp_redux.root / (STD + ".lmsk")
    lmsk.write_text("# masks for %s\n4200 4210 # a\n4300 4310\n" % STD)
    run(run_primitive, synthetic_ccddata, 0)
    # source: analytic @ c0fed5c - plan section 13 B13 / the .lmsk format the
    # primitive reads: '#' lines are comments, other lines "w0 w1 [# comment]"
    masks = [ln.split("#")[0].split() for ln in lmsk.read_text().splitlines()
             if ln.strip() and not ln.strip().startswith("#")]
    assert [[float(v) for v in m] for m in masks] == [[4200.0, 4210.0],
                                                      [4300.0, 4310.0]]


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=(
    "B32: the shipped line-mask files (data/stds/bd02d3375.lmsk, "
    "feige110.lmsk) are never found: the 'then in stds directory' lookup "
    "uses strip_fname(stdfile) + '.lmsk', and strip_fname returns the bare "
    "stem, so it looks in the current directory again (MakeInvsens.py:391, "
    "since eaa516c)"))
def test_shipped_line_mask_file_is_used(run_primitive, synthetic_ccddata,
                                        tmp_redux, caplog):
    # source: analytic @ c0fed5c - code comment and docstring intent at
    # MakeInvsens.py:386-395 - "first in local directory, then in stds
    # directory"
    caplog.set_level(logging.INFO, logger="KCWI")
    run(run_primitive, synthetic_ccddata, 0, std="bd02d3375")
    assert any("Using line mask file" in r.getMessage()
               for r in caplog.records)


# ------------------------------------------------- physical scale (T.8)
# Added by T.8 (suite quality check) in a fresh session, after a deliberate
# breakage of the photon-energy constant in MakeInvsens went unnoticed: the
# prompt tests above never look at the numbers in the output files.
#
# Specification (class docstring, units of the output files): the standard's
# reference spectrum F(lambda) (erg / s / cm^2 / A, data/stds) is divided by
# the observed spectrum in e- / s (summed over the star, per bin of dw A) to
# give the inverse sensitivity (row 0 of ``*_invsens.fits``, in erg / s /
# cm^2 / A per e- / s per bin); the effective area (row 0 of ``*_ea.fits``,
# cm^2 / A) is the observed rate divided by the photon rate the star delivers
# per cm^2: F * lambda / (h c) * dw photons / s / cm^2.  The Keck telescope
# area is 760000 cm^2, so a system of throughput eta observes
#     counts / s = eta * 760000 * F * lambda / (h c) * dw
# and the effective area must come out as eta * 760000.
AREA_KECK = 760000.0                    # cm^2 (plan/primitive: Keck)


def throughput_cube(synthetic_ccddata, eta, expt=100.0):
    from astropy import constants as const
    from astropy import units as u
    from kcwidrp.core.kcwi_pkg_resources import get_resource_path
    frame = synthetic_ccddata(
        IMTYPE="OBJECT", OBJECT=STD, TARGNAME=STD, DARCOR=True,
        XPOSURE=expt, CRVAL3=WALL[0], CD3_3=1.0, CRPIX3=1.0,
        WAVGOOD0=WGOOD[0], WAVGOOD1=WGOOD[1], WAVALL0=WALL[0],
        WAVALL1=WALL[1], DARPADY=2, SKYCOR=True, AIRMASS=0.0, DICHFRAC=1.0,
        BGRATNAM="BL")
    w = WALL[0] + np.arange(NW)
    with fits.open(get_resource_path("kcwidrp", "data/stds/%s.fits" % STD)
                   ) as hdul:
        swl = hdul[1].data["WAVELENGTH"]
        sflx = hdul[1].data["FLUX"]
    flux = np.interp(w, swl, sflx)
    # source: independent - astropy.constants (CODATA) for h c
    hc = (const.h * const.c).to(u.erg * u.AA).value
    rate = eta * AREA_KECK * flux * w / hc * 1.0       # dw = 1 A
    cube = np.zeros((NW, NY, NSL))
    cube[:, 10, 2] = rate * expt
    return KCCDData(cube, meta=frame.header, unit="electron",
                    uncertainty=np.sqrt(cube + 1.0)), w, flux, rate, hc


def test_effective_area_is_throughput_times_telescope_area(
        run_primitive, synthetic_ccddata, tmp_redux):
    eta = 0.30
    ccd, w, flux, rate, hc = throughput_cube(synthetic_ccddata, eta)
    run_primitive(MakeInvsens, ccd, config_overrides={"plot_level": 0})
    with fits.open(tmp_redux.redux / "kb260101_00001_ea.fits") as hdul:
        fearea = hdul[0].data[1].copy()
        earea = hdul[0].data[0].copy()
    good = (w > 4100.0) & (w < 4500.0)
    # source: analytic - see the block comment above: eta * 760000 cm^2
    # Tolerance: 1 %.  The star's spectrum is built with a linear
    # interpolation of the reference table and the primitive uses a cubic
    # one (largest difference measured: 6e-4 of the value, plus float32
    # storage); the fitted row is smoother still (5e-6).  1 % is 17 times the
    # spread and far below the factor-10 error of a wrong photon constant.
    np.testing.assert_allclose(earea[good], eta * AREA_KECK, rtol=0.01)
    np.testing.assert_allclose(fearea[good], eta * AREA_KECK, rtol=0.01)


def test_invsens_is_reference_flux_over_observed_rate(
        run_primitive, synthetic_ccddata, tmp_redux):
    ccd, w, flux, rate, hc = throughput_cube(synthetic_ccddata, 0.30)
    run_primitive(MakeInvsens, ccd, config_overrides={"plot_level": 0})
    with fits.open(tmp_redux.redux / "kb260101_00001_invsens.fits") as hdul:
        invsens = hdul[0].data[0].copy()
        fitted = hdul[0].data[1].copy()
        observed = hdul[0].data[2].copy()
    good = (w > 4100.0) & (w < 4500.0)
    # source: analytic - observed spectrum in e-/s per bin is the rate put
    # into the cube; the inverse sensitivity is F / observed.  Tolerance as
    # above (1 %).
    np.testing.assert_allclose(observed[good], rate[good], rtol=1e-6)
    np.testing.assert_allclose(invsens[good], flux[good] / rate[good],
                               rtol=0.01)
    np.testing.assert_allclose(fitted[good], flux[good] / rate[good],
                               rtol=0.01)
