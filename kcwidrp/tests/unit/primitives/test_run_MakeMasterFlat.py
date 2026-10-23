"""MakeMasterFlat BM-ledge peak selection run through ``run_primitive`` (T.4).

(The T.3 file ``test_MakeMasterFlat.py`` covers bm_ledge_position.)

Specification (class docstring, plan section 13 B1, inventory P24/P25): for
the BM grating the flat shows a "ledge" (a step in the spectrum of the
reference slice) near ``bm_ledge_position(cwave, dich)``.  The primitive
smooths the reference-slice spectrum within +-25 A of that wavelength, fits
a B-spline, takes its derivative on a 0.5 A grid, and finds the peaks of the
derivative.  One peak: that is the ledge.  Several or none: at plot_level
>= 2 the user is asked

    ``Peak? <int>: ``   (an index into the trimmed grid, xvals[n] =
                        ledge - 15 + 0.5 n A)

after the ledge figure's ``Next? <cr>: ``; below plot_level 2 the highest
peak is meant to be chosen automatically.  Data at wavelengths >= the chosen
ledge position ``apk`` are then divided by the step ratio measured across
``apk`` (fits at apk-8 ... apk-3 and apk+3 ... apk+8 A), which removes the
ledge from the master-flat input.

What is measured: the reference-slice data handed to the master-flat
B-spline fit, captured at that call (the first Bspline.iterfit call with
other than the ledge fit's 100 breakpoints).  The capture stops the run
there (``ReachedMasterFit``); everything after it (blue/red slice fits, the
ratio flat) is outside the prompt site and not run.  Comparing the captured
data with the input flat shows where the correction starts (apk, to the
0.1 A row spacing) and that the ledge is gone.

Input: a synthetic BLUE BM dome flat (internal flats would add the
vignetting fits, which are not needed to reach the site) with wavelength,
slice and position maps: 800 rows at 0.1 A / row around the ledge, 40
positions, all in the reference slice 9.  The flat is 10000 e- with
logistic steps (width 0.3 A) at chosen wavelengths.  Figures are built at
every plot level (bokeh_plot returns at once: bokeh_session is None);
save_plot is replaced by a no-op (PNG export through a browser) although
the run stops before its first call.
"""
from types import SimpleNamespace

import numpy as np
import pytest

from kcwidrp.core.bspline import Bspline as real_bspline
from kcwidrp.primitives import MakeMasterFlat as make_master_flat_module
from kcwidrp.primitives.MakeMasterFlat import MakeMasterFlat, \
    bm_ledge_position
from kcwidrp.primitives.kcwi_file_primitives import KCCDData, \
    kcwi_fits_writer

HEAD = "c0fed5c"
NROW, NCOL, DROW = 800, 40, 0.1
BCWAVE = 4500.0
# source: characterization @ b0e2550 (T.3) - the documented ledge formula;
# the fixture's header has RCWAVE, so the dichroic is in
LEDGE = bm_ledge_position(BCWAVE, True)          # 5143.899 A
EXTRA = {"new_type": "MDOME", "stack_type": "SDOME", "min_files": 1}
NEXT = "Next? <cr>: "
PEAK = "Peak? <int>: "


class ReachedMasterFit(Exception):
    """Raised at the master-flat B-spline fit, after the ledge block."""


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
    # PNG export needs a browser; not part of what is tested
    monkeypatch.setattr(make_master_flat_module, "save_plot",
                        lambda *a, **k: None)
    return got


def steps(wave, step_list):
    """10000 e- plus logistic steps [(wavelength, height), ...]."""
    y = np.full(wave.shape, 10000.0)
    for w, h in step_list:
        y += h / (1.0 + np.exp(-(wave - w) / 0.3))
    return y


def setup_flat(synthetic_ccddata, tmp_redux, step_list):
    """Write maps, the stacked flat and proc-table rows; return
    (frame to run on, {row wavelength: flat value})."""
    rows = np.arange(NROW)
    wave_row = (LEDGE - 40.0 + DROW * rows).astype(np.float32)
    wave = np.repeat(wave_row[:, None], NCOL, axis=1)
    flat = steps(wave_row.astype(float), step_list).astype(np.float32)
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


def run_to_master_fit(run_primitive, frame, plot_level, answers=None):
    with pytest.raises(ReachedMasterFit):
        run_primitive(MakeMasterFlat, frame, action_extra=EXTRA,
                      config_overrides={"plot_level": plot_level},
                      answers=answers)


def correction(capture, original):
    """(first corrected wavelength, ratio applied there) from the capture."""
    orig = np.array([original[x] for x in capture.x])
    changed = capture.y != orig
    if not changed.any():
        return None, 1.0
    first = capture.x[changed].min()
    return first, float(np.median(orig[changed] / capture.y[changed]))


def ledge_left(capture, at):
    """Step ratio still present across ``at`` (mean 3-8 A above / below)."""
    up = (capture.x >= at + 3) & (capture.x <= at + 8)
    lo = (capture.x >= at - 8) & (capture.x <= at - 3)
    return capture.y[up].mean() / capture.y[lo].mean()


BIG = (LEDGE + 5.0, 2000.0)       # the ledge: +20 %
SMALL = (LEDGE - 8.0, 300.0)      # a second, lower derivative peak


# --------------------------------------------------------------- passing
def test_single_peak_ledge_is_removed(run_primitive, synthetic_ccddata,
                                      tmp_redux, capture):
    frame, original = setup_flat(synthetic_ccddata, tmp_redux, [BIG])
    run_to_master_fit(run_primitive, frame, 1)
    first, ratio = correction(capture, original)
    # source: analytic @ c0fed5c - the derivative of a single step peaks at the
    # step; apk lies on the 0.5 A grid and the correction starts at the first
    # row (0.1 A spacing) at or above it: within 0.6 A of the step.
    assert abs(first - BIG[0]) <= 0.6
    # source: analytic @ c0fed5c - the step is 12000 / 10000; after the
    # correction no step is left.  Tolerance: the fits' and smoothing's
    # residuals, 1 %.
    assert ratio == pytest.approx(1.2, rel=0.01)
    assert ledge_left(capture, BIG[0]) == pytest.approx(1.0, rel=0.01)
    assert run_primitive.last.prompts == []


def test_interactive_peak_choice(run_primitive, synthetic_ccddata,
                                 tmp_redux, capture):
    # source: characterization @ c0fed5c (inventory P24, P25): with two
    # peaks at plot_level 2 the ledge figure's Next?, then Peak?, then the
    # ledge-region figure's Next?; answer n puts apk at xvals[n] =
    # LEDGE - 15 + 0.5 n, so 40 -> LEDGE + 5
    frame, original = setup_flat(synthetic_ccddata, tmp_redux, [SMALL, BIG])
    run_to_master_fit(run_primitive, frame, 2, ["", "40", ""])
    assert run_primitive.last.prompts == [NEXT, PEAK, NEXT]
    first, _ = correction(capture, original)
    # apk = LEDGE + 5 exactly; first row at or above it (0.1 A rows)
    assert LEDGE + 5.0 <= first < LEDGE + 5.0 + DROW + 1e-3
    assert ledge_left(capture, BIG[0]) == pytest.approx(1.0, rel=0.01)


def test_interactive_peak_reasks_non_numeric(run_primitive,
                                             synthetic_ccddata, tmp_redux,
                                             capture):
    # source: characterization @ c0fed5c (inventory P25): empty, negative
    # and decimal answers are asked again
    frame, _ = setup_flat(synthetic_ccddata, tmp_redux, [SMALL, BIG])
    run_to_master_fit(run_primitive, frame, 2,
                      ["", "", "-3", "4.5", "40", ""])
    assert run_primitive.last.prompts == [NEXT] + [PEAK] * 4 + [NEXT]


# ----------------------------------------------------------- known bugs
def test_several_peaks_highest_is_chosen(run_primitive, synthetic_ccddata,
                                         tmp_redux, capture):
    frame, original = setup_flat(synthetic_ccddata, tmp_redux, [SMALL, BIG])
    run_to_master_fit(run_primitive, frame, 1)
    first, _ = correction(capture, original)
    # source: analytic @ c0fed5c - plan section 13 B1 - the highest derivative
    # peak is the big step; tolerance as in test_single_peak_ledge_is_removed
    assert first is not None and abs(first - BIG[0]) <= 0.6
    assert ledge_left(capture, BIG[0]) == pytest.approx(1.0, rel=0.01)


def test_no_peak_does_not_fail(run_primitive, synthetic_ccddata, tmp_redux,
                               capture):
    # source: analytic @ c0fed5c - plan section 13 B1b - the primitive must
    # carry on to the master-flat fit; with no ledge peak the maintainer chose
    # (2026-10-23, 1.3.4) to skip the ledge correction, so nothing is corrected
    frame, original = setup_flat(synthetic_ccddata, tmp_redux, [])
    run_to_master_fit(run_primitive, frame, 1)
    assert capture.x is not None
    first, _ = correction(capture, original)
    assert first is None


@pytest.mark.parametrize("bad", [
    pytest.param("999", marks=pytest.mark.xfail(
        strict=True, raises=IndexError, reason=(
            "B1c: the peak index is not bounds-checked; 999 >= len(xvals) "
            "raises IndexError at apk = xvals[ipk] "
            "(MakeMasterFlat.py:469-479)"))),
    pytest.param("\u00b2", marks=pytest.mark.xfail(
        strict=True, raises=ValueError, reason=(
            "B1c: '\u00b2' (superscript two).isnumeric() is True, so the "
            "loop accepts it and int() raises ValueError "
            "(MakeMasterFlat.py:471-474)"))),
], ids=["out-of-range", "superscript-two"])
def test_invalid_peak_index_is_asked_again(run_primitive, synthetic_ccddata,
                                           tmp_redux, capture, bad):
    # source: analytic @ c0fed5c - plan section 8.2 / 13 B1c - invalid input
    # re-prompts; then 40 (LEDGE + 5) is used
    frame, original = setup_flat(synthetic_ccddata, tmp_redux, [SMALL, BIG])
    run_to_master_fit(run_primitive, frame, 2, ["", bad, "40", ""])
    assert run_primitive.last.prompts == [NEXT, PEAK, PEAK, NEXT]
    first, _ = correction(capture, original)
    assert LEDGE + 5.0 <= first < LEDGE + 5.0 + DROW + 1e-3
