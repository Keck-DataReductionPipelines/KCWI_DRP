"""FluxCalibrate run through ``run_primitive`` (task T.4, numeric case).

Specification (class docstring, plan T.4): the inverse sensitivity written
by MakeInvsens (``redux/<std frame>_invsens.fits``, primary array of three
rows: raw invsens, fitted invsens, observed spectrum; wavelength axis
CRVAL1 + i * CDELT1) is found through the proc table (TYPE INVSENS, same
STATEID; of several, the nearest in MJD).  Its fitted row (row 1), in
erg / s / cm^2 / A per (e- / s), multiplies every spaxel's spectrum
divided by the exposure time XPOSURE; the result is stored in units of
1e-16 erg / s / cm^2 / A ("flam16"), so

    out = in * invsens * 1e16 / XPOSURE * 10**(0.4 * k(lambda) * AIRMASS)

where the last factor is the atmospheric extinction correction with the
extinction curve k in mag / airmass (``kcwi_correct_extin``; with AIRMASS =
0 it is exactly 1).  STDCOR, MSFILE and MSIMNO record the calibration.

The invsens file is built here (as MakeInvsens writes it) and entered in the
proc table, so the lookup is the real one.
"""
import logging

import numpy as np
import pytest
from astropy.io import fits

from kcwidrp.core.kcwi_pkg_resources import get_resource_path
from kcwidrp.primitives.FluxCalibrate import FluxCalibrate
from kcwidrp.primitives.kcwi_file_primitives import KCCDData

HEAD = "c0fed5c"
NW, NY, NX = 50, 3, 4
CRVAL3, CD3_3 = 4000.0, 2.0
STD = "kb260101_00005"


def invsens_curve():
    # A smooth, positive curve (arbitrary; values typical of KCWI BL,
    # ~1e-17 erg/s/cm^2/A per e-/s)
    w = CRVAL3 + CD3_3 * np.arange(NW)
    return 1.0e-17 * (1.0 + 0.5 * np.sin((w - CRVAL3) / 30.0))


def register_invsens(synthetic_ccddata, tmp_redux, curve, *, frameno=5,
                     stem=STD, mjd=61041.0, crval1=CRVAL3, cdelt1=CD3_3):
    """Write redux/<stem>_invsens.fits and add its proc-table row."""
    data = np.vstack([curve * 1.1, curve, np.ones_like(curve)])
    hdr = fits.Header({"CRVAL1": crval1, "CDELT1": cdelt1, "CRPIX1": 1.0,
                       "FRAMENO": frameno})
    fits.PrimaryHDU(data, hdr).writeto(
        tmp_redux.redux / (stem + "_invsens.fits"))
    std = synthetic_ccddata(FRAMENO=frameno, OFNAME=stem + ".fits", MJD=mjd)
    tmp_redux.proctab.update_proctab(frame=std, suffix="invsens",
                                     newtype="INVSENS",
                                     filename=stem + ".fits")
    return stem + "_invsens.fits"


def object_cube(synthetic_ccddata, **hdr):
    frame = synthetic_ccddata(CRVAL3=CRVAL3, CD3_3=CD3_3, CRPIX3=1.0,
                              FRAMENO=1, **hdr)
    rng = np.random.default_rng(3)
    cube = rng.uniform(100.0, 1000.0, size=(NW, NY, NX))
    unc = np.sqrt(cube)
    ccd = KCCDData(cube, meta=frame.header, unit="electron",
                   uncertainty=unc)
    return ccd, cube.copy(), unc.copy()


def test_output_is_input_times_invsens_over_exposure_time(
        run_primitive, synthetic_ccddata, tmp_redux):
    curve = invsens_curve()
    name = register_invsens(synthetic_ccddata, tmp_redux, curve)
    ccd, cube, unc = object_cube(synthetic_ccddata, XPOSURE=300.0,
                                 AIRMASS=0.0)
    out = run_primitive(FluxCalibrate, ccd).ccddata
    # source: analytic @ c0fed5c (docstring formula; AIRMASS 0 -> extinction
    # factor 1)
    scale = curve * 1.0e16 / 300.0
    expected = cube * scale[:, None, None]
    # Tolerance: the output is stored as float32 by the writer (1 ulp =
    # 1.2e-7 relative); the float64 products agree to ~1e-16.
    np.testing.assert_allclose(out.data, expected, rtol=1.2e-7, atol=0)
    np.testing.assert_allclose(out.uncertainty.array,
                               unc * scale[:, None, None], rtol=1.2e-7,
                               atol=0)
    hdr = out.header
    # source: analytic @ c0fed5c - docstring
    assert hdr["STDCOR"] is True
    assert hdr["MSFILE"] == name
    assert hdr["MSIMNO"] == 5
    assert out.unit.to_string() == "1e+16 erg / (Angstrom s cm2)"
    assert (tmp_redux.redux / "kb260101_00001_icubes.fits").exists()


def test_extinction_correction_uses_airmass(run_primitive, synthetic_ccddata,
                                            tmp_redux):
    curve = invsens_curve()
    register_invsens(synthetic_ccddata, tmp_redux, curve)
    ccd, cube, _ = object_cube(synthetic_ccddata, XPOSURE=100.0,
                               AIRMASS=1.3)
    out = run_primitive(FluxCalibrate, ccd).ccddata
    # source: published @ c0fed5c - SNfactory Mauna Kea extinction curve (Buton
    # et al. 2013, A&A 549, A8) as shipped in data/extin/snfext.fits
    # (mag/airmass on a 2 A grid).  The cube's planes (4000 + 2j A) are table
    # nodes, so the curve is read off the table directly; mag -> flux 10**(0.4
    # k X).
    with fits.open(get_resource_path("kcwidrp",
                                     "data/extin/snfext.fits")) as hdul:
        lam = np.asarray(hdul[1].data["LAMBDA"], dtype=float)
        ext = np.asarray(hdul[1].data["EXT"], dtype=float)
    w = CRVAL3 + CD3_3 * np.arange(NW)
    k = np.array([ext[np.flatnonzero(lam == wl)[0]] for wl in w])
    expected = cube * (curve * 1.0e16 / 100.0 *
                       10.0 ** (0.4 * k * 1.3))[:, None, None]
    # Tolerance: float32 storage (1.2e-7) plus the cubic spline through the
    # table evaluated at its own nodes (exact up to ~1e-15 relative).
    np.testing.assert_allclose(out.data, expected, rtol=2e-7, atol=0)
    assert out.header["EXTCOR"] is True


def test_zero_exposure_time_falls_back_to_ttime(run_primitive,
                                                synthetic_ccddata, tmp_redux):
    # source: characterization @ c0fed5c - XPOSURE <= 0 -> TTIME is used
    curve = invsens_curve()
    register_invsens(synthetic_ccddata, tmp_redux, curve)
    ccd, cube, _ = object_cube(synthetic_ccddata, XPOSURE=0.0, TTIME=50.0,
                               AIRMASS=0.0)
    out = run_primitive(FluxCalibrate, ccd).ccddata
    # Tolerance: float32 storage, 1 ulp
    np.testing.assert_allclose(
        out.data, cube * (curve * 1.0e16 / 50.0)[:, None, None],
        rtol=1.2e-7, atol=0)


def test_nearest_invsens_in_mjd_is_used(run_primitive, synthetic_ccddata,
                                        tmp_redux):
    # source: analytic @ c0fed5c - docstring of the proc-table search - nearest
    # in MJD
    curve = invsens_curve()
    register_invsens(synthetic_ccddata, tmp_redux, 3.0 * curve, frameno=5,
                     stem="kb260101_00005", mjd=61040.0)
    near = register_invsens(synthetic_ccddata, tmp_redux, curve, frameno=6,
                            stem="kb260101_00006", mjd=61041.1)
    ccd, cube, _ = object_cube(synthetic_ccddata, XPOSURE=10.0, AIRMASS=0.0,
                               MJD=61041.2)
    out = run_primitive(FluxCalibrate, ccd).ccddata
    assert out.header["MSFILE"] == near
    np.testing.assert_allclose(
        out.data, cube * (curve * 1.0e16 / 10.0)[:, None, None],
        rtol=1.2e-7, atol=0)


def test_invsens_on_other_grid_is_resampled(run_primitive, synthetic_ccddata,
                                            tmp_redux, caplog):
    # The resampling branch logs a malformed INFO record (B31, see
    # below), which pytest's log handler turns into a TypeError; keep the
    # KCWI logger at WARNING so this test checks the arithmetic whatever
    # level an earlier test left behind.
    caplog.set_level(logging.WARNING, logger="KCWI")
    # source: analytic @ c0fed5c - an invsens tabulated on a grid offset by
    # half a pixel is interpolated (cubic) onto the cube's planes; for a linear
    # curve the cubic interpolant is exact.
    w_cube = CRVAL3 + CD3_3 * np.arange(NW)
    w_inv = (CRVAL3 - 1.0) + CD3_3 * np.arange(NW + 1)
    lin = lambda w: 1.0e-17 * (1.0 + (w - CRVAL3) / 1000.0)  # noqa: E731
    register_invsens(synthetic_ccddata, tmp_redux, lin(w_inv),
                     crval1=CRVAL3 - 1.0)
    ccd, cube, _ = object_cube(synthetic_ccddata, XPOSURE=1.0, AIRMASS=0.0)
    out = run_primitive(FluxCalibrate, ccd).ccddata
    # Tolerance: float32 storage plus cubic-spline rounding (~1e-15)
    np.testing.assert_allclose(
        out.data, cube * (lin(w_cube) * 1.0e16)[:, None, None], rtol=2e-7,
        atol=0)


def test_no_invsens_leaves_cube_uncalibrated(run_primitive,
                                             synthetic_ccddata):
    # source: analytic @ c0fed5c - docstring - without an INVSENS row the
    # precondition fails
    ccd, cube, _ = object_cube(synthetic_ccddata, XPOSURE=1.0, AIRMASS=0.0)
    run_primitive(FluxCalibrate, ccd)
    assert np.array_equal(ccd.data, cube)
    assert "STDCOR" not in ccd.header or ccd.header["STDCOR"] is not True



@pytest.mark.xfail(strict=True, raises=TypeError, reason=(
    "B31: FluxCalibrate.py:98 calls self.logger.info(w0, mcw0, dw, "
    "...) with a float as the message and seven extra arguments, so "
    "formatting the record fails ('not all arguments converted'); in a "
    "reduction logging at INFO this prints a '--- Logging error ---' "
    "traceback whenever the invsens grid differs from the cube's"))
def test_resampling_message_is_logged(run_primitive, synthetic_ccddata,
                                      tmp_redux, caplog):
    # source: analytic @ c0fed5c - docstring of logging.Logger.info - msg %
    # args must format. pytest's log handler re-raises a formatting error
    # instead of printing it, which is what makes the defect visible here.
    caplog.set_level(logging.INFO)
    w_inv = (CRVAL3 - 1.0) + CD3_3 * np.arange(NW + 1)
    register_invsens(synthetic_ccddata, tmp_redux, 1.0e-17 + 0.0 * w_inv,
                     crval1=CRVAL3 - 1.0)
    ccd, _, _ = object_cube(synthetic_ccddata, XPOSURE=1.0, AIRMASS=0.0)
    run_primitive(FluxCalibrate, ccd)
    assert any("resampling standard" in r.getMessage()
               for r in caplog.records)
