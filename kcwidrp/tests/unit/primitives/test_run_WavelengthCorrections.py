"""WavelengthCorrections run through ``run_primitive`` (task T.4, numeric).

(The T.3 file ``test_WavelengthCorrections.py`` covers a2v_conversion, the
velocity value written to VCORR and get_wav_axis on their own.)

Specification (class docstring, plan T.4): the primitive reads the
``redux/<frame>_icube.fits`` cube, applies the radial velocity correction of
``radial_velocity_correction`` (default "heliocentric") and, with
``air_to_vacuum`` (default True), the air-to-vacuum conversion, and writes
``redux/<frame>_icubew.fits``.  The cube keeps its wavelength grid (CRVAL3,
CD3_3, CRPIX3 = 1: plane j is at CRVAL3 + j * CD3_3, FITS WCS paper I);
the data are resampled so that light observed at air wavelength lambda_air
lands at its corrected vacuum wavelength

    lambda = n(lambda_air) * lambda_air * (1 + v/c),

with v the correction velocity (astropy ``SkyCoord.radial_velocity_correction``
for the target, DATE-BEG and the Keck location in kcwi.cfg; astropy's
documented convention: true velocity = measured + v + measured * v / c, i.e.
1 + z_true = (1 + z_measured)(1 + v/c)) and n the refractive index of
standard air.  CTYPE3 becomes WAVE, VCORR records v.

Expected values (independent of the primitive's arithmetic):
* air -> vacuum for [O III] 5006.843 A: the published vacuum wavelength
  5008.239 A (SDSS / NIST line lists, Morton 1991 IAU conversion);
* air -> vacuum on a grid: Ciddor (1996, Appl. Opt. 35, 1566) eq. 1 for
  standard air, inverted here by fixed-point iteration;
* v: astropy (the same coordinates, time and location).
"""
import numpy as np
import pytest
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits

from kcwidrp.primitives.kcwi_file_primitives import KCCDData, \
    kcwi_fits_writer
from kcwidrp.primitives.WavelengthCorrections import WavelengthCorrections

HEAD = "c0fed5c"
C_KMS = 299792.458                 # km/s, exact (SI)
NAME = "kb260101_00001"            # OFNAME stem of the synthetic frame
TARGET = dict(TARGRA=150.1, TARGDEC=2.2, **{"DATE-BEG":
                                            "2026-01-15T10:00:00.000"})
# source: published @ c0fed5c (see docstring)
O3_AIR, O3_VAC = 5006.843, 5008.239


def n_ciddor(lam_air):
    """Refractive index of standard air, Ciddor (1996) eq. 1.

    # source: published @ c0fed5c - Ciddor (1996) eq. 1 (15 C, 101325 Pa, dry,
    # 450 ppm CO2); sigma in inverse micrometres.
    """
    s2 = (1.0e4 / np.asarray(lam_air, dtype=float)) ** 2
    return 1.0 + 1.0e-8 * (5792105.0 / (238.0185 - s2) +
                           167917.0 / (57.362 - s2))


def vac_to_air(lam_vac):
    """Invert lambda_vac = n(lambda_air) * lambda_air by iteration."""
    lam = np.asarray(lam_vac, dtype=float)
    air = lam / 1.000277
    for _ in range(20):           # converges to < 1e-12 A in ~4 steps
        air = lam / n_ciddor(air)
    return air


def velocity(kcwi_config):
    # source: independent @ c0fed5c - astropy, Keck location from kcwi.cfg
    cfg = kcwi_config()
    keck = EarthLocation.from_geodetic(lat=cfg.latitude, lon=cfg.longitude,
                                       height=cfg.altitude)
    target = SkyCoord(TARGET["TARGRA"], TARGET["TARGDEC"], unit="deg",
                      obstime=TARGET["DATE-BEG"])
    return target.radial_velocity_correction(
        kind="heliocentric", location=keck).to_value("km/s")


def write_icube(synthetic_ccddata, tmp_redux, spectrum, crval3, cd3_3):
    """Write redux/<NAME>_icube.fits holding ``spectrum`` in 2x2 spaxels."""
    frame = synthetic_ccddata(CRVAL3=crval3, CD3_3=cd3_3, CRPIX3=1.0,
                              CTYPE3="AWAV", **TARGET)
    del frame.header["VCORR"]           # a fresh cube has no correction yet
    cube = np.repeat(np.asarray(spectrum, dtype=float)[:, None, None], 4,
                     axis=1).reshape(len(spectrum), 2, 2)
    ccd = KCCDData(cube, meta=frame.header, unit="electron")
    ccd.header["HISTORY"] = "synthetic icube"
    kcwi_fits_writer(ccd, output_file=NAME + ".fits",
                     output_dir=str(tmp_redux.redux), suffix="icube")
    return frame


def run(run_primitive, frame, **config):
    out = run_primitive(WavelengthCorrections, frame,
                        config_overrides=config or None)
    return out.ccddata


def grid(cube_hdr):
    # source: published @ c0fed5c - FITS WCS paper I, CRPIX3 = 1
    n = cube_hdr["NAXIS3"]
    return cube_hdr["CRVAL3"] + cube_hdr["CD3_3"] * np.arange(n)


def test_line_lands_at_vacuum_velocity_corrected_wavelength(
        run_primitive, synthetic_ccddata, tmp_redux, kcwi_config):
    crval, cd, nw = 4900.0, 0.5, 440
    w = crval + cd * np.arange(nw)
    sigma = 1.5                                          # Angstrom
    line = 1000.0 * np.exp(-0.5 * ((w - O3_AIR) / sigma) ** 2)
    frame = write_icube(synthetic_ccddata, tmp_redux, line, crval, cd)
    out = run(run_primitive, frame)
    spec = np.asarray(out.data[:, 0, 0], dtype=np.float64)
    wave = grid(out.header)
    win = np.abs(wave - O3_VAC) < 10 * sigma
    centre = np.sum(wave[win] * spec[win]) / np.sum(spec[win])
    v = velocity(kcwi_config)
    # source: published @ c0fed5c (O3_VAC) and independent (v, astropy)
    expected = O3_VAC * (1.0 + v / C_KMS)
    # Tolerance: 0.01 A.  The published wavelengths are rounded to 0.001 A
    # and the Morton (1991) conversion differs from Ciddor's (which the DRP
    # uses) by < 0.002 A here; the first moment of a 3-px-sigma Gaussian
    # after two cubic resamplings is good to ~1e-4 A.  Insensitive to B19
    # on purpose: B19 moves the line by about CD3_3 * ((n - 1) +
    # lambda dn/dlambda + v/c) ~ 2e-4 A, see the strict xfail below.
    assert centre == pytest.approx(expected, abs=0.01)
    # The correction is not negligible compared with the tolerance:
    # source: analytic @ c0fed5c - vacuum shift 1.4 A, velocity shift |v|/c *
    # 5008 A
    assert abs(expected - O3_AIR) > 1.0
    hdr = out.header
    # source: analytic @ c0fed5c - docstring
    assert hdr["CTYPE3"] == "WAVE"
    # source: analytic @ c0fed5c - docstring
    assert hdr["VCORRTYP"] == "heliocentric"
    # Tolerance: same astropy call, 1 m/s (as T.3)
    assert hdr["VCORR"] == pytest.approx(v, abs=1e-3)
    # the grid itself is not changed
    assert (hdr["CRVAL3"], hdr["CD3_3"]) == (crval, cd)
    assert (tmp_redux.redux / (NAME + "_icubew.fits")).exists()


def ramp_case(run_primitive, synthetic_ccddata, tmp_redux, kcwi_config):
    """Spectrum whose value is the observed air wavelength of each pixel.

    Cubic interpolation reproduces a (nearly) linear function, so after the
    corrections every output pixel holds the observed air wavelength of the
    light that landed in it.  Returns (got, expected) for the inner pixels.
    """
    crval, cd, nw = 4900.0, 0.5, 440
    w = crval + cd * np.arange(nw)
    frame = write_icube(synthetic_ccddata, tmp_redux, w - 5000.0, crval, cd)
    out = run(run_primitive, frame)
    got = np.asarray(out.data[:, 0, 0], dtype=np.float64) + 5000.0
    v = velocity(kcwi_config)
    # source: published @ c0fed5c (Ciddor eq. 1, inverted) and independent (v)
    expected = vac_to_air(grid(out.header) / (1.0 + v / C_KMS))
    inner = slice(20, -20)      # keep away from the extrapolated ends
    return got[inner], expected[inner], grid(out.header)[inner]


def test_ramp_maps_pixels_to_corrected_wavelengths(run_primitive,
                                                   synthetic_ccddata,
                                                   tmp_redux, kcwi_config):
    got, expected, wave = ramp_case(run_primitive, synthetic_ccddata,
                                    tmp_redux, kcwi_config)
    # Tolerance: 2e-3 A, chosen to be insensitive to B19 (an error of about
    # 1.7e-4 A on this 0.5 A/px grid, see below) and to float32 storage of
    # values ~100 (8e-6); everything else (formula coefficient rounding,
    # order of the two corrections, spline error) is < 1e-5 A.
    np.testing.assert_allclose(got, expected, rtol=0, atol=2e-3)
    # source: analytic @ c0fed5c - the mapping really moves light: the air
    # wavelength in a pixel is (n - 1) * 5000 = 1.4 A plus v/c * 5000 A (|v| <
    # 31 km/s, i.e. < 0.52 A) below the pixel's grid wavelength.
    assert np.all(wave - got > 0.8) and np.all(wave - got < 2.0)


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=(
    "B19: get_wav_axis puts CRVAL3 on 0-based pixel CRPIX3 (one pixel low), "
    "and air2vac/heliocentric resample on that axis, so light lands about "
    "CD3_3 * ((n - 1) + lambda dn/dlambda + v/c) ~ 1.7e-4 A (0.5 A/px) away "
    "from its corrected wavelength"))
def test_ramp_maps_pixels_exactly(run_primitive, synthetic_ccddata,
                                  tmp_redux, kcwi_config):
    got, expected, _ = ramp_case(run_primitive, synthetic_ccddata,
                                 tmp_redux, kcwi_config)
    # Tolerance: 2e-5 A.  With a correct axis the remaining differences are
    # float32 storage of values ~100 (<= 4e-6 A), the code's 1.67918e-3 vs
    # Ciddor's published 1.67917e-3 (~1e-6 A), applying v before rather than
    # after the air index (~1e-8 A) and the spline error (< 1e-8 A).
    # Measured at c0fed5c: max |got - expected| = 1.70e-4 A (B19); with
    # get_wav_axis patched to the FITS convention, 6.9e-6 A.
    np.testing.assert_allclose(got, expected, rtol=0, atol=2e-5)


def test_no_velocity_correction_only_air_to_vacuum(run_primitive,
                                                   synthetic_ccddata,
                                                   tmp_redux):
    crval, cd, nw = 4900.0, 0.5, 440
    w = crval + cd * np.arange(nw)
    sigma = 1.5
    line = 1000.0 * np.exp(-0.5 * ((w - O3_AIR) / sigma) ** 2)
    frame = write_icube(synthetic_ccddata, tmp_redux, line, crval, cd)
    out = run(run_primitive, frame, radial_velocity_correction="none")
    spec = np.asarray(out.data[:, 0, 0], dtype=np.float64)
    wave = grid(out.header)
    win = np.abs(wave - O3_VAC) < 10 * sigma
    centre = np.sum(wave[win] * spec[win]) / np.sum(spec[win])
    # source: published @ c0fed5c - [O III] vacuum wavelength; tolerance as
    # above
    assert centre == pytest.approx(O3_VAC, abs=0.01)
    # source: analytic @ c0fed5c - docstring
    assert out.header["VCORR"] == 0.0
    assert out.header["VCORRTYP"] == "none"


def test_missing_icube_fails_precondition(run_primitive, synthetic_ccddata,
                                          tmp_redux):
    # source: analytic @ c0fed5c - docstring of _pre_condition - no icube file,
    # nothing written
    frame = synthetic_ccddata()
    run_primitive(WavelengthCorrections, frame)
    assert not (tmp_redux.redux / (NAME + "_icubew.fits")).exists()
