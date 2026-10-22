"""Characterization tests for WavelengthCorrections (task T.3):
a2v_conversion, the heliocentric/barycentric velocity correction, and
get_wav_axis.

These are methods of a framework primitive, but they only use ``self.logger``,
``self.config.instrument`` (latitude, longitude, altitude), ``self.get_wav_axis``
and ``self.action.args.cwave``.  The tests call the functions on the class
with a ``types.SimpleNamespace`` standing in for ``self`` -- no framework,
queue or context is created.  Not isolated: nothing here needed the whole
primitive; ``air2vac`` (the cube loop around a2v_conversion) and ``_perform``
are not covered by T.3.
"""
import logging
from types import SimpleNamespace

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits

from kcwidrp.primitives.WavelengthCorrections import WavelengthCorrections as WC

LOG = logging.getLogger("test_wavelength_corrections")
C_KMS = 299792.458          # speed of light, km/s (exact, SI)


def _stub(config=None, cwave=4000.0):
    stub = SimpleNamespace(logger=LOG, config=config,
                           action=SimpleNamespace(
                               args=SimpleNamespace(cwave=cwave)))
    stub.get_wav_axis = lambda header: WC.get_wav_axis(stub, header)
    return stub


# ------------------------------------------------------------ a2v_conversion
# (air Angstrom, vacuum Angstrom)
# source: published - standard-air wavelengths of strong lines and their
# vacuum values as listed in the SDSS rest-wavelength line list and the NIST
# Atomic Spectra Database, both of which use the IAU conversion of Morton
# (1991, ApJS 77, 119) based on Edlen (1966):
#   [O II]  3726.032 -> 3727.092 ; 3728.815 -> 3729.875
#   [O III] 5006.843 -> 5008.239
#   Na I D2 5889.951 -> 5891.583
#   H-alpha 6562.80  -> 6564.614
LINES = [(3726.032, 3727.092), (3728.815, 3729.875), (5006.843, 5008.239),
         (5889.951, 5891.583), (6562.80, 6564.614)]


def test_known_lines_to_one_hundredth_angstrom():
    air = np.array([a for a, _ in LINES]) * u.AA
    vac = WC.a2v_conversion(_stub(), air)
    # Tolerance: 0.01 Angstrom (the task's stated accuracy).  The tabulated
    # values are rounded to 0.001 and Morton's formula (Edlen) differs from
    # the Ciddor (1996) formula used here by < 0.002 A in this range.
    np.testing.assert_allclose(vac.to_value(u.AA), [v for _, v in LINES],
                               rtol=0, atol=0.01)


def test_output_is_quantity_in_angstrom_and_longer_than_input():
    # source: analytic - n(air) > 1, so lambda_vac = n * lambda_air > lambda_air
    air = np.array([4000.0, 5000.0, 7000.0]) * u.AA
    vac = WC.a2v_conversion(_stub(), air)
    assert vac.unit == u.AA
    assert np.all(vac.value > air.value)


def test_nanometre_input_is_converted():
    # source: published - [O III] 500.6843 nm air -> 5008.239 A vacuum (same
    # tolerance as above)
    vac = WC.a2v_conversion(_stub(), np.array([500.6843, 500.6843]) * u.nm)
    np.testing.assert_allclose(vac.to_value(u.AA), 5008.239, rtol=0,
                               atol=0.01)


def test_below_2000_angstrom_is_unchanged():
    # source: characterization @ b0e2550 - the conversion is only applied from
    # 2000 A upward; 1999 A comes back identical, 2000 A is shifted
    air = np.array([1500.0, 1999.0, 2000.0, 2500.0]) * u.AA
    vac = WC.a2v_conversion(_stub(), air).value
    assert vac[0] == 1500.0 and vac[1] == 1999.0
    assert vac[2] > 2000.0 and vac[3] > 2500.0


def test_shift_grows_monotonically_with_wavelength():
    # source: analytic - (lambda_vac - lambda_air) = (n - 1) lambda, which
    # increases with lambda over the optical range since n - 1 changes by
    # < 10 % over 3500-9000 A while lambda grows 2.6x
    air = np.linspace(3500.0, 9000.0, 50) * u.AA
    d = WC.a2v_conversion(_stub(), air).value - air.value
    assert np.all(np.diff(d) > 0)


# ------------------------------------------- heliocentric / barycentric
KECK = dict(lat=19.82656, lon=-155.4742, height=4160)   # kcwi.cfg values
ARGS = dict(TARGRA=150.1, TARGDEC=2.2)


def _cfg(kcwi_config):
    return SimpleNamespace(instrument=kcwi_config())


def _obj(date, ra=150.1, dec=2.2, crval3=4000.0, cd33=1.0, **extra):
    hdr = fits.Header({"TARGRA": ra, "TARGDEC": dec, "DATE-BEG": date,
                       "CRVAL3": crval3, "CD3_3": cd33})
    for k, v in extra.items():
        hdr[k] = v
    return SimpleNamespace(data=np.zeros((3, 2, 2)), header=hdr)


@pytest.mark.parametrize("date", ["2026-01-15T10:00:00.000",
                                  "2023-05-01T08:30:00.000"])
@pytest.mark.parametrize("mode", ["heliocentric", "barycentric"])
def test_velocity_matches_astropy(kcwi_config, date, mode):
    # source: independent - astropy SkyCoord.radial_velocity_correction for the
    # same coordinates, time and Keck location (19.82656 N, 155.4742 W, 4160 m).
    obj = _obj(date)
    WC.heliocentric(_stub(_cfg(kcwi_config)), obj, mode, resample=False)
    target = SkyCoord(150.1, 2.2, unit="deg", obstime=date)
    keck = EarthLocation.from_geodetic(**KECK)
    expected = target.radial_velocity_correction(
        kind=mode, location=keck).to_value("km/s")
    # Tolerance: 1 m/s = 1e-3 km/s (task requirement); same astropy call.
    assert obj.header["VCORR"] == pytest.approx(expected, rel=0, abs=1e-3)
    assert obj.header["VCORRTYP"] == mode


def test_velocity_is_bounded_by_earth_orbital_speed(kcwi_config):
    # source: analytic - the helio/barycentric correction cannot exceed the
    # Earth's orbital speed (~29.8 km/s, up to 30.3 at perihelion) plus the
    # site rotation speed (0.46 km/s * cos(lat) = 0.44 km/s): < 31 km/s
    for date in ["2026-01-15T10:00:00.000", "2025-07-15T10:00:00.000"]:
        for mode in ("heliocentric", "barycentric"):
            obj = _obj(date, ra=0.0, dec=0.0)
            WC.heliocentric(_stub(_cfg(kcwi_config)), obj, mode,
                            resample=False)
            assert abs(obj.header["VCORR"]) < 31.0


def test_velocity_is_small_toward_the_ecliptic_pole(kcwi_config):
    # source: analytic - toward the north ecliptic pole (RA 18h, Dec +66.56
    # deg) the Earth's orbital velocity is perpendicular to the line of sight,
    # leaving the rotation term (<= 0.44 km/s) and the small eccentricity /
    # lunar / barycentric offsets: |v| < 1 km/s.
    obj = _obj("2026-04-10T10:00:00.000", ra=270.0, dec=66.56)
    WC.heliocentric(_stub(_cfg(kcwi_config)), obj, "heliocentric",
                    resample=False)
    assert abs(obj.header["VCORR"]) < 1.0


def test_wavelength_zero_point_scaled_without_resampling():
    # source: analytic - with an explicit vcorr = 30 km/s and resample=False,
    # CRVAL3 and CD3_3 are multiplied by (1 + v/c).
    obj = _obj("2026-01-15T10:00:00.000", crval3=4000.0, cd33=1.5)
    WC.heliocentric(_stub(), obj, "heliocentric", resample=False, vcorr=30.0)
    f = 1.0 + 30.0 / C_KMS
    # Tolerance: one multiplication and one division in float64, 2 ulp.
    assert obj.header["CRVAL3"] == pytest.approx(4000.0 * f, rel=4e-16, abs=0)
    assert obj.header["CD3_3"] == pytest.approx(1.5 * f, rel=4e-16, abs=0)
    assert obj.header["VCORR"] == 30.0


def test_existing_vcorr_is_rolled_back():
    # source: analytic - a header already corrected by 5 km/s, corrected to
    # 12 km/s, is scaled by (1 + (12 - 5)/c)
    obj = _obj("2026-01-15T10:00:00.000", crval3=4000.0, VCORR=5.0)
    WC.heliocentric(_stub(), obj, "heliocentric", resample=False, vcorr=12.0)
    # Tolerance: as above, 2 ulp.
    assert obj.header["CRVAL3"] == pytest.approx(
        4000.0 * (1.0 + 7.0 / C_KMS), rel=4e-16, abs=0)
    assert obj.header["VCORR"] == 12.0


def _cube_header(nw=20, crval=4000.0, cd=1.0):
    return fits.Header({"NAXIS": 3, "NAXIS1": 2, "NAXIS2": 2, "NAXIS3": nw,
                        "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN",
                        "CTYPE3": "AWAV", "CRVAL3": crval, "CRPIX3": 1.0,
                        "CD3_3": cd})


def test_resampled_constant_cube_is_unchanged():
    # source: analytic - a spectrum that is constant in wavelength is still
    # constant after being re-sampled on the corrected grid (cubic
    # interpolation reproduces constants); the header records VCORR.
    hdr = _cube_header()
    obj = SimpleNamespace(data=np.full((20, 2, 2), 7.0), header=hdr)
    WC.heliocentric(_stub(), obj, "heliocentric", resample=True, vcorr=10.0)
    # Tolerance: cubic spline of a constant in float64, error ~1e-14.
    np.testing.assert_allclose(obj.data, 7.0, rtol=1e-12, atol=0)
    assert obj.header["VCORR"] == 10.0


def test_resampled_mask_flags_pixels_shifted_off_the_grid():
    # source: analytic - mask cubes are resampled with fill value 128 outside
    # the shifted grid.  At +50 km/s the grid shifts by 4000 * 50/c = 0.67 A,
    # less than one 1 A pixel, so only the first pixel (below the shifted
    # grid's start) is flagged; every other pixel keeps its 0.
    hdr = _cube_header(nw=20, cd=1.0)
    obj = SimpleNamespace(data=np.zeros((20, 2, 2)), header=hdr)
    WC.heliocentric(_stub(), obj, "heliocentric", mask=True, resample=True,
                    vcorr=50.0)
    assert np.all(obj.data[0] == 128)
    assert np.all(obj.data[1:] == 0)


# ------------------------------------------------------------- get_wav_axis
def test_wav_axis_length_and_step():
    # source: analytic - 5 pixels with CD3_3 = 1.5
    axis = WC.get_wav_axis(_stub(), _hdr5())
    assert len(axis) == 5
    np.testing.assert_allclose(np.diff(axis), 1.5, rtol=0, atol=1e-12)


def _hdr5():
    hdr = _cube_header(nw=5, crval=4000.0, cd=1.5)
    return hdr


@pytest.mark.xfail(strict=True, reason="B19: get_wav_axis ignores the "
                   "FITS 1-based pixel convention: wav0 + (i - CRPIX)*dw "
                   "puts CRVAL3 on 0-based pixel CRPIX3, not CRPIX3 - 1 "
                   "(WavelengthCorrections.py get_wav_axis; the DRP writes "
                   "CRPIX3 = 1 with CRVAL3 the first pixel)")
def test_wav_axis_first_pixel_is_crval_when_crpix_is_one():
    # source: published - FITS WCS paper I: for 1-based pixel p the world
    # coordinate is CRVAL + (p - CRPIX)*CDELT; with CRPIX3 = 1 the first
    # pixel (p = 1, 0-based index 0) has wavelength CRVAL3 = 4000.
    axis = WC.get_wav_axis(_stub(), _hdr5())
    np.testing.assert_allclose(axis, 4000.0 + 1.5 * np.arange(5), rtol=0,
                               atol=1e-9)
