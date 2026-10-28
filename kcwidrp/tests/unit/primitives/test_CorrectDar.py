"""Characterization tests for CorrectDar.atm_disper (task T.3).

atm_disper(w0, w1, airmass, ...) returns the atmospheric dispersion in
arcsec of wavelength w1 relative to w0: 206265 * (n(w0) - n(w1)) * tan(z) with
sec z = airmass.

Expected values come from the published dispersion formula for standard air
(15 C, 101325 Pa, dry, 450 ppm CO2) -- Ciddor (1996, Appl. Opt. 35, 1566),
eq. 1, which is the Peck & Reeder (1972, JOSA 62, 958) / Edlen form:

    n - 1 = 1e-8 * ( 5792105/(238.0185 - s^2) + 167917/(57.362 - s^2) ),

s = 1/lambda in um^-1 -- evaluated in the test with plain arithmetic, whereas
the code under test calls the ``ref_index.ciddor`` package (full Ciddor
algorithm with temperature, pressure, humidity and CO2).  The two agree for
standard air.  The geometry (tan z for sec z = airmass) is exact trigonometry.
"""
import math

import pytest

from kcwidrp.primitives.CorrectDar import atm_disper

ARCSEC_PER_RAD = 206265.0   # the constant the DRP uses (206264.806 exact)


def _n_standard_air(lam_angstrom):
    # source: published - Ciddor (1996) eq. 1 / Peck & Reeder (1972)
    s2 = (1.0e4 / lam_angstrom) ** 2
    return 1.0 + 1.0e-8 * (5792105.0 / (238.0185 - s2)
                           + 167917.0 / (57.362 - s2))


STD = dict(temperature=15.0, pressure_pa=101325.0, humidity=0.0, co2=450.0)


@pytest.mark.parametrize("w0, w1, airmass", [
    (5000.0, 4000.0, 1.5),
    (5000.0, 7000.0, 1.5),
    (4500.0, 3500.0, 2.0),
])
def test_standard_air_matches_published_formula(w0, w1, airmass):
    z = math.acos(1.0 / airmass)
    expected = ARCSEC_PER_RAD * (_n_standard_air(w0) - _n_standard_air(w1)) \
        * math.tan(z)
    got = atm_disper(w0, w1, airmass, **STD)
    # Tolerance: Ciddor's full algorithm at standard air and eq. 1 differ only
    # by rounding of the published coefficients; measured relative difference
    # 6e-11, so 1e-9 leaves a factor ~15.
    assert got == pytest.approx(expected, rel=1e-9, abs=0)


def test_bluer_wavelength_is_displaced_negatively():
    # source: analytic - refractive index rises toward the blue, so
    # (n0 - n1) < 0 for w1 < w0 and the offset is negative; positive for w1 > w0
    assert atm_disper(5000.0, 4000.0, 1.5) < 0.0
    assert atm_disper(5000.0, 7000.0, 1.5) > 0.0


def test_zero_offset_at_reference_wavelength():
    # source: analytic - w1 == w0 gives n0 - n1 = 0 exactly
    assert atm_disper(5000.0, 5000.0, 1.5) == 0.0


def test_zero_at_zenith():
    # source: analytic - airmass 1 is z = 0, tan z = 0
    assert atm_disper(4000.0, 5000.0, 1.0) == 0.0


def test_antisymmetric_in_wavelength_order():
    # source: analytic - swapping w0 and w1 changes the sign exactly
    a = atm_disper(4000.0, 6000.0, 1.3)
    b = atm_disper(6000.0, 4000.0, 1.3)
    assert a == -b


def test_scales_with_tan_zenith_distance():
    # source: analytic - airmass 2.0 -> z = 60 deg, tan z = sqrt(3);
    # airmass 1.25 -> cos z = 0.8, sin z = 0.6, tan z = 0.75.
    hi = atm_disper(4000.0, 6000.0, 2.0)
    lo = atm_disper(4000.0, 6000.0, 1.25)
    # Tolerance: acos/tan in float64, each correct to ~1 ulp; 1e-12 relative.
    assert hi / lo == pytest.approx(math.sqrt(3.0) / 0.75, rel=1e-12, abs=0)


def test_mauna_kea_conditions_scale_with_air_density():
    # source: analytic - (n - 1) is proportional to the air density, i.e. to
    # P/T for an ideal gas, so the dispersion at 10 C and 61100 Pa (dry) is
    # the standard-air value times (61100/101325)*(288.15/283.15).
    z = math.acos(1.0 / 1.5)
    std = ARCSEC_PER_RAD * (_n_standard_air(5000.0) - _n_standard_air(4000.0)) \
        * math.tan(z)
    expected = std * (61100.0 / 101325.0) * (288.15 / 283.15)
    got = atm_disper(5000.0, 4000.0, 1.5, temperature=10.0, pressure_pa=61100.0,
                     humidity=0.0, co2=450.0)
    # Tolerance: air is not exactly an ideal gas (compressibility departs from
    # 1 by ~1e-3 at 1 atm, less at 0.6 atm); measured difference 1.3e-4, so
    # 5e-4 relative.
    assert got == pytest.approx(expected, rel=5e-4, abs=0)
