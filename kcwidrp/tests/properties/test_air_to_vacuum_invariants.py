"""Properties of the air-to-vacuum wavelength conversion (task T.5).

WavelengthCorrections.a2v_conversion(wave) converts standard-air wavelengths
to vacuum: lambda_vac = n(lambda) * lambda_air.  Over the KCWI range
(3000-11000 A) dry standard air has normal dispersion and no absorption, so:

1. the conversion is strictly increasing (it is a change of variable on the
   same spectrum: a redder line stays redder);
2. vacuum wavelengths are longer than air wavelengths (n > 1);
3. (vac - air) / air = n - 1 lies in the physical range for standard air.

Range in 3: n - 1 from Edlen (1966, Metrologia 2, 71), eq. for standard
air (15 C, 101325 Pa, dry, 300 ppm CO2):

    (n - 1) * 1e8 = 8342.13 + 2406030 / (130 - s^2) + 15997 / (38.9 - s^2)

with s = 1 / lambda in um^-1, evaluated at the drawn wavelength: 2.9155e-4
at 3000 A, 2.7389e-4 at 11000 A.  The conversion must agree with it to 0.1 %
of (n - 1), see the tolerance comment.  This is a different published
formula from the one in the code (Ciddor 1996), so it is not a restatement.

The method only uses ``self.logger``; it is called on the class with a
SimpleNamespace standing in for ``self``, as in the T.3 tests.
T.3 already pins 5 known lines and a 3-point monotonicity check; here the
properties are checked on random wavelength arrays.
"""
import logging
from types import SimpleNamespace

import numpy as np
from astropy import units as u
from hypothesis import given
from hypothesis import strategies as st

from kcwidrp.primitives.WavelengthCorrections import WavelengthCorrections as WC

STUB = SimpleNamespace(logger=logging.getLogger("test_a2v_properties"))
LO, HI = 3000.0, 11000.0    # KCWI blue + red wavelength range, Angstrom

wavelength = st.floats(min_value=LO, max_value=HI, allow_nan=False)


def _edlen_n_minus_1(lam_angstrom):
    # source: published - Edlen (1966), standard air (module docstring)
    s2 = (1.0e4 / lam_angstrom) ** 2
    return 1.0e-8 * (8342.13 + 2406030.0 / (130.0 - s2)
                     + 15997.0 / (38.9 - s2))


@st.composite
def increasing_wavelengths(draw):
    """1-30 strictly increasing wavelengths in [LO, HI].

    Neighbours are at least 1e-6 A apart: d(lambda_vac)/d(lambda_air) is
    about 1.0003, and the rounding error of the conversion is a few ulp of
    lambda (ulp(1e4) = 1.8e-12 A), so steps of 1e-6 A are resolved with a
    margin of ~1e5 while still probing very close pairs.
    """
    start = draw(wavelength)
    steps = draw(st.lists(st.floats(min_value=1e-6, max_value=800.0),
                          min_size=0, max_size=29))
    waves = [start]
    for s in steps:
        if waves[-1] + s > HI:
            break
        waves.append(waves[-1] + s)
    return np.array(waves)


@given(air=increasing_wavelengths())
def test_conversion_is_strictly_increasing(air):
    # source: analytic - a change of variable with d(vac)/d(air) = n + air *
    # dn/d(air) > 0 (normal dispersion, |air dn/dair| ~ 1e-5 << n)
    vac = WC.a2v_conversion(STUB, air * u.AA).to_value(u.AA)
    assert np.all(np.diff(vac) > 0)


@given(air=increasing_wavelengths())
def test_vacuum_is_longer_than_air(air):
    # source: analytic - the refractive index of air is greater than 1
    vac = WC.a2v_conversion(STUB, air * u.AA).to_value(u.AA)
    assert np.all(vac > air)


@given(air=increasing_wavelengths())
def test_shift_is_within_standard_air_refractivity(air):
    vac = WC.a2v_conversion(STUB, air * u.AA).to_value(u.AA)
    got = vac / air - 1.0
    expected = _edlen_n_minus_1(air)   # source: published - Edlen (1966)
    # Tolerance: 0.1 % of (n - 1), i.e. 3e-7 in n.  Standard-air formulas
    # differ by more than rounding: Edlen (1966) assumes 300 ppm CO2, Ciddor
    # (1996) 450 ppm (+0.008 % in n - 1), and evaluating the formula at the
    # air rather than the vacuum wavelength shifts n - 1 by up to 0.01 % at
    # 3000 A.  0.1 % is ten times that spread, and 10x smaller than any change
    # of a coefficient that matters (1 % of n - 1 = 0.03 A at 11000 A).
    np.testing.assert_allclose(got, expected, rtol=1e-3, atol=0)
