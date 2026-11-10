"""Properties of the differential atmospheric refraction (task T.5).

CorrectDar.atm_disper(w0, w1, airmass, temperature, pressure_pa, humidity,
co2) returns the displacement (arcsec) of an image at wavelength w1 relative
to one at w0.  For a plane-parallel atmosphere the refraction is
R(lambda) = (n(lambda) - 1) tan z (e.g. Filippenko 1982, PASP 94, 715,
eq. 1), so the differential refraction is (n(w0) - n(w1)) tan z with
sec z = airmass.  Hence, for any physical atmosphere:

1. the displacement of a wavelength relative to itself is zero;
2. swapping the two wavelengths reverses the displacement;
3. its size grows with airmass (tan z grows with sec z for z in [0, 90));
4. its size grows with the wavelength separation, on either side of the
   reference: n(lambda) is strictly decreasing in lambda over 3000-11000 A
   (normal dispersion of air and water vapour; no absorption bands are
   modelled), so |n(w0) - n(w1)| grows with |w1 - w0|.

Atmospheres are drawn over the range the DRP sees and the Ciddor (1996)
model is valid for: -10..25 C, 55-105 kPa (Maunakea summit is ~61 kPa),
0-100 % relative humidity, 300-500 ppm CO2.  T.3 pins single values of 1-3
with fixed numbers; here they hold for random wavelengths and atmospheres.
"""
from hypothesis import given
from hypothesis import strategies as st

from kcwidrp.primitives.CorrectDar import atm_disper

wave = st.floats(min_value=3000.0, max_value=11000.0)
airmass = st.floats(min_value=1.0, max_value=3.0)
atmosphere = st.fixed_dictionaries({
    "temperature": st.floats(min_value=-10.0, max_value=25.0),
    "pressure_pa": st.floats(min_value=55000.0, max_value=105000.0),
    "humidity": st.floats(min_value=0.0, max_value=100.0),
    "co2": st.floats(min_value=300.0, max_value=500.0),
})


@given(w=wave, x=airmass, atm=atmosphere)
def test_no_displacement_at_the_reference_wavelength(w, x, atm):
    # source: analytic - n(w) - n(w) = 0
    assert atm_disper(w, w, x, **atm) == 0.0


@given(w0=wave, w1=wave, x=airmass, atm=atmosphere)
def test_swapping_wavelengths_reverses_the_displacement(w0, w1, x, atm):
    # source: analytic - (n0 - n1) tan z = -(n1 - n0) tan z.  Exact: IEEE
    # subtraction is antisymmetric and negation is exact.
    assert atm_disper(w1, w0, x, **atm) == -atm_disper(w0, w1, x, **atm)


@given(w0=wave, w1=wave, x1=airmass, dx=st.floats(min_value=0.01,
                                                  max_value=2.0),
       atm=atmosphere)
def test_displacement_grows_with_airmass(w0, w1, x1, dx, atm):
    # source: analytic - tan(arcsec X) = sqrt(X^2 - 1) is strictly increasing
    # for X >= 1.  Airmasses at least 0.01 apart and wavelengths at least
    # 10 A apart, so the two values differ far above rounding.
    if abs(w1 - w0) < 10.0:
        w1 = w0 + 10.0 if w0 < 10990.0 else w0 - 10.0
    lo = abs(atm_disper(w0, w1, x1, **atm))
    hi = abs(atm_disper(w0, w1, x1 + dx, **atm))
    assert hi > lo


@st.composite
def three_wavelengths(draw):
    """a < b < c in 3000-11000 A, neighbours at least 1 A apart.

    1 A: dn/dlambda ~ 1e-9 per A at 11000 A, ~1e7 times the float64
    resolution of n, so the comparisons below are far above rounding.
    """
    a = draw(st.floats(min_value=3000.0, max_value=10998.0))
    b = draw(st.floats(min_value=a + 1.0, max_value=10999.0))
    c = draw(st.floats(min_value=b + 1.0, max_value=11000.0))
    return a, b, c


@given(abc=three_wavelengths(), redward=st.booleans(),
       x=st.floats(min_value=1.01, max_value=3.0), atm=atmosphere)
def test_displacement_grows_with_wavelength_separation(abc, redward, x, atm):
    # source: analytic - n(lambda) strictly decreasing (module docstring).
    # Reference at one end; the middle wavelength is nearer to it than the
    # far end, on the red side (reference a) or the blue side (reference c).
    a, b, c = abc
    ref, near, far = (a, b, c) if redward else (c, b, a)
    assert abs(atm_disper(ref, far, x, **atm)) > \
        abs(atm_disper(ref, near, x, **atm))
