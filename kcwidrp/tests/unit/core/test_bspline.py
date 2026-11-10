"""Characterization tests for kcwidrp.core.bspline (task T.3).

Smooth analytic functions are fitted with iterfit() and evaluated with
Bspline.value().  The expected values are the functions themselves.
"""
import numpy as np
import pytest

from kcwidrp.core.bspline import Bspline


def test_cubic_polynomial_is_reproduced_exactly():
    # source: analytic - a cubic B-spline (nord=4) space contains every cubic
    # polynomial, so a least-squares fit of an exact cubic recovers it.
    x = np.linspace(0.0, 10.0, 400)
    y = 2.0 + 0.5 * x - 0.03 * x**2 + 0.001 * x**3
    sset, outmask = Bspline.iterfit(x, y, fullbkpt=np.linspace(0, 10, 11))
    yfit, _ = sset.value(x)
    # Tolerance: normal-equation solve in float64 on 400 points, measured max
    # error 2e-13 on values of order 1-7; 1e-10 leaves a margin of 500.
    np.testing.assert_allclose(yfit, y, rtol=0, atol=1e-10)
    assert np.all(outmask)


def test_sine_fit_within_spline_error_bound():
    # source: analytic - sin(x) on [0, 2 pi] with 40 equal intervals.
    x = np.linspace(0.0, 2 * np.pi, 500)
    y = np.sin(x)
    sset, outmask = Bspline.iterfit(x, y, fullbkpt=np.linspace(0, 2*np.pi, 41))
    yfit, _ = sset.value(x)
    # Tolerance: cubic-spline error ~ c*h^4*max|f''''| with h = 0.157,
    # |sin''''| <= 1 -> a few 1e-5 (measured 5e-5 at the ends); 1e-4 still
    # fails for a wrong spline order or knot placement.
    np.testing.assert_allclose(yfit, y, rtol=0, atol=1e-4)


def test_evaluation_off_the_data_points():
    # source: analytic - the fit evaluated at x not used in the fit follows
    # the function, same bound as above
    x = np.linspace(0.0, 2 * np.pi, 500)
    sset, _ = Bspline.iterfit(x, np.sin(x),
                              fullbkpt=np.linspace(0, 2*np.pi, 41))
    xe = np.array([0.1234, 1.0, 3.3, 6.0])
    ye, _ = sset.value(xe)
    # Tolerance: same bound as test_sine_fit_within_spline_error_bound.
    np.testing.assert_allclose(ye, np.sin(xe), rtol=0, atol=1e-4)


def test_outlier_is_rejected_and_fit_unaffected():
    # source: analytic - one point 5 above sin(x) with invvar for sigma=0.01
    # is rejected (outmask False) and the fit still follows sin(x).
    x = np.linspace(0.0, 2 * np.pi, 500)
    y = np.sin(x)
    y[250] += 5.0
    sset, outmask = Bspline.iterfit(
        x, y, invvar=np.full(500, 1.0 / 0.01**2), upper=3, lower=3,
        fullbkpt=np.linspace(0, 2*np.pi, 41))
    assert not outmask[250]
    yfit, _ = sset.value(x)
    clean = np.arange(500) != 250
    # Tolerance: same spline-error bound as the sine test.
    np.testing.assert_allclose(yfit[clean], np.sin(x)[clean], rtol=0,
                               atol=1e-4)


def test_mismatched_sizes_rejected():
    # source: analytic - documented: dimensions of xdata and ydata must agree
    with pytest.raises(ValueError):
        Bspline.iterfit(np.arange(10.0), np.arange(9.0))
