"""Characterization tests for kcwidrp.core.bspline (task T.3).

Smooth analytic functions are fitted with iterfit() and evaluated with
Bspline.value().  The expected values are the functions themselves.
"""
import numpy as np
import pytest
from numpy.linalg import LinAlgError

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


# Banded Cholesky fallback (#231, fixed in #232).  When scipy's
# cholesky_banded raises LinAlgError, cholesky_band redoes the decomposition
# column by column to find the failing column, which fit() then masks.

BW = 4


def to_band(a, n):
    """Dense symmetric matrix -> padded lower-band form used by Bspline."""
    band = np.zeros((BW, n + BW))
    for r in range(BW):
        band[r, 0:n - r] = np.diag(a, -r)
    return band


def random_spd(n, seed=0):
    """Symmetric positive definite matrix of bandwidth BW (L @ L.T)."""
    rng = np.random.default_rng(seed)
    low = np.zeros((n, n))
    for r in range(BW):
        idx = np.arange(n - r)
        low[idx + r, idx] = rng.uniform(1, 2, n - r) if r == 0 else \
            rng.uniform(-0.5, 0.5, n - r)
    return low @ low.T


def force_fallback(monkeypatch):
    def raise_linalg(*args, **kwargs):
        raise LinAlgError("forced")
    monkeypatch.setattr(Bspline, "cholesky_banded", raise_linalg)


def test_cholesky_fallback_matches_scipy(monkeypatch):
    # source: independent - scipy.linalg.cholesky_banded (the normal path)
    # on the same positive definite band matrix
    n = 50
    band = to_band(random_spd(n), n)
    err, expected = Bspline.cholesky_band(band)
    assert err == -1
    force_fallback(monkeypatch)
    err, result = Bspline.cholesky_band(band)
    assert err == -1
    assert result.shape == band.shape
    # Tolerance: two float64 factorizations of a well-conditioned 50 x 50
    # matrix (diagonal 1-5, off-diagonals below 2) differ only in round-off,
    # ~1e-15; 1e-12 still fails if any earlier term is not subtracted.
    np.testing.assert_allclose(result, expected, rtol=1e-12, atol=1e-12)


def test_cholesky_not_positive_definite_returns_index():
    # source: analytic - with a[20, 21] = 5 and diagonal 2, the leading
    # minor of order 22 contains the block [[2, 5], [5, 2]] with determinant
    # 4 - 25 < 0, so the pivot of column 21 is not positive.  The diagonal
    # is positive, so the cheap "bad entries" check passes first.
    n = 50
    a = np.eye(n) * 2.0
    a[21, 20] = a[20, 21] = 5.0
    err, out = Bspline.cholesky_band(to_band(a, n))
    assert isinstance(err, np.ndarray)
    assert err.tolist() == [21]


def test_maskpoints_accepts_an_int():
    # source: analytic - documented: err may be an int; -1 means the error
    # points were masked, which clears at least one breakpoint
    x = np.linspace(0, 1, 200)
    bkpt = np.linspace(0, 1, 21)
    sset = Bspline.Bspline(x, fullbkpt=np.concatenate(
        [[0] * 3, bkpt, [1] * 3]))
    assert sset.maskpoints(10) == -1
    assert not sset.mask.all()


def test_iterfit_survives_cholesky_fallback(monkeypatch):
    # source: independent - the same fit with scipy's factorization
    rng = np.random.default_rng(1)
    x = np.sort(rng.uniform(0, 10, 2000))
    y = np.sin(x) + rng.normal(0, 0.01, x.size)
    bkpt = np.linspace(x.min(), x.max(), 41)
    expected, _ = Bspline.iterfit(x, y, fullbkpt=bkpt)
    force_fallback(monkeypatch)
    sset, _ = Bspline.iterfit(x, y, fullbkpt=bkpt)
    # Tolerance: the two factorizations differ in round-off only; the fit
    # values (order 1) agree to ~1e-13, and 1e-8 still fails if the fallback
    # masks a breakpoint the scipy path keeps.
    np.testing.assert_allclose(sset.value(x)[0], expected.value(x)[0],
                               rtol=1e-8, atol=1e-10)
