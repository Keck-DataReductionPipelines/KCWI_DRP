import numpy as np
from numpy.linalg import LinAlgError

from kcwidrp.core.bspline import Bspline

BW = 4


def to_band(a, n):
    """Dense symmetric matrix -> padded lower-band form used by Bspline."""
    band = np.zeros((BW, n + BW))
    for r in range(BW):
        band[r, 0:n - r] = np.diag(a, -r)
    return band


def random_spd(n, seed=0):
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


def test_fallback_matches_scipy(monkeypatch):
    n = 50
    band = to_band(random_spd(n), n)
    err, expected = Bspline.cholesky_band(band)
    assert err == -1
    force_fallback(monkeypatch)
    err, result = Bspline.cholesky_band(band)
    assert err == -1
    assert result.shape == band.shape
    np.testing.assert_allclose(result, expected, rtol=1e-12, atol=1e-12)


def test_not_positive_definite_returns_index():
    # Positive diagonal, so the cheap "bad entries" check passes, but the
    # oversized off-diagonal makes the matrix indefinite.
    n = 50
    a = np.eye(n) * 2.0
    a[21, 20] = a[20, 21] = 5.0
    err, out = Bspline.cholesky_band(to_band(a, n))
    assert isinstance(err, np.ndarray)
    assert err.tolist() == [21]


def test_maskpoints_accepts_int():
    x = np.linspace(0, 1, 200)
    bkpt = np.linspace(0, 1, 21)
    sset = Bspline.Bspline(x, fullbkpt=np.concatenate(
        [[0] * 3, bkpt, [1] * 3]))
    assert sset.maskpoints(10) == -1
    assert not sset.mask.all()


def test_iterfit_survives_fallback(monkeypatch):
    rng = np.random.default_rng(1)
    x = np.sort(rng.uniform(0, 10, 2000))
    y = np.sin(x) + rng.normal(0, 0.01, x.size)
    bkpt = np.linspace(x.min(), x.max(), 41)
    expected, _ = Bspline.iterfit(x, y, fullbkpt=bkpt)
    force_fallback(monkeypatch)
    sset, _ = Bspline.iterfit(x, y, fullbkpt=bkpt)
    np.testing.assert_allclose(sset.value(x)[0], expected.value(x)[0],
                               rtol=1e-8, atol=1e-10)
