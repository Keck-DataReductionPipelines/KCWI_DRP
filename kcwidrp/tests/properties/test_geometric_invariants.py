"""Property: the estimated inverse transform undoes the forward one (T.5).

AsymmetricPolynomialTransform has no explicit inverse (``inverse`` raises and
says to estimate one by exchanging source and destination).  SolveGeom does
exactly that, with orders (2, 4):

    tform = estimate_transform('asympolynomial', src, dst, order=(2, 4))
    invtf = estimate_transform('asympolynomial', dst, src, order=(2, 4))

Invariant: for a well-conditioned map f and control points covering the
detector domain, invtf(tform(x)) = x for every x in the domain, not only at
the control points.

"Well-conditioned" coefficients are generated as
    f(p) = A p + t + q(p),
* A = R(theta) diag(s1, s2) R(phi): any rotation, any shear and scales
  s1, s2 in [0.5, 2], so cond(A) <= 4 and |A^-1| <= 2 (no fold, no
  collapse, pixel scale changes of at most 2x as between binnings);
* t puts the image in the positive quadrant, 1-500 px from the axes, like
  the detector and output-pixel frames SolveGeom uses (the estimator scales
  by the median |x| and |y|, so a frame centred on 0 is outside its domain);
* q is a smooth quadratic distortion with each component bounded by
  D <= 0.5 px over the domain.  With q = 0 the inverse is affine and exactly
  representable; with q != 0 it is not polynomial, and the second-order
  residual |A^-1 Dq A^-1 q| <= 2 * (8 D / L) * 2 * sqrt(2) * D ~ 45 D^2 / L
  = 0.003 px for D = 0.5, L = 4096 (an upper bound: the fit absorbs most of
  it).

Domain: x, y in [1, 4096], the unbinned KCWI detector.  Control points: a
15 x 15 grid over the domain (SolveGeom fits a few hundred bar/arc points).
T.3 checks one affine round trip on a 7 x 7 grid near the origin.
"""
import numpy as np
from hypothesis import given
from hypothesis import strategies as st

from kcwidrp.core.geometric import estimate_transform

L = 4096.0
_g = np.linspace(1.0, L, 15)
_gx, _gy = np.meshgrid(_g, _g)
GRID = np.column_stack([_gx.ravel(), _gy.ravel()])

# Tolerance for the round trip, in pixels: 0.01 px.
# - Rounding: the estimator solves the normal equations with an explicit
#   inverse, so its error is ~cond(W)^2 * eps * L; cond(W) = 9e3 for the
#   15 x 15 grid scaled by its median, giving ~1e-4 px (worst seen over 600
#   random maps: 5e-4 px, on transformed grids whose cond is larger).
# - Truncation: <= 0.003 px for the distortion bound above.
# 0.01 px is 3x the sum and 1/10 of the ~0.1 px accuracy a geometric
# solution needs to place slices; a wrong inverse is off by >= 1 px.
TOL_PX = 0.01


def _rot(a):
    return np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])


@st.composite
def maps(draw, distortion=True):
    """A well-conditioned map f (module docstring)."""
    angle = st.floats(min_value=0.0, max_value=2 * np.pi)
    scale = st.floats(min_value=0.5, max_value=2.0)
    A = _rot(draw(angle)) @ np.diag([draw(scale), draw(scale)]) \
        @ _rot(draw(angle))
    margin = np.array([draw(st.floats(1.0, 500.0)),
                       draw(st.floats(1.0, 500.0))])
    D = draw(st.floats(0.0, 0.5)) if distortion else 0.0
    c = np.array([draw(st.floats(-1.0, 1.0)) for _ in range(6)])
    cx, cy = c[:3], c[3:]
    nx = max(np.abs(cx).sum(), 1e-12)
    ny = max(np.abs(cy).sum(), 1e-12)

    def raw(p):
        u = (p - L / 2) / (L / 2)          # u in [-1, 1] over the domain
        quad = np.column_stack([u[:, 0] ** 2, u[:, 0] * u[:, 1],
                                u[:, 1] ** 2])
        q = D * np.column_stack([quad @ cx / nx, quad @ cy / ny])
        return p @ A.T + q

    # Translate so the image of the domain starts `margin` px from the axes.
    corners = np.array([[1, 1], [1, L], [L, 1], [L, L]], dtype=float)
    lo = np.vstack([raw(GRID), raw(corners)]).min(axis=0)
    t = margin - lo

    return _Map(lambda p: raw(p) + t,
                "A=%s, t=%s, D=%.3g px" % (np.round(A, 4).tolist(),
                                          np.round(t, 2).tolist(), D))


class _Map:
    """A callable map with a readable repr (for hypothesis reports)."""

    def __init__(self, func, text):
        self._func, self._text = func, text

    def __call__(self, p):
        return self._func(p)

    def __repr__(self):
        return "Map(%s)" % self._text


points = st.lists(st.tuples(st.floats(1.0, L), st.floats(1.0, L)),
                  min_size=1, max_size=20).map(np.array)


def _round_trip_error(f, pts):
    dst = f(GRID)
    tform = estimate_transform("asympolynomial", GRID, dst, order=(2, 4))
    invtf = estimate_transform("asympolynomial", dst, GRID, order=(2, 4))
    return np.abs(invtf(tform(pts)) - pts).max()


@given(f=maps(distortion=False), pts=points)
def test_affine_round_trip_over_detector(f, pts):
    # source: analytic - an affine map's inverse is affine, so it is in the
    # fitted polynomial family and the round trip is the identity.
    assert _round_trip_error(f, pts) <= TOL_PX


@given(f=maps(distortion=True), pts=points)
def test_distorted_round_trip_over_detector(f, pts):
    # source: analytic - inverse of a near-affine map, truncation bound in
    # the module docstring.
    assert _round_trip_error(f, pts) <= TOL_PX
