"""Characterization tests for FitCenter.pascal_shift (task T.3).

pascal_shift(coefficients, x0) takes polynomial coefficients in numpy's
*descending* order (highest power first, as np.polyfit returns) and returns
the coefficients of the same polynomial after the reference point is moved
to x0, q(t) = p(t - x0)  [pascal_shift's own docstring: "shift the polynomial
coefficients to a new reference point"].

Expected values: numpy.polynomial's composition P(coefs)(t - x0) -- an
independent implementation.
"""
import numpy as np
import pytest
from numpy.polynomial import Polynomial

from kcwidrp.primitives.FitCenter import pascal_shift


def _numpy_shift(desc_coefs, x0):
    # source: independent - numpy.polynomial composition p(t - x0)
    p = Polynomial(list(reversed(desc_coefs)))           # ascending order
    q = p(Polynomial([-x0, 1.0]))
    asc = list(q.coef) + [0.0] * (len(desc_coefs) - len(q.coef))
    return list(reversed(asc))


@pytest.mark.parametrize("coefs", [
    [2.0, 3.0],                                    # linear
    [1.0, -2.0, 0.5],                              # quadratic
    [0.5, 1.0, -2.0, 4.0],                         # cubic
    [1e-3, -2e-2, 0.3, -1.0, 7.0],                 # quartic
    [1e-5, 2e-4, -3e-3, 4e-2, -0.5, 6.0, 10.0],    # 7 coefficients (maximum)
])
@pytest.mark.parametrize("x0", [1.0, -2.5, 4000.0])
def test_matches_numpy_polynomial_shift(coefs, x0):
    got = pascal_shift(coefficients=list(coefs), x0=x0)
    expected = _numpy_shift(coefs, x0)
    assert len(got) == len(coefs)
    # Tolerance: both sides sum binomial-weighted products of coefficients and
    # powers of x0 in float64; terms can cancel, so a coefficient that is
    # nearly zero is only accurate to ~1e-16 of the largest term that went
    # into it (atol = 1e-12 * that scale, 4 orders above float64 epsilon times
    # the <=64 binomial weights); rtol 1e-9 for the rest.
    scale = max(abs(c) * max(1.0, abs(x0)) ** i
                for i, c in enumerate(reversed(coefs)))
    np.testing.assert_allclose(got, expected, rtol=1e-9, atol=1e-12 * scale)


def test_constant_term_is_polynomial_value_at_minus_x0():
    # source: analytic - q(0) = p(-x0).  p = 2t^2 - 3t + 1, x0 = 2 -> p(-2) =
    # 8 + 6 + 1 = 15, the last coefficient (descending order).
    got = pascal_shift([2.0, -3.0, 1.0], 2.0)
    assert got[-1] == 15.0


def test_hand_computed_quadratic():
    # source: analytic - p(t) = 2t^2 - 3t + 1, q(t) = p(t - 2) =
    # 2(t^2 - 4t + 4) - 3t + 6 + 1 = 2t^2 - 11t + 15
    assert pascal_shift([2.0, -3.0, 1.0], 2.0) == [2.0, -11.0, 15.0]


def test_linear_coefficient_unchanged_slope():
    # source: analytic - shifting a line keeps its slope:
    # 3t + 2 -> 3(t - 1) + 2 = 3t - 1
    assert pascal_shift([3.0, 2.0], 1.0) == [3.0, -1.0]


def test_zero_shift_returns_input():
    # source: analytic - x0 = 0 is the identity (the code warns and returns the
    # input unchanged)
    c = [1.0, 2.0, 3.0]
    assert pascal_shift(c, 0.0) == c


def test_no_coefficients_returns_none():
    # source: characterization @ b0e2550 - empty input prints an error and
    # returns None
    assert pascal_shift([], 1.0) is None
    assert pascal_shift(None, 1.0) is None
