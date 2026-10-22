"""Characterization tests for StackFlats.inv_median (task T.3).

inv_median(a) = 1 / median(a), the scaling used when combining flats.
All expected values are hand-computed.
"""
import numpy as np
import pytest

from kcwidrp.primitives.StackFlats import inv_median


def test_odd_length():
    # source: analytic - median of [5, 1, 4] is 4, inverse 0.25 (exact)
    assert inv_median(np.array([5.0, 1.0, 4.0])) == 0.25


def test_even_length_averages_middle_pair():
    # source: analytic - median of [1, 2, 6, 10] is (2 + 6)/2 = 4, inverse 0.25
    assert inv_median(np.array([10.0, 1.0, 6.0, 2.0])) == 0.25


def test_two_dimensional_uses_all_elements():
    # source: analytic - 2x2 array 2,4,8,16: median (4 + 8)/2 = 6 -> 1/6
    # Tolerance: 1/6 is not representable; compare to the correctly rounded
    # value 1/6 within one ulp (2.8e-17).
    assert inv_median(np.array([[2.0, 4.0], [8.0, 16.0]])) == \
        pytest.approx(1.0 / 6.0, rel=2e-16, abs=0)


def test_outlier_does_not_move_median():
    # source: analytic - an outlier of 1e9 leaves median 2 -> 0.5
    assert inv_median(np.array([2.0, 2.0, 2.0, 1e9, 2.0])) == 0.5


def test_zero_median_is_infinite():
    # source: analytic - the median of [0, 0, 7] is 0; 1.0/np.float64(0)
    # is inf (with a numpy warning)
    with np.errstate(divide="ignore"):
        assert inv_median(np.array([0.0, 0.0, 7.0])) == np.inf
