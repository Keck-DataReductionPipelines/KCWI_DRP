"""Characterization tests for kcwidrp.core.kcwi_plotting.get_plot_lims
(task T.3).  Hand-computed; these survive the Track A plotting rewrite as long
as get_plot_lims keeps its contract: (min - pad*range, max + pad*range).

Values are chosen so every result is exactly representable in binary floating
point, so the comparisons are exact.
"""
import numpy as np

from kcwidrp.core.kcwi_plotting import get_plot_lims


def test_default_padding_is_five_percent():
    # source: analytic - range 20, pad 0.05*20 = 1: (0 - 1, 20 + 1)
    assert get_plot_lims(np.array([0.0, 20.0, 5.0])) == (-1.0, 21.0)


def test_explicit_padding():
    # source: analytic - range 8, pad 0.25*8 = 2: (0 - 2, 8 + 2)
    assert get_plot_lims([0.0, 8.0, 4.0], padding=0.25) == (-2.0, 10.0)


def test_zero_padding_gives_data_range():
    # source: analytic - no padding: exactly (min, max)
    assert get_plot_lims(np.array([3.0, -7.0, 1.0]), padding=0.0) == (-7.0, 3.0)


def test_negative_values_and_order_independent():
    # source: analytic - range 4 (from -6 to -2), pad 0.5*4 = 2: (-8, 0)
    assert get_plot_lims(np.array([-2.0, -6.0, -4.0]), padding=0.5) == (-8.0, 0.0)


def test_nan_values_are_ignored():
    # source: analytic - NaNs do not affect the range: data 0..8, pad 2
    d = np.array([np.nan, 0.0, 8.0, np.nan])
    assert get_plot_lims(d, padding=0.25) == (-2.0, 10.0)


def test_two_dimensional_input_uses_whole_array():
    # source: analytic - range over all elements: 2..10, pad 0.25*8 = 2
    d = np.array([[2.0, 3.0], [10.0, 4.0]])
    assert get_plot_lims(d, padding=0.25) == (0.0, 12.0)


def test_constant_data_clipped_to_unit_range():
    # source: analytic - a range under 1e-4 with clip=True is replaced by 1,
    # so with padding 0.5 the limits are c -/+ 0.5
    assert get_plot_lims(np.array([4.0, 4.0, 4.0]), padding=0.5) == (3.5, 4.5)


def test_constant_data_without_clip_is_degenerate():
    # source: analytic - clip=False keeps the zero range: (c, c)
    assert get_plot_lims(np.array([4.0, 4.0]), padding=0.5, clip=False) == \
        (4.0, 4.0)


def test_tiny_range_is_also_clipped():
    # source: analytic - range 5e-5 < 1e-4 counts as no range: (0 - 0.5, 5e-5 + 0.5)
    lo, hi = get_plot_lims(np.array([0.0, 5e-5]), padding=0.5)
    assert lo == -0.5
    # Tolerance: 5e-5 + 0.5 is not exactly representable; the sum is rounded
    # to the nearest double (relative 1e-16).
    assert abs(hi - 0.50005) < 1e-15
