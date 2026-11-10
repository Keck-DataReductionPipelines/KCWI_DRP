"""Characterization tests for MakeMasterFlat.bm_ledge_position (task T.3).

bm_ledge_position(cwave, dich) is a linear fit of the pixel/Angstrom position
of the "ledge" seen in BM-grating flats against the central wavelength, with
a different intercept with and without the dichroic.  No independent source
for the coefficients exists, so the table below is a characterization of the
code at b0e2550 (hand-checked: 0.240742 * 4000 = 962.968 and so on).
"""
import pytest

from kcwidrp.primitives.MakeMasterFlat import bm_ledge_position

# (cwave, dich, expected)
# source: characterization @ b0e2550
TABLE = [
    (0.0, False, 4044.56),
    (0.0, True, 4060.56),
    (4000.0, False, 5007.528),
    (4000.0, True, 5023.528),
    (4500.0, False, 5127.899),
    (4500.0, True, 5143.899),
    (5000.0, False, 5248.27),
    (5000.0, True, 5264.27),
    (4600.5, False, 5152.093571),
    (4600.5, True, 5168.093571),
]


@pytest.mark.parametrize("cwave, dich, expected", TABLE)
def test_ledge_position_table(cwave, dich, expected):
    # Tolerance: the expected values are the exact decimal products/sums; the
    # code adds float64 numbers that are not exactly representable, so the
    # last bit differs (~1e-12 at 5e3); 1e-9 absolute.
    assert bm_ledge_position(cwave, dich) == pytest.approx(expected, rel=0,
                                                           abs=1e-9)


def test_dichroic_shifts_ledge_by_sixteen():
    # source: characterization @ b0e2550 - the two intercepts differ by 16
    assert bm_ledge_position(4500.0, True) - \
        bm_ledge_position(4500.0, False) == pytest.approx(16.0, abs=1e-9)


def test_ledge_moves_with_cwave_at_fixed_slope():
    # source: characterization @ b0e2550 - slope 0.240742 per Angstrom
    d = bm_ledge_position(5000.0, False) - bm_ledge_position(4000.0, False)
    assert d == pytest.approx(240.742, abs=1e-9)
