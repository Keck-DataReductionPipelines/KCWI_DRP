"""Integer header values compare exactly in tier 4 (written from the spec only).

Two Python integers (not bool) agree only if equal, whatever rtol/atol. The
tolerance |va - vb| <= atol + rtol*|vb| applies only when at least one value is
a float; two NaN floats agree; everything else is compared with ==.
"""
import pytest
from astropy.io import fits

from kcwidrp.tests.regression.compare import compare_outputs, numbers_agree

PROC = "kcwib.proc"
TOL = {"*_icubes.fits": {"rtol": 1e-3, "atol": 0.0, "reason": "test"}}


@pytest.mark.parametrize("va, vb, rtol, atol, expected", [
    (4096, 4097, 1e-3, 0, False),         # source: spec, ints are exact
    (4096, 4096, 1e-3, 0, True),          # source: spec, equal ints
    (4096.0, 4097.0, 1e-3, 0, True),      # source: |1| <= 1e-3*4097 = 4.097
    (4096.0, 4106.0, 1e-3, 0, False),     # source: |10| > 1e-3*4106 = 4.106
    (4096, 4097.0, 1e-3, 0, True),        # source: int vs float, 1 <= 4.097
    (4096, 4106.0, 1e-3, 0, False),       # source: int vs float, 10 > 4.106
    (4096.0, 4097, 1e-3, 0, True),        # source: float vs int, 1 <= 4.097
    (4096.0, 4106, 1e-3, 0, False),       # source: float vs int, 10 > 4.106
    (5, 6, 0.0, 10, False),               # source: spec, atol ignored for ints
    (5, 6, 1.0, 10, False),               # source: spec, rtol ignored for ints
    (5.0, 6.0, 0.0, 10, True),            # source: |1| <= 10, floats
    (float("nan"), float("nan"), 0.0, 0.0, True),  # source: spec, NaN agree
    ("abc", "abc", 1e-3, 0, True),        # source: spec, == on strings
    ("abc", "abd", 1e-3, 0, False),       # source: spec, == on strings
    (True, True, 1e-3, 0, True),          # source: spec, bools use ==
    (True, False, 1e-3, 0, False),        # source: spec, bools use ==
])
def test_numbers_agree(va, vb, rtol, atol, expected):
    assert numbers_agree(va, vb, rtol, atol) is expected


def make_tree(root, **cards):
    (root / "redux").mkdir(parents=True)
    (root / "plots").mkdir()
    hdu = fits.PrimaryHDU()
    for key, val in cards.items():
        hdu.header[key] = val
    hdu.writeto(root / "redux" / "x_icubes.fits")
    (root / PROC).write_text("| FRAMENO | TYPE | filename |\n"
                             "|       1 | BIAS | KB.1.fits |\n")


def compare(tmp_path, run_cards, ref_cards):
    make_tree(tmp_path / "run", **run_cards)
    make_tree(tmp_path / "ref", **ref_cards)
    return compare_outputs(tmp_path / "run", tmp_path / "ref", PROC,
                           tolerances=TOL)


def test_integer_card_off_by_one_fails_tier4_naming_keyword(tmp_path):
    out = compare(tmp_path, {"NPLANES": 4097}, {"NPLANES": 4096})
    assert len(out.headers) == 1                  # source: spec, ints exact
    assert "NPLANES" in out.headers[0]
    assert out.structure == [] and out.pixels == []


def test_equal_integer_card_passes(tmp_path):
    out = compare(tmp_path, {"NPLANES": 4096}, {"NPLANES": 4096})
    assert out.headers == []                      # source: analytic


def test_float_card_within_tolerance_passes(tmp_path):
    out = compare(tmp_path, {"CWAVE": 4097.0}, {"CWAVE": 4096.0})
    assert out.headers == []                      # source: 1 <= 1e-3*4096


def test_float_card_outside_tolerance_fails(tmp_path):
    out = compare(tmp_path, {"CWAVE": 4106.0}, {"CWAVE": 4096.0})
    assert len(out.headers) == 1                  # source: 10 > 1e-3*4096
    assert "CWAVE" in out.headers[0]
