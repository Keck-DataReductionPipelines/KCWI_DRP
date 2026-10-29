"""Path rule for header values and pickled strings (written from the spec only).

A string is a path only if it starts with "/" (or is a Windows absolute path)
or ends in a known file suffix; a slash alone (a unit, an OBJECT name) does
not make a path. Paths compare by file name, every other string exactly.
"""
import pickle

import pytest
from astropy.io import fits

from kcwidrp.tests.regression.compare import (_is_path, compare_outputs,
                                              diff_objects)

PROC = "kcwib.proc"
UNIT = "1e+16 erg / (Angstrom s cm2)"
OBJ = "DOME/Dark/32.5 sec."


@pytest.mark.parametrize("value, expected", [
    (UNIT, False),                                # source: spec example
    (OBJ, False),                                 # source: spec example
    ("data/defect_TUP_2x2.dat", True),            # source: spec example
    ("/s/sdata1400/kcwi4/2024apr06", True),       # source: spec example
    ("kb_arc.FITS", True),                        # source: spec, case-insens.
    ("redux/x_icube.fits", True),                 # source: spec, suffix
    ("C:\\data\\kb_arc.fits", True),              # source: spec, Windows
    ("C:\\data\\dir", True),                      # source: spec, Windows abs.
    ("", False),                                  # source: spec, not a path
    ("star", False),                              # source: spec, no / or suffix
    ("  /a/b  ", True),                           # source: spec, strip first
    (3.5, False),                                 # source: spec, non-string
    (None, False),                                # source: spec, non-string
])
def test_is_path(value, expected):
    assert _is_path(value) is expected


def make_tree(root, **cards):
    (root / "redux").mkdir(parents=True)
    (root / "plots").mkdir()
    hdu = fits.PrimaryHDU()
    hdu.header["BUNIT"] = UNIT
    hdu.header["OBJECT"] = OBJ
    for key, val in cards.items():
        hdu.header[key] = val
    hdu.writeto(root / "redux" / "x_icube.fits")
    (root / PROC).write_text("| FRAMENO | TYPE | filename |\n"
                             "|       1 | BIAS | KB.1.fits |\n")


def compare(tmp_path, run_cards, ref_cards):
    make_tree(tmp_path / "run", **run_cards)
    make_tree(tmp_path / "ref", **ref_cards)
    return compare_outputs(tmp_path / "run", tmp_path / "ref", PROC)


def test_identical_headers_have_no_differences(tmp_path):
    out = compare(tmp_path, {"ARCFL": "/a/kb_arc.fits"},
                  {"ARCFL": "/a/kb_arc.fits"})
    assert out.headers == [] and out.structure == []  # source: analytic


def test_bunit_change_fails_tier4_naming_bunit(tmp_path):
    out = compare(tmp_path, {"BUNIT": UNIT.replace("e+16", "e+17")}, {})
    assert len(out.headers) == 1                  # source: spec, one per card
    assert "BUNIT" in out.headers[0]
    assert out.structure == [] and out.pixels == []


def test_object_with_slashes_change_fails_tier4_naming_object(tmp_path):
    out = compare(tmp_path, {"OBJECT": "DOME/Dark/30.0 sec."}, {})
    assert len(out.headers) == 1                  # source: spec
    assert "OBJECT" in out.headers[0]


def test_absolute_path_in_another_directory_passes(tmp_path):
    out = compare(tmp_path, {"ARCFL": "/run/dir/kb_arc.fits"},
                  {"ARCFL": "/ref/other/kb_arc.fits"})
    assert out.headers == [] and out.structure == []  # source: spec


def test_windows_path_in_another_directory_passes(tmp_path):
    out = compare(tmp_path, {"ARCFL": "C:\\run\\kb_arc.fits"},
                  {"ARCFL": "/ref/kb_arc.fits"})
    assert out.headers == []                      # source: spec, \ as /


def test_path_with_different_file_name_fails_naming_arcfl(tmp_path):
    out = compare(tmp_path, {"ARCFL": "/run/dir/kb_arc.fits"},
                  {"ARCFL": "/run/dir/kb_other.fits"})
    assert len(out.headers) == 1                  # source: spec
    assert "ARCFL" in out.headers[0]


def test_relative_path_with_different_file_name_fails(tmp_path):
    out = compare(tmp_path, {"BPFILE": "data/defect_A_2x2.dat"},
                  {"BPFILE": "data/defect_B_2x2.dat"})
    assert len(out.headers) == 1                  # source: spec
    assert "BPFILE" in out.headers[0]


def test_ignored_cards_still_ignored(tmp_path):
    out = compare(tmp_path, {"HISTORY": "run one", "COMMENT": "a"},
                  {"HISTORY": "run two", "COMMENT": "b"})
    assert out.headers == [] and out.structure == []  # source: spec


def test_diff_objects_unit_string_with_slash_is_exact():
    a = {"unit": UNIT}
    b = {"unit": UNIT.replace("e+16", "e+17")}
    diffs = diff_objects(a, b)
    assert len(diffs) == 1 and "unit" in diffs[0]     # source: spec


def test_diff_objects_path_in_another_directory_is_equal():
    a = {"f": "/run/dir/kb_arc.fits"}
    b = {"f": "/ref/other/kb_arc.fits"}
    assert diff_objects(a, b) == []               # source: spec


def test_diff_objects_path_with_different_file_name_differs():
    diffs = diff_objects({"f": "/d/kb_arc.fits"}, {"f": "/d/kb_x.fits"})
    assert len(diffs) == 1 and "f" in diffs[0]    # source: spec
