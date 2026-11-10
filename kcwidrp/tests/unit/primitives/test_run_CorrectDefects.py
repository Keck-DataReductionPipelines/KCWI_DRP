"""CorrectDefects run through ``run_primitive`` (task T.4, detector level).

Specification (class docstring): a defect list is looked up in the package
data directory from the CCD amplifier mode and the x and y binning
(``data/defect_<AMPMODE>_<xbin>x<ybin>.dat``); every listed bad area is
cleaned and flagged; BPFILE records the list used and NBPCLEAN the number of
pixels cleaned.  A defect file has the columns ``X0 X1 Y0 Y1``: 1-based,
inclusive column and row ranges.

Three groups of tests:

* the choice of file, with the package-data lookup replaced by a recorder so
  the requested name is seen directly;
* what is cleaned and flagged, with a small defect file written to tmp_path;
* the real data files, with the expected pixel sets built in the test from an
  independent parse of the file (``numpy.loadtxt``).
"""
import numpy as np
import pytest

from kcwidrp.primitives import CorrectDefects as defects_module
from kcwidrp.primitives.CorrectDefects import CorrectDefects


class Lookup:
    """Stands in for get_resource_path: records requests, maps to tmp files."""

    def __init__(self, tmp_path, files=None):
        self.requested = []
        self.tmp_path = tmp_path
        self.files = files or {}       # relative name -> text

    def __call__(self, package, path):
        self.requested.append((package, path))
        target = self.tmp_path / ("lookup_" + path.replace("/", "_"))
        if path in self.files:
            target.write_text(self.files[path])
        return str(target)


@pytest.fixture
def lookup(monkeypatch, tmp_path):
    lk = Lookup(tmp_path)
    monkeypatch.setattr(defects_module, "get_resource_path", lk)
    return lk


# -- which defect file ---------------------------------------------------------

@pytest.mark.parametrize("camera,ampmode,binning", [
    ("BLUE", "ALL", (1, 1)), ("BLUE", "ALL", (2, 2)),
    ("BLUE", "TUP", (1, 1)), ("BLUE", "TUP", (2, 2)),
    ("RED", "L2U2L1U1", (1, 1)), ("RED", "L2U2L1U1", (2, 2)),
    ("RED", "L2U2", (1, 1)), ("RED", "L2U2", (2, 2)),
])
def test_file_name_from_ampmode_and_binning(run_primitive, synthetic_ccddata,
                                            lookup, camera, ampmode,
                                            binning):
    ccd = synthetic_ccddata(camera=camera, ampmode=ampmode, binning=binning)
    run_primitive(CorrectDefects, ccd)
    # source: analytic @ cf88fb8 - docstring: file chosen by ampmode and the
    # x and y binning; the package is kcwidrp.
    assert lookup.requested == [
        ("kcwidrp", "data/defect_%s_%dx%d.dat" % (ampmode, binning[0],
                                                  binning[1]))]


def test_unequal_binning_in_name(run_primitive, synthetic_ccddata, lookup):
    # source: analytic @ cf88fb8 - x and y binning enter separately, in that
    # order (the fixture's frame is only read for its header).
    ccd = synthetic_ccddata(camera="BLUE", ampmode="ALL", binning=(1, 2))
    run_primitive(CorrectDefects, ccd)
    assert lookup.requested[0][1] == "data/defect_ALL_1x2.dat"


def test_ampmode_whitespace_is_stripped(run_primitive, synthetic_ccddata,
                                        lookup):
    # source: characterization @ cf88fb8 - an AMPMODE card padded with blanks
    # ("TBO ") still selects defect_TBO_...; the code calls .strip().
    ccd = synthetic_ccddata(camera="BLUE", ampmode="ALL", AMPMODE="TBO  ")
    run_primitive(CorrectDefects, ccd)
    assert lookup.requested[0][1] == "data/defect_TBO_2x2.dat"


@pytest.mark.parametrize("nasname,numopen,suffix", [
    ("Mask", 2, "_nas"), ("Mask", 1, ""), ("Open", 2, ""), ("Open", 1, ""),
])
def test_nod_and_shuffle_suffix(run_primitive, synthetic_ccddata, lookup,
                                nasname, numopen, suffix):
    # source: characterization @ cf88fb8 - the "_nas" suffix is not in the
    # docstring; the code adds it when the nod-and-shuffle mask is in
    # (BNASNAM contains 'Mask') and more than one shuffle window is open
    # (NUMOPEN > 1).
    ccd = synthetic_ccddata(camera="BLUE", ampmode="ALL", BNASNAM=nasname,
                            NUMOPEN=numopen)
    run_primitive(CorrectDefects, ccd)
    assert lookup.requested[0][1] == "data/defect_ALL_2x2%s.dat" % suffix


def test_missing_defect_file_is_reported_not_fatal(run_primitive,
                                                   synthetic_ccddata, lookup):
    # source: characterization @ cf88fb8 - no such list: nothing is changed,
    # BPCLEAN False, NBPCLEAN 0, BPFILE left as it was.
    ccd = synthetic_ccddata(camera="BLUE", ampmode="ALL", BNASNAM="Mask",
                            NUMOPEN=2)
    before = ccd.data.copy()
    bpfile = ccd.header.get("BPFILE")
    out = run_primitive(CorrectDefects, ccd)
    assert np.array_equal(out.ccddata.data, before)
    assert not out.ccddata.flags.any()
    hdr = out.ccddata.header
    assert hdr["BPCLEAN"] is False
    assert hdr["NBPCLEAN"] == 0
    assert hdr.get("BPFILE") == bpfile


# -- what is cleaned -----------------------------------------------------------

def test_cleans_and_flags_listed_pixels(run_primitive, synthetic_ccddata,
                                        monkeypatch, tmp_path):
    lk = Lookup(tmp_path, {"data/defect_ALL_2x2.dat":
                           "X0 X1 Y0 Y1\n10 12 5 6\n20 20 30 30\n"})
    monkeypatch.setattr(defects_module, "get_resource_path", lk)
    ccd = synthetic_ccddata(camera="BLUE", ampmode="ALL")
    ny, nx = ccd.data.shape
    rr, cc = np.mgrid[0:ny, 0:nx]
    ccd.data = 100.0 + 3.0 * cc + 1000.0 * rr          # known, varies per row
    original = ccd.data.copy()
    out = run_primitive(CorrectDefects, ccd)
    # source: independent @ cf88fb8 - defect file rows (1-based, inclusive):
    #   cols 10-12, rows 5-6  -> 0-based rows 4-5, cols 9-11 (6 pixels)
    #   col 20, row 30        -> 0-based row 29, col 19 (1 pixel)
    bad = np.zeros((ny, nx), dtype=bool)
    bad[4:6, 9:12] = True
    bad[29, 19] = True
    expected_flags = np.where(bad, 2, 0).astype(np.uint8)
    assert np.array_equal(out.ccddata.flags, expected_flags)
    # source: characterization @ cf88fb8 - a bad pixel is replaced by the
    # median of the two good pixels on each side in its row.  On the ramp
    # 100 + 3 col + 1000 row the neighbours at cols 7,8,12,13 (0-based) have
    # median 100 + 3 (8 + 12)/2 + 1000 row = 130 + 1000 row; for the single
    # pixel at col 19 the neighbours are cols 17,18,20,21: median
    # 100 + 3 (18 + 20)/2 + 1000*29 = 157 + 29000.
    expected = original.copy()
    expected[4, 9:12] = 130.0 + 1000.0 * 4
    expected[5, 9:12] = 130.0 + 1000.0 * 5
    expected[29, 19] = 157.0 + 1000.0 * 29
    assert np.array_equal(out.ccddata.data, expected)
    hdr = out.ccddata.header
    assert hdr["NBPCLEAN"] == 7                    # source: independent
    assert hdr["BPCLEAN"] is True                  # source: independent
    assert hdr["BPFILE"] == "data/defect_ALL_2x2.dat"


def test_flags_are_added_to_existing_flags(run_primitive, synthetic_ccddata,
                                           monkeypatch, tmp_path):
    # source: characterization @ cf88fb8 - the defect flag value 2 is added
    # to what the pixel already carries (a saturated pixel, 8, becomes 10).
    lk = Lookup(tmp_path, {"data/defect_ALL_2x2.dat":
                           "X0 X1 Y0 Y1\n10 10 5 5\n"})
    monkeypatch.setattr(defects_module, "get_resource_path", lk)
    ccd = synthetic_ccddata(camera="BLUE", ampmode="ALL")
    ccd.flags[4, 9] = 8
    ccd.flags[0, 0] = 8
    out = run_primitive(CorrectDefects, ccd)
    assert out.ccddata.flags[4, 9] == 10
    assert out.ccddata.flags[0, 0] == 8


def test_flags_created_when_absent(run_primitive, synthetic_ccddata,
                                   monkeypatch, tmp_path):
    # source: characterization @ cf88fb8 - a frame without a flags array gets
    # a uint8 zero array of the data shape.
    lk = Lookup(tmp_path, {"data/defect_ALL_2x2.dat":
                           "X0 X1 Y0 Y1\n10 10 5 5\n"})
    monkeypatch.setattr(defects_module, "get_resource_path", lk)
    ccd = synthetic_ccddata(camera="BLUE", ampmode="ALL")
    ccd.flags = None
    out = run_primitive(CorrectDefects, ccd)
    assert out.ccddata.flags.dtype == np.uint8
    assert out.ccddata.flags.shape == ccd.data.shape
    assert out.ccddata.flags.sum() == 2
    assert out.ccddata.flags[4, 9] == 2


# -- the real defect files -----------------------------------------------------

# camera, ampmode, binning, file (None: the package has no such list)
REAL = [
    ("BLUE", "ALL", (1, 1), "defect_ALL_1x1.dat"),
    ("BLUE", "ALL", (2, 2), "defect_ALL_2x2.dat"),
    ("BLUE", "TUP", (1, 1), "defect_TUP_1x1.dat"),
    ("BLUE", "TUP", (2, 2), "defect_TUP_2x2.dat"),
    ("RED", "L2U2", (1, 1), "defect_L2U2_1x1.dat"),
    ("RED", "L2U2", (2, 2), "defect_L2U2_2x2.dat"),
    ("RED", "L2U2L1U1", (1, 1), "defect_L2U2L1U1_1x1.dat"),
    ("RED", "L2U2L1U1", (2, 2), None),
]
# Pixels listed in each file: sum of (X1-X0+1)*(Y1-Y0+1), by hand from the
# files at cf88fb8.
REAL_COUNT = {
    "defect_ALL_1x1.dat": 189,          # 1 x (2245-2057+1)
    "defect_ALL_2x2.dat": 94,           # 1 x (1122-1029+1)
    "defect_TUP_1x1.dat": 2319,         # 37 + 6*6 + 2246
    "defect_TUP_2x2.dat": 1182,         # 44 + 5*3 + 1123
    "defect_L2U2_1x1.dat": 3070,        # 5 x 614
    "defect_L2U2_2x2.dat": 306,         # 1 x 306
    "defect_L2U2L1U1_1x1.dat": 1457,    # 1 x 1457
}


@pytest.mark.parametrize("camera,ampmode,binning,fname", REAL,
                         ids=["%s-%s-%dx%d" % (c, a, b[0], b[1])
                              for c, a, b, _ in REAL])
def test_real_defect_file_flags_exactly_its_pixels(
        run_primitive, synthetic_ccddata, camera, ampmode, binning, fname):
    ccd = synthetic_ccddata(camera=camera, ampmode=ampmode, binning=binning)
    if fname is None:
        # source: characterization @ cf88fb8 - no defect list is shipped for
        # RED L2U2L1U1 at 2x2: nothing is cleaned.
        out = run_primitive(CorrectDefects, ccd)
        assert out.ccddata.header["BPCLEAN"] is False
        assert out.ccddata.header["NBPCLEAN"] == 0
        return
    from kcwidrp.core.kcwi_pkg_resources import get_resource_path
    table = np.loadtxt(get_resource_path("kcwidrp", "data/" + fname),
                       skiprows=1, dtype=int, ndmin=2)
    ny = int(table[:, 3].max()) + 8
    nx = int(table[:, 1].max()) + 8
    ccd.data = np.ones((ny, nx), dtype=np.float32)   # trimmed-frame size
    ccd.flags = np.zeros((ny, nx), dtype=np.uint8)
    # source: independent @ cf88fb8 - X0 X1 Y0 Y1 are 1-based inclusive
    expected = np.zeros((ny, nx), dtype=np.uint8)
    for x0, x1, y0, y1 in table:
        expected[y0 - 1:y1, x0 - 1:x1] = 2
    out = run_primitive(CorrectDefects, ccd)
    assert np.array_equal(out.ccddata.flags, expected)
    assert out.ccddata.header["NBPCLEAN"] == REAL_COUNT[fname]
    assert int(expected.sum() // 2) == REAL_COUNT[fname]   # hand count agrees
    assert out.ccddata.header["BPFILE"] == "data/" + fname
    assert out.ccddata.header["BPCLEAN"] is True
