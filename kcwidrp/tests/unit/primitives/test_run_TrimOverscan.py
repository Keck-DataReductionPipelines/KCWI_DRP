"""TrimOverscan run through ``run_primitive`` (task T.4, detector level).

Specification (class docstring): the image is cut to the data sections
(DSECn) of the amplifiers and reassembled in the physical amplifier layout;
the raw section keywords ASECn/BSECn/CSECn/DSECn (and TSECn on RED) are
replaced by ATSECn, the amplifier's section in the trimmed image.  A BIAS
frame is also written as ``*_intb.fits`` and added to the proc table.

The expected layout is written down here from the camera's amplifier
geometry (the fixture's block comment, taken from real headers): a "left"
amplifier fills the left half of the trimmed image and a "right" one the
right half; "lower"/"upper" amplifiers the lower/upper half of the rows,
"full" all rows.  The fixture's frames have 32 data columns per amplifier and
64 data rows per amplifier row (128 for the 2-amp modes).
"""
import re

import numpy as np
import pytest
from astropy.io import fits

from kcwidrp.primitives.TrimOverscan import TrimOverscan

ND = 32          # data columns per amp in the fixture's default frame
# camera, ampmode -> {amp: (side, band)}, and the trimmed frame (rows, cols)
LAYOUTS = {
    ("BLUE", "ALL"): ({1: ("left", "lower"), 2: ("right", "lower"),
                       3: ("left", "upper"), 4: ("right", "upper")}, 128),
    ("BLUE", "TUP"): ({1: ("left", "full"), 2: ("right", "full")}, 128),
    ("RED", "L2U2L1U1"): ({0: ("left", "upper"), 1: ("left", "lower"),
                           2: ("right", "upper"), 3: ("right", "lower")},
                          128),
    ("RED", "L2U2"): ({1: ("left", "full"), 3: ("right", "full")}, 128),
}
CASES = [(c, a, b) for (c, a) in LAYOUTS for b in ((1, 1), (2, 2))]
CASE_IDS = ["%s-%s-%dx%d" % (c, a, b[0], b[1]) for c, a, b in CASES]


def rect(header, key):
    """(row slice, column slice) of a '[x0:x1,y0:y1]' header section."""
    m = re.fullmatch(r"\[(\d+):(\d+)(?::-?\d+)?,(\d+):(\d+)(?::-?\d+)?\]",
                     header[key])
    x0, x1, y0, y1 = map(int, m.groups())
    return (slice(min(y0, y1) - 1, max(y0, y1)),
            slice(min(x0, x1) - 1, max(x0, x1)))


def slot(side, band, nrows):
    """Row and column slices of an amp in the trimmed image."""
    cols = slice(0, ND) if side == "left" else slice(ND, 2 * ND)
    half = nrows // 2
    rows = {"lower": slice(0, half), "upper": slice(half, nrows),
            "full": slice(0, nrows)}[band]
    return rows, cols


def filled(synthetic_ccddata, camera, ampmode, binning, **kw):
    """Frame whose DSEC rectangles hold a distinct ramp per amp.

    Returns (ccd, {amp: pattern array}); pixels outside the data sections
    are NaN so any leak of overscan into the output shows.
    """
    ccd = synthetic_ccddata(camera=camera, ampmode=ampmode, binning=binning,
                            **kw)
    data = np.full(ccd.data.shape, np.nan)
    flags = np.zeros(ccd.data.shape, dtype=np.uint8)
    patterns = {}
    for k, ia in enumerate(LAYOUTS[(camera, ampmode)][0]):
        drow, dcol = rect(ccd.header, "DSEC%d" % ia)
        shape = (drow.stop - drow.start, dcol.stop - dcol.start)
        patt = 1000.0 * (k + 1) + np.arange(shape[0] * shape[1]).reshape(
            shape)
        data[drow, dcol] = patt
        flags[drow.start, dcol.start] = 8 + k       # marks each amp's corner
        patterns[ia] = patt
    ccd.data = data
    ccd.flags = flags
    return ccd, patterns


@pytest.mark.parametrize("camera,ampmode,binning", CASES, ids=CASE_IDS)
def test_trim_shape_and_position(run_primitive, synthetic_ccddata, camera,
                                 ampmode, binning):
    ccd, patterns = filled(synthetic_ccddata, camera, ampmode, binning)
    layout, nrows = LAYOUTS[(camera, ampmode)]
    out = run_primitive(TrimOverscan, ccd)
    # source: independent @ cf88fb8 - two amps of 32 columns side by side,
    # 128 data rows in total (fixture geometry, see module docstring).
    assert out.ccddata.shape == (nrows, 2 * ND)
    assert out.ccddata.flags.shape == (nrows, 2 * ND)
    expected = np.full((nrows, 2 * ND), np.nan)
    exp_flags = np.zeros((nrows, 2 * ND), dtype=np.uint8)
    for k, (ia, (side, band)) in enumerate(layout.items()):
        rows, cols = slot(side, band, nrows)
        expected[rows, cols] = patterns[ia]
        exp_flags[rows.start, cols.start] = 8 + k
    # source: independent @ cf88fb8 - amp data land unflipped in their slot;
    # no NaN (overscan/prescan/gap) pixel enters the output.
    assert np.array_equal(out.ccddata.data, expected)
    # source: independent @ cf88fb8 - flags are trimmed with the data
    assert np.array_equal(out.ccddata.flags, exp_flags)
    # source: characterization @ cf88fb8 - the output array is float32
    assert out.ccddata.data.dtype == np.float32
    assert out.ccddata.flags.dtype == np.uint8


@pytest.mark.parametrize("camera,ampmode,binning", CASES, ids=CASE_IDS)
def test_trim_header_sections(run_primitive, synthetic_ccddata, camera,
                              ampmode, binning):
    ccd, _ = filled(synthetic_ccddata, camera, ampmode, binning)
    layout, nrows = LAYOUTS[(camera, ampmode)]
    out = run_primitive(TrimOverscan, ccd)
    hdr = out.ccddata.header
    for ia, (side, band) in layout.items():
        rows, cols = slot(side, band, nrows)
        # source: independent @ cf88fb8 - ATSECn is the amp's 1-based
        # '[x0:x1,y0:y1]' range in the trimmed image
        assert hdr["ATSEC%d" % ia] == "[%d:%d,%d:%d]" % (
            cols.start + 1, cols.stop, rows.start + 1, rows.stop)
        for stem in ("ASEC", "BSEC", "CSEC", "DSEC", "TSEC"):
            assert "%s%d" % (stem, ia) not in hdr   # source: docstring
    assert hdr["NAXIS1"] == 2 * ND                  # source: independent
    assert hdr["NAXIS2"] == nrows                   # source: independent
    assert hdr["OSCANTRM"] is True                  # source: docstring


def test_bias_trim_writes_intb_and_proctab(run_primitive, synthetic_ccddata,
                                           tmp_redux):
    # source: independent @ cf88fb8 - docstring: "If the input image is a
    # bias frame, writes out a *_intb.fits file" (here under redux/, the
    # configured output_directory) and the frame enters the proc table with
    # TYPE BIAS and suffix intb.
    ccd, _ = filled(synthetic_ccddata, "BLUE", "ALL", (2, 2), IMTYPE="BIAS")
    out = run_primitive(TrimOverscan, ccd)
    written = tmp_redux.redux / "kb260101_00001_intb.fits"
    assert written.is_file()
    with fits.open(written) as hdul:
        assert hdul[0].data.shape == (128, 64)
    tab = tmp_redux.proctab.proctab
    assert len(tab) == 1
    assert tab["TYPE"][0] == "BIAS"
    assert tab["SUFF"][0] == "intb"
    assert tmp_redux.proctab_file.is_file()      # proc table written
    assert out.ccddata.header["IMTYPE"] == "BIAS"


def test_non_bias_trim_writes_nothing(run_primitive, synthetic_ccddata,
                                      tmp_redux):
    # source: independent @ cf88fb8 - docstring: only biases are written
    # (and saveintims is False by default)
    ccd, _ = filled(synthetic_ccddata, "BLUE", "ALL", (2, 2))
    run_primitive(TrimOverscan, ccd)
    assert list(tmp_redux.redux.iterdir()) == []
    assert len(tmp_redux.proctab.proctab) == 0
