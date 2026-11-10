"""SubtractOverscan run through ``run_primitive`` (task T.4, detector level).

Specification (the class docstring and the plan): for each amplifier the
overscan offset is the median of the overscan section (BSECn, minus
``oscanbuf`` columns at each end), truncated to an integer, and it is
subtracted from the amplifier's data section (DSECn).  The overscan fit also
yields a read-noise estimate OSCNRNn.

Inputs are built so the right answer follows from arithmetic, not from the
code: each amplifier's overscan rectangle holds a known constant or a linear
ramp, the data rectangle a known pattern, and everything else is left alone.
The section rectangles are parsed here from the header strings (1-based,
inclusive, possibly descending), independently of ``parse_imsec``.
"""
import re

import numpy as np
import pytest

from kcwidrp.primitives import SubtractOverscan as subtract_overscan_module
from kcwidrp.primitives.SubtractOverscan import SubtractOverscan

# (camera, ampmode, binning, amp numbers as the header numbers them)
MODES = [
    ("BLUE", "ALL", (2, 2), (1, 2, 3, 4)),
    ("BLUE", "ALL", (1, 1), (1, 2, 3, 4)),
    ("BLUE", "TUP", (2, 2), (1, 2)),
    ("BLUE", "TUP", (1, 1), (1, 2)),
    ("RED", "L2U2L1U1", (2, 2), (0, 1, 2, 3)),
    ("RED", "L2U2L1U1", (1, 1), (0, 1, 2, 3)),
    ("RED", "L2U2", (2, 2), (1, 3)),
    ("RED", "L2U2", (1, 1), (1, 3)),
]
pytestmark = pytest.mark.filterwarnings(
    "ignore:Polyfit may be poorly conditioned")   # order-7 fit, see below

MODE_IDS = ["%s-%s-%dx%d" % (c, a, b[0], b[1]) for c, a, b, _ in MODES]


def rect(header, key):
    """(row slice, column slice) of a '[x0:x1,y0:y1]' header section.

    1-based inclusive FITS ranges (either direction, optional ':stride' on
    RED) become 0-based python slices.
    """
    m = re.fullmatch(r"\[(\d+):(\d+)(?::-?\d+)?,(\d+):(\d+)(?::-?\d+)?\]",
                     header[key])
    x0, x1, y0, y1 = map(int, m.groups())
    return (slice(min(y0, y1) - 1, max(y0, y1)),
            slice(min(x0, x1) - 1, max(x0, x1)))


def build(synthetic_ccddata, camera, ampmode, binning, amps, *, overscan_fn,
          pattern_base=5000.0):
    """Frame with a known overscan and data pattern per amplifier.

    ``overscan_fn(k, nrows, ncols)`` returns the amplifier's overscan block.
    Returns (ccd, expected_osval list, rows slices, dcols slices).
    """
    ccd = synthetic_ccddata(camera=camera, ampmode=ampmode, binning=binning)
    hdr = ccd.header
    data = np.full(ccd.data.shape, 77.0)         # unread/gap rows, untouched
    expected = data.copy()
    info = []
    for k, ia in enumerate(amps):
        brow, bcol = rect(hdr, "BSEC%d" % ia)
        drow, dcol = rect(hdr, "DSEC%d" % ia)
        assert brow == drow                       # layout sanity
        nrows = brow.stop - brow.start
        ncols = bcol.stop - bcol.start
        block = overscan_fn(k, nrows, ncols)
        data[brow, bcol] = block
        yy = np.arange(nrows)[:, None]
        patt = pattern_base + 10.0 * k + (yy % 5) + np.zeros(
            (1, dcol.stop - dcol.start))
        data[drow, dcol] = patt
        info.append((ia, brow, bcol, drow, dcol, patt))
    ccd.data = data
    return ccd, info, expected


def check_subtracted(out, before, info, osvals):
    """Data rectangles lowered by osval, every other pixel unchanged."""
    expected = before.copy()
    for (ia, brow, bcol, drow, dcol, patt), osval in zip(info, osvals):
        expected[drow, dcol] = patt - osval
    assert np.array_equal(out.ccddata.data, expected)


@pytest.mark.parametrize("camera,ampmode,binning,amps", MODES, ids=MODE_IDS)
def test_constant_overscan_removed_per_amp(run_primitive, synthetic_ccddata,
                                           camera, ampmode, binning, amps):
    # source: independent @ cf88fb8 - the overscan of amp k is the constant
    # 1000 + 25*k, so the median, and the offset to subtract, is that
    # constant; each data rectangle ends up at pattern - constant.
    consts = [1000.0 + 25.0 * k for k in range(len(amps))]
    ccd, info, _ = build(
        synthetic_ccddata, camera, ampmode, binning, amps,
        overscan_fn=lambda k, nr, nc: np.full((nr, nc), consts[k]))
    before = ccd.data.copy()
    out = run_primitive(SubtractOverscan, ccd)
    check_subtracted(out, before, info, consts)
    hdr = out.ccddata.header
    for k, ia in enumerate(amps):
        # source: independent @ cf88fb8 - OSCNVALn records the offset removed
        assert hdr["OSCNVAL%d" % ia] == consts[k]
        # source: analytic @ cf88fb8 - a constant overscan leaves no
        # residual about its polynomial fit, so the noise estimate is 0.
        # Tolerance: the order-7 polyfit of the 2-amp modes is ill
        # conditioned (RankWarning), leaving rounding residuals of a few
        # 1e-3 e- (seen: 0.003); 0.01 e- is 1/300 of a typical 3 e- noise.
        assert abs(hdr["OSCNRN%d" % ia]) < 0.01
    assert hdr["OSCANSUB"] is True      # source: analytic @ cf88fb8


def _sloped_case(run_primitive, synthetic_ccddata, camera, ampmode, binning,
                 amps, row_slope, col_slope):
    # Amp k's overscan block is v(r, c) = base_k + row_slope r + col_slope c
    # (r, c = 0-based row / column offsets in the block).  The offset to
    # subtract is the median of the block with `oscanbuf` columns removed at
    # both ends, over all rows.  An arithmetic progression is symmetric about
    # its centre, so the median over a symmetric set of rows and columns is
    # the mean: base + row_slope (nr-1)/2 + col_slope (nc-1)/2.
    bases = [900.0 + 40.0 * k for k in range(len(amps))]

    def ramp(k, nr, nc):
        r = np.arange(nr)[:, None]
        c = np.arange(nc)[None, :]
        return bases[k] + row_slope * r + col_slope * c

    ccd, info, _ = build(synthetic_ccddata, camera, ampmode, binning, amps,
                         overscan_fn=ramp)
    before = ccd.data.copy()
    out = run_primitive(SubtractOverscan, ccd)
    osvals = []
    for k, (ia, brow, bcol, drow, dcol, patt) in enumerate(info):
        nr = brow.stop - brow.start
        nc = bcol.stop - bcol.start
        osvals.append(bases[k] + row_slope * (nr - 1) / 2.0
                      + col_slope * (nc - 1) / 2.0)
    return out, before, info, osvals


@pytest.mark.parametrize("camera,ampmode,binning,amps", MODES, ids=MODE_IDS)
def test_overscan_sloped_along_rows_removed_per_amp(
        run_primitive, synthetic_ccddata, camera, ampmode, binning, amps):
    # source: analytic @ cf88fb8 - slope of 2 DN per row, none along columns:
    # offset = base + 2 (nr-1)/2 = base + (nr-1), an integer, so the
    # truncation to int changes nothing.
    out, before, info, osvals = _sloped_case(
        run_primitive, synthetic_ccddata, camera, ampmode, binning, amps,
        row_slope=2.0, col_slope=0.0)
    check_subtracted(out, before, info, osvals)
    for k, ia in enumerate(amps):
        assert out.ccddata.header["OSCNVAL%d" % ia] == osvals[k]
        # source: analytic @ cf88fb8 - the per-row overscan vector is a
        # straight line, inside the fit's polynomial family (order 2 or 7),
        # so the residual about the fit is rounding noise.
        # Tolerance: as above, 0.01 e- covers the ill-conditioned order-7 fit.
        assert abs(out.ccddata.header["OSCNRN%d" % ia]) < 0.01


@pytest.mark.xfail(strict=True, reason="B29: SubtractOverscan.py:57-58 "
                   "overscan window is [bsec_x0+oscanbuf : bsec_x1-oscanbuf) "
                   "with bsec_x1 inclusive, so the margin is oscanbuf columns "
                   "at the start but oscanbuf+1 at the end")
@pytest.mark.parametrize("camera,ampmode,binning,amps", MODES, ids=MODE_IDS)
def test_overscan_sloped_along_columns_uses_symmetric_margin(
        run_primitive, synthetic_ccddata, camera, ampmode, binning, amps):
    # source: analytic @ cf88fb8 - slope of 2 DN per column, none along rows:
    # with `oscanbuf` columns excluded at both ends the window is symmetric
    # about the block centre and the median is base + 2 (nc-1)/2.  The code
    # excludes one column more at the far end, centring the window half a
    # column low, so it subtracts 1 DN less than this value (slope 2 DN per
    # column).
    out, before, info, osvals = _sloped_case(
        run_primitive, synthetic_ccddata, camera, ampmode, binning, amps,
        row_slope=0.0, col_slope=2.0)
    check_subtracted(out, before, info, osvals)


def test_non_integer_offset_is_truncated(run_primitive, synthetic_ccddata):
    # source: characterization @ cf88fb8 - the offset is taken as
    # int(median): an overscan of 1000.75 DN subtracts 1000, not 1000.75.
    # (A reduction therefore keeps up to 1 DN of per-amp bias; the plan does
    # not list it as a bug.)
    amps = (1, 2, 3, 4)
    ccd, info, _ = build(
        synthetic_ccddata, "BLUE", "ALL", (2, 2), amps,
        overscan_fn=lambda k, nr, nc: np.full((nr, nc), 1000.75))
    before = ccd.data.copy()
    out = run_primitive(SubtractOverscan, ccd)
    check_subtracted(out, before, info, [1000.0] * 4)
    assert out.ccddata.header["OSCNVAL1"] == 1000


def test_narrow_overscan_is_left_alone(run_primitive, synthetic_ccddata):
    # source: analytic @ cf88fb8 - the docstring-level rule: with an overscan
    # no wider than minoscanpix (75 for BLUE) the amp is not fitted, so no
    # pixel changes and OSCANSUB is False.  A 60 px overscan is below 75.
    ccd = synthetic_ccddata(camera="BLUE", binning=(2, 2), overscan=60)
    before = ccd.data.copy()
    out = run_primitive(SubtractOverscan, ccd)
    assert np.array_equal(out.ccddata.data, before)
    assert out.ccddata.header["OSCANSUB"] is False
    assert "OSCNVAL1" not in out.ccddata.header


def test_plot_level_3_with_viewer_off(run_primitive, synthetic_ccddata,
                                      monkeypatch):
    # source: independent @ cf88fb8 - B10 (plan s13): at plot_level 3 with
    # bokeh_session None the primitive must run to the end.  It asks
    # "Next?" once per amp and once after the grid (4 amps -> 5 answers).
    # save_plot is replaced because it needs a browser to export a PNG.
    saved = []
    monkeypatch.setattr(subtract_overscan_module, "save_plot",
                        lambda fig, filename=None: saved.append(filename))
    ccd, info, _ = build(
        synthetic_ccddata, "BLUE", "ALL", (2, 2), (1, 2, 3, 4),
        overscan_fn=lambda k, nr, nc: np.full((nr, nc), 1000.0))
    before = ccd.data.copy()
    out = run_primitive(SubtractOverscan, ccd,
                        config_overrides={"plot_level": 3},
                        answers=[""] * 5)
    assert run_primitive.last.context.bokeh_session is None
    check_subtracted(out, before, info, [1000.0] * 4)
    # source: independent @ cf88fb8 - filename "oscan_%05d.png" % FRAMENO
    assert saved == ["oscan_%05d.png" % ccd.header["FRAMENO"]]
