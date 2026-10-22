"""FlagSaturation run through ``run_primitive`` (task T.4, detector level).

Specification (class docstring): pixels with more than 60,000 counts get the
value 8 in the flags array; SATFLAG becomes True and NSATFLAG the number of
saturated pixels flagged.
"""
import numpy as np
import pytest

from kcwidrp.primitives.FlagSaturation import FlagSaturation


def frame(synthetic_ccddata, **kw):
    ccd = synthetic_ccddata(**kw)
    ccd.data = np.full(ccd.data.shape, 1000.0)
    return ccd


def test_flags_set_exactly_at_saturated_pixels(run_primitive,
                                               synthetic_ccddata):
    ccd = frame(synthetic_ccddata)
    sat = [(3, 5), (10, 0), (20, 40), (-1, -1)]
    for r, c in sat:
        ccd.data[r, c] = 65535.0
    ccd.data[4, 4] = 60001.0     # just above the threshold
    ccd.data[5, 5] = 60000.0     # at the threshold: not saturated
    ccd.data[6, 6] = 59999.0
    ccd.data[7, 7] = np.nan      # no comparison is true for NaN
    ccd.data[8, 8] = -5000.0
    expected = np.zeros(ccd.data.shape, dtype=np.uint8)
    for r, c in sat + [(4, 4)]:
        expected[r, c] = 8
    out = run_primitive(FlagSaturation, ccd)
    # source: analytic @ cf88fb8 - docstring: flag value 8 where counts
    # > 60000 (strictly greater), nowhere else
    assert np.array_equal(out.ccddata.flags, expected)
    assert out.ccddata.flags.dtype == np.uint8
    assert out.ccddata.header["SATFLAG"] is True       # source: docstring


@pytest.mark.parametrize("camera,ampmode", [("BLUE", "ALL"), ("RED", "L2U2")])
def test_no_saturation_gives_clear_flags(run_primitive, synthetic_ccddata,
                                         camera, ampmode):
    # source: analytic @ cf88fb8 - a 1000 DN frame has nothing above 60000
    out = run_primitive(FlagSaturation,
                        frame(synthetic_ccddata, camera=camera,
                              ampmode=ampmode))
    assert not out.ccddata.flags.any()
    assert out.ccddata.header["NSATFLAG"] == 0


def test_zero_pixels_are_also_flagged(run_primitive, synthetic_ccddata):
    # source: characterization @ cf88fb8 - the code also adds 8 where
    # data == 0 (a dead/bleed marker); the docstring only mentions > 60000.
    # Pinned so a change is noticed; if it is intended the docstring should
    # say so.
    ccd = frame(synthetic_ccddata)
    ccd.data[2, 3] = 0.0
    out = run_primitive(FlagSaturation, ccd)
    expected = np.zeros(ccd.data.shape, dtype=np.uint8)
    expected[2, 3] = 8
    assert np.array_equal(out.ccddata.flags, expected)


def test_existing_flags_are_replaced(run_primitive, synthetic_ccddata):
    # source: characterization @ cf88fb8 - flags already in the frame are
    # discarded: the primitive stores a fresh array, it does not OR into the
    # old one.
    ccd = frame(synthetic_ccddata)
    ccd.flags[1, 1] = 2
    out = run_primitive(FlagSaturation, ccd)
    assert out.ccddata.flags[1, 1] == 0


@pytest.mark.xfail(strict=True, reason="B28: FlagSaturation.py:41,45 "
                   "NSATFLAG is np.sum of the tuple np.where returns, i.e. "
                   "the sum of the row and column indices, not the pixel "
                   "count")
def test_nsatflag_is_the_pixel_count(run_primitive, synthetic_ccddata):
    ccd = frame(synthetic_ccddata)
    for r, c in [(3, 5), (10, 7), (20, 40)]:
        ccd.data[r, c] = 65535.0
    out = run_primitive(FlagSaturation, ccd)
    # source: analytic @ cf88fb8 - docstring: "the count of saturated
    # pixels": three pixels were saturated.
    assert out.ccddata.header["NSATFLAG"] == 3
