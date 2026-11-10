"""CorrectGain run through ``run_primitive`` (task T.4, detector level).

Specification (class docstring): the image is divided into amplifier regions
by ATSECn and each region is multiplied by its GAINn; GAINCOR becomes True and
BUNIT 'electron'.  BLUE amplifiers are numbered 1-4, RED 0-3, and the gain is
looked up by that header number.

The input is a trimmed image (as TrimOverscan leaves it) built here: a
constant per amplifier region, so the expected output is constant * gain.
Regions are laid out as in the fixture's amplifier geometry; ATSECn is
written into the header by this test, independent of TrimOverscan.
"""
import numpy as np
import pytest

from kcwidrp.primitives.CorrectGain import CorrectGain

# camera, ampmode -> {amp number: ATSEC string} for a 6 row x 8 column image
# (4-amp: 3 rows per half, 4 columns per half)
ATSEC = {
    ("BLUE", "ALL"): {1: "[1:4,1:3]", 2: "[5:8,1:3]", 3: "[1:4,4:6]",
                      4: "[5:8,4:6]"},
    ("BLUE", "TUP"): {1: "[1:4,1:6]", 2: "[5:8,1:6]"},
    ("RED", "L2U2L1U1"): {0: "[1:4,4:6]", 1: "[1:4,1:3]", 2: "[5:8,4:6]",
                          3: "[5:8,1:3]"},
    ("RED", "L2U2"): {1: "[1:4,1:6]", 3: "[5:8,1:6]"},
}
CASES = list(ATSEC)
IDS = ["%s-%s" % c for c in CASES]


def region(sec):
    x, y = sec.strip("[]").split(",")
    x0, x1 = map(int, x.split(":"))
    y0, y1 = map(int, y.split(":"))
    return slice(y0 - 1, y1), slice(x0 - 1, x1)


@pytest.mark.parametrize("camera,ampmode", CASES, ids=IDS)
def test_each_amp_region_scaled_by_its_header_gain(run_primitive,
                                                   synthetic_ccddata, camera,
                                                   ampmode):
    ccd = synthetic_ccddata(camera=camera, ampmode=ampmode)
    amps = list(ATSEC[(camera, ampmode)])
    # distinct, exactly representable gains (multiples of 1/8), not the
    # fixture's, so the test cannot pass by using a wrong amp's gain
    gains = {ia: 1.0 + (k + 1) / 8.0 for k, ia in enumerate(amps)}
    data = np.zeros((6, 8))
    expected = np.zeros((6, 8))
    for ia in amps:
        ccd.header["ATSEC%d" % ia] = ATSEC[(camera, ampmode)][ia]
        ccd.header["GAIN%d" % ia] = gains[ia]
        rows, cols = region(ATSEC[(camera, ampmode)][ia])
        data[rows, cols] = 100.0 * (ia + 1)
        # source: analytic @ cf88fb8 - docstring: region * GAINn
        expected[rows, cols] = 100.0 * (ia + 1) * gains[ia]
    ccd.data = data
    ccd.flags = np.zeros((6, 8), dtype=np.uint8)
    out = run_primitive(CorrectGain, ccd)
    assert np.array_equal(out.ccddata.data, expected)
    hdr = out.ccddata.header
    assert hdr["GAINCOR"] is True               # source: docstring
    assert hdr["BUNIT"] == "electron"           # source: docstring
    assert out.ccddata.unit == "electron"       # source: docstring


def test_gain_by_amp_number_blue_is_one_based(run_primitive,
                                              synthetic_ccddata):
    # source: independent @ cf88fb8 - BLUE has amps 1-4: GAIN0 exists in no
    # BLUE header, and a GAIN0 card must be ignored.
    ccd = synthetic_ccddata(camera="BLUE", ampmode="ALL")
    ccd.header["GAIN0"] = 99.0
    for ia, sec in ATSEC[("BLUE", "ALL")].items():
        ccd.header["ATSEC%d" % ia] = sec
        ccd.header["GAIN%d" % ia] = 2.0
    ccd.data = np.ones((6, 8))
    out = run_primitive(CorrectGain, ccd)
    assert np.array_equal(out.ccddata.data, np.full((6, 8), 2.0))


def test_gain_by_amp_number_red_is_zero_based(run_primitive,
                                              synthetic_ccddata):
    # source: independent @ cf88fb8 - RED has amps 0-3 and GAIN0 is used; a
    # GAIN4 card (nonexistent amp) must be ignored.
    ccd = synthetic_ccddata(camera="RED", ampmode="L2U2L1U1")
    ccd.header["GAIN4"] = 99.0
    for ia, sec in ATSEC[("RED", "L2U2L1U1")].items():
        ccd.header["ATSEC%d" % ia] = sec
        ccd.header["GAIN%d" % ia] = 0.5 * (ia + 1)
    ccd.data = np.ones((6, 8))
    out = run_primitive(CorrectGain, ccd)
    expected = np.ones((6, 8))
    for ia, sec in ATSEC[("RED", "L2U2L1U1")].items():
        rows, cols = region(sec)
        expected[rows, cols] = 0.5 * (ia + 1)    # source: analytic
    assert np.array_equal(out.ccddata.data, expected)


def test_gain_applied_to_trimmed_output_of_fixture_frame(run_primitive,
                                                         synthetic_ccddata):
    # source: independent @ cf88fb8 - chained after TrimOverscan, the header's
    # own GAINn are applied.  The fixture's BLUE ALL gains are the post-ingest
    # values of real frame kb260115_00011 (amp 1-4: 0.157, 0.160, 0.158,
    # 0.158), so a flat 1000 DN image becomes 1000 * gain in each quadrant.
    from kcwidrp.primitives.TrimOverscan import TrimOverscan
    ccd = synthetic_ccddata(camera="BLUE", ampmode="ALL")
    ccd.data = np.full(ccd.data.shape, 1000.0)
    trimmed = run_primitive(TrimOverscan, ccd)
    out = run_primitive(CorrectGain, trimmed)
    gains = {1: 0.157, 2: 0.160, 3: 0.158, 4: 0.158}
    d = out.ccddata.data
    assert d.dtype == np.float32
    # Tolerance: data are float32 (TrimOverscan), so 1000 * gain is rounded
    # to float32 relative precision 6e-8; rtol 1e-6 leaves a margin of ~16.
    np.testing.assert_allclose(d[:64, :32], 1000.0 * gains[1], rtol=1e-6)
    np.testing.assert_allclose(d[:64, 32:], 1000.0 * gains[2], rtol=1e-6)
    np.testing.assert_allclose(d[64:, :32], 1000.0 * gains[3], rtol=1e-6)
    np.testing.assert_allclose(d[64:, 32:], 1000.0 * gains[4], rtol=1e-6)
