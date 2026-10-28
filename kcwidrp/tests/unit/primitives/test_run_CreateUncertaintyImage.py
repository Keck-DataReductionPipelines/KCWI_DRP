"""CreateUncertaintyImage run through ``run_primitive`` (task T.4).

Physics (class docstring: "starts with pure Poisson noise ... uses BIASRNn
with ATSECn to apply the appropriate readnoise addition to each amplifier
region"; the input is in electrons).  For a pixel with N electrons in an
amplifier with read noise RN electrons the standard deviation is

    sigma = sqrt(Var_shot + Var_read) = sqrt(N + RN**2)

(shot noise is Poisson with variance N; the two sources add in variance).

To keep the comparison exact the test picks N a perfect square and
N + RN**2 a perfect square, so the right answer is an integer that no
floating-point formula can get wrong by rounding:

    N = 16, RN = 3  ->  sqrt(25)  = 5
    N =  9, RN = 4  ->  sqrt(25)  = 5
    N = 144, RN = 5 ->  sqrt(169) = 13
    N = 64, RN = 6  ->  sqrt(100) = 10
"""
import numpy as np
import pytest
from astropy.nddata import StdDevUncertainty

from kcwidrp.primitives.CreateUncertaintyImage import CreateUncertaintyImage

# (electrons, read noise, sigma) per amp, in header amp order
CASES_AMP = [(16.0, 3.0, 5.0), (9.0, 4.0, 5.0), (144.0, 5.0, 13.0),
             (64.0, 6.0, 10.0)]
# 6 x 8 trimmed image; ATSEC strings per camera/mode (same geometry as
# test_run_CorrectGain)
ATSEC = {
    ("BLUE", "ALL"): {1: "[1:4,1:3]", 2: "[5:8,1:3]", 3: "[1:4,4:6]",
                      4: "[5:8,4:6]"},
    ("BLUE", "TUP"): {1: "[1:4,1:6]", 2: "[5:8,1:6]"},
    ("RED", "L2U2L1U1"): {0: "[1:4,4:6]", 1: "[1:4,1:3]", 2: "[5:8,4:6]",
                          3: "[5:8,1:3]"},
    ("RED", "L2U2"): {1: "[1:4,1:6]", 3: "[5:8,1:6]"},
}
MODES = list(ATSEC)
IDS = ["%s-%s" % m for m in MODES]


def region(sec):
    x, y = sec.strip("[]").split(",")
    x0, x1 = map(int, x.split(":"))
    y0, y1 = map(int, y.split(":"))
    return slice(y0 - 1, y1), slice(x0 - 1, x1)


def prepared(synthetic_ccddata, camera, ampmode, counts_sign=1.0):
    ccd = synthetic_ccddata(camera=camera, ampmode=ampmode)
    data = np.zeros((6, 8))
    sigma = np.zeros((6, 8))
    for k, (ia, sec) in enumerate(ATSEC[(camera, ampmode)].items()):
        n, rn, sg = CASES_AMP[k]
        ccd.header["ATSEC%d" % ia] = sec
        ccd.header["BIASRN%d" % ia] = rn
        rows, cols = region(sec)
        data[rows, cols] = counts_sign * n
        sigma[rows, cols] = sg
    ccd.data = data
    ccd.flags = None
    ccd.mask = None
    return ccd, sigma


@pytest.mark.parametrize("camera,ampmode", MODES, ids=IDS)
def test_uncertainty_is_poisson_plus_readnoise_per_amp(
        run_primitive, synthetic_ccddata, camera, ampmode):
    ccd, sigma = prepared(synthetic_ccddata, camera, ampmode)
    out = run_primitive(CreateUncertaintyImage, ccd)
    unc = out.ccddata.uncertainty
    assert isinstance(unc, StdDevUncertainty)
    # source: analytic @ cf88fb8 - sigma = sqrt(N + RN^2) per pixel, with the
    # (N, RN, sigma) triples in the module docstring; exact (integers).
    assert np.array_equal(unc.array, sigma)
    assert out.ccddata.header["UNCSTD"] is True        # source: docstring


@pytest.mark.parametrize("camera,ampmode", MODES, ids=IDS)
def test_without_readnoise_keywords_poisson_only(run_primitive,
                                                 synthetic_ccddata, camera,
                                                 ampmode):
    ccd, _ = prepared(synthetic_ccddata, camera, ampmode)
    for ia in ATSEC[(camera, ampmode)]:
        del ccd.header["BIASRN%d" % ia]
    out = run_primitive(CreateUncertaintyImage, ccd)
    # source: analytic @ cf88fb8 - docstring: "if BIASRNn keywords present":
    # with none, sigma = sqrt(N): 4, 3, 12, 8 for N = 16, 9, 144, 64.
    expected = np.zeros((6, 8))
    for k, (ia, sec) in enumerate(ATSEC[(camera, ampmode)].items()):
        expected[region(sec)] = np.sqrt(CASES_AMP[k][0])
    assert np.array_equal(out.ccddata.uncertainty.array, expected)


def test_negative_counts_use_absolute_value(run_primitive,
                                            synthetic_ccddata):
    # source: characterization @ cf88fb8 - for N < 0 (bias-subtracted noise)
    # the shot-noise term uses |N|, so -16 e- with RN 3 gives 5, like +16.
    ccd, sigma = prepared(synthetic_ccddata, "BLUE", "ALL", counts_sign=-1.0)
    out = run_primitive(CreateUncertaintyImage, ccd)
    assert np.array_equal(out.ccddata.uncertainty.array, sigma)


def test_nan_pixel_gives_nan_uncertainty(run_primitive, synthetic_ccddata):
    # source: analytic @ cf88fb8 - sqrt(NaN + RN^2) is NaN, and only there.
    ccd, sigma = prepared(synthetic_ccddata, "BLUE", "ALL")
    ccd.data[1, 1] = np.nan
    sigma[1, 1] = np.nan
    out = run_primitive(CreateUncertaintyImage, ccd)
    assert np.array_equal(out.ccddata.uncertainty.array, sigma,
                          equal_nan=True)


def test_missing_flags_and_mask_are_created(run_primitive, synthetic_ccddata):
    # source: characterization @ cf88fb8 - absent flags and mask become
    # uint8 zero arrays of the data shape.
    ccd, _ = prepared(synthetic_ccddata, "BLUE", "ALL")
    out = run_primitive(CreateUncertaintyImage, ccd)
    assert out.ccddata.flags.dtype == np.uint8
    assert out.ccddata.flags.shape == (6, 8)
    assert not out.ccddata.flags.any()
    assert out.ccddata.mask is not None
    assert not np.asarray(out.ccddata.mask).any()
    assert np.asarray(out.ccddata.mask).shape == (6, 8)


def test_input_data_untouched(run_primitive, synthetic_ccddata):
    # source: analytic @ cf88fb8 - an uncertainty image does not change data
    ccd, _ = prepared(synthetic_ccddata, "BLUE", "ALL")
    before = ccd.data.copy()
    out = run_primitive(CreateUncertaintyImage, ccd)
    assert np.array_equal(out.ccddata.data, before)
