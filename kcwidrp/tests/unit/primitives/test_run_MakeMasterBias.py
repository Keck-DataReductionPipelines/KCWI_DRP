"""MakeMasterBias run through ``run_primitive`` (task T.4, detector level).

Specification (class docstring): biases that went through TrimOverscan
(``redux/<stem>_intb.fits`` plus a BIAS row in the proc table) are stacked
once ``min_files`` of them are in the proc table: method 'average' with an
upper 2-sigma clip (centre = median, deviation = MAD-based sigma) through
``ccdproc.combine``.  The master is written to ``redux/mbias_<CCDCFG>.fits``
and entered in the proc table as MBIAS.  The read noise of each amplifier,
BIASRNn, comes from the difference of the 2nd and 3rd bias.

The inputs are built the way the pipeline builds them: synthetic BIAS frames
run through TrimOverscan (which writes the ``_intb`` files and proc-table
rows), whose data are constants per frame, so averages and spreads follow
from arithmetic.
"""
import numpy as np
import pytest
from astropy.io import fits

from kcwidrp.primitives.MakeMasterBias import MakeMasterBias
from kcwidrp.primitives.TrimOverscan import TrimOverscan

LEVELS = [98.0, 99.0, 100.0, 101.0, 102.0]


def make_biases(run_primitive, synthetic_ccddata, levels, patterns=None,
                hot=None):
    """Run TrimOverscan on one BIAS frame per level; return the last output.

    ``patterns[i]`` (+1 / -1 checkerboard amplitude d) is added to frame i
    as d * (-1)**(row+col); ``hot`` = (frame index, row, col, value) sets one
    pixel of one frame.
    """
    out = None
    for i, level in enumerate(levels):
        ccd = synthetic_ccddata(camera="BLUE", ampmode="ALL", binning=(2, 2),
                                IMTYPE="BIAS", FRAMENO=i + 1,
                                OFNAME="kb260101_%05d.fits" % (i + 1),
                                TTIME=0.0, EXPTIME=0.0)
        ny, nx = ccd.data.shape
        rr, cc = np.mgrid[0:ny, 0:nx]
        ccd.data = np.full((ny, nx), level)
        # A raw frame has no uncertainty or mask yet (CreateUncertaintyImage
        # makes them after trimming); the fixture's full-size ones would not
        # match the trimmed data when the _intb file is read back.
        ccd.uncertainty = None
        ccd.mask = None
        if patterns and patterns.get(i):
            ccd.data = ccd.data + patterns[i] * (-1.0) ** (rr + cc)
        if hot is not None and hot[0] == i:
            ccd.data[hot[1], hot[2]] = hot[3]
        out = run_primitive(TrimOverscan, ccd)
    return out


def stack(run_primitive, last, min_files=5):
    return run_primitive(MakeMasterBias, last,
                         action_extra={"min_files": min_files,
                                       "new_type": "MBIAS"})


def read_master(tmp_redux, out):
    path = tmp_redux.redux / out.name
    with fits.open(path) as hdul:
        return (hdul[0].data.copy(), hdul["UNCERT"].data.copy(),
                hdul[0].header.copy())


def test_average_of_five_biases_and_its_uncertainty(run_primitive,
                                                    synthetic_ccddata,
                                                    tmp_redux):
    last = make_biases(run_primitive, synthetic_ccddata, LEVELS)
    out = stack(run_primitive, last)
    data, unc, hdr = read_master(tmp_redux, out)
    # source: analytic @ cf88fb8 - mean of 98..102 is 100; nothing is clipped
    # (102 is 2 above the median 100, the clip limit is 2 * 1.4826 * MAD(1)
    # = 2.97 above)
    # Tolerance: written as float32: relative 6e-8.
    np.testing.assert_allclose(data, 100.0, rtol=1e-6)
    # source: analytic @ cf88fb8 - ccdproc's average uncertainty is the
    # standard error: population std / sqrt(N).  Population variance of
    # 98..102 is (4+1+0+1+4)/5 = 2, so sigma = sqrt(2) / sqrt(5) = sqrt(0.4).
    # Tolerance: float32 storage, relative 6e-8.
    np.testing.assert_allclose(unc, np.sqrt(0.4), rtol=1e-6)
    assert data.shape == (128, 64)                      # source: analytic
    # source: analytic @ cf88fb8 - docstring / header keywords
    assert hdr["NSTACK"] == 5
    assert hdr["STCKMETH"] == "average"
    assert hdr["STCKSIGU"] == 2.0
    assert hdr["IMTYPE"] == "MBIAS"
    for i in range(5):
        assert hdr["STACKF%d" % (i + 1)] == "kb260101_%05d_intb.fits" % (
            i + 1)


def test_master_name_and_proctab_entry(run_primitive, synthetic_ccddata,
                                       tmp_redux):
    last = make_biases(run_primitive, synthetic_ccddata, LEVELS)
    ccdcfg = last.ccddata.header["CCDCFG"]
    out = stack(run_primitive, last)
    # source: analytic @ cf88fb8 - docstring: written as mbias_<CCDCFG>.fits
    # (the name SubtractBias looks for)
    assert out.name == "mbias_%s.fits" % ccdcfg
    assert (tmp_redux.redux / out.name).is_file()
    tab = tmp_redux.proctab.proctab
    mrows = tab[tab["TYPE"] == "MBIAS"]
    assert len(mrows) == 1                     # source: docstring
    assert mrows["SUFF"][0] == "mbias"
    assert tmp_redux.proctab_file.is_file()    # proc table rewritten


def test_too_few_biases_does_nothing(run_primitive, synthetic_ccddata,
                                     tmp_redux):
    # source: analytic @ cf88fb8 - docstring: stacks "based on the config
    # parameter min_files"; with 4 biases and min_files 5 the precondition
    # fails, apply() returns None and no master is written.
    last = make_biases(run_primitive, synthetic_ccddata, LEVELS[:4])
    out = stack(run_primitive, last, min_files=5)
    assert out is None
    assert not list(tmp_redux.redux.glob("mbias_*"))
    assert not (tmp_redux.proctab.proctab["TYPE"] == "MBIAS").any()


def test_exact_minimum_number_is_enough(run_primitive, synthetic_ccddata,
                                        tmp_redux):
    # source: analytic @ cf88fb8 - ">= min_files": 3 biases, min_files 3
    last = make_biases(run_primitive, synthetic_ccddata, [99.0, 100.0, 101.0])
    out = stack(run_primitive, last, min_files=3)
    data, _, hdr = read_master(tmp_redux, out)
    np.testing.assert_allclose(data, 100.0, rtol=1e-6)  # source: analytic
    assert hdr["NSTACK"] == 3


def test_hot_pixel_in_one_bias_is_rejected(run_primitive, synthetic_ccddata,
                                           tmp_redux):
    # source: analytic @ cf88fb8 - four biases at 100 and one at 100 with a
    # 1000 DN hot pixel: the median at that pixel is 100 and the MAD is 0,
    # so the upper clip limit is 100 and 1000 is rejected; the average of
    # the remaining four is 100.  (The frames are trimmed first, so the hot
    # pixel is placed in a data region of the raw frame: raw (10, 40) is
    # inside the first amplifier's data section DSEC1 for BLUE ALL 2x2.)
    last = make_biases(run_primitive, synthetic_ccddata, [100.0] * 5,
                       hot=(2, 10, 40, 1000.0))
    out = stack(run_primitive, last)
    data, _, _ = read_master(tmp_redux, out)
    np.testing.assert_allclose(data, 100.0, rtol=1e-6)


def test_read_noise_from_second_and_third_bias(run_primitive,
                                               synthetic_ccddata,
                                               tmp_redux):
    # source: analytic @ cf88fb8 - docstring: BIASRNn from the difference of
    # the 2nd and 3rd bias.  Frames 2 and 3 carry +/- d/2 checkerboards
    # (d = 4 DN), so their difference is +/- d on every pixel: its std is d.
    # Per amp: RN = std(diff) * gain / sqrt(2) [two independent reads add in
    # variance] with the post-ingest gains of the fixture (0.157, 0.160,
    # 0.158, 0.158 e-/DN, from real frame kb260115_00011).
    # Tolerance: the header value is rounded to 3 decimals (+-5e-4), and the
    # code divides by 1.414 rather than sqrt(2) (a 1.5e-4 relative
    # difference, 7e-5 e- here): 1e-3 e- covers both.
    d = 4.0
    last = make_biases(run_primitive, synthetic_ccddata, LEVELS,
                       patterns={1: d / 2.0, 2: -d / 2.0})
    out = stack(run_primitive, last)
    _, _, hdr = read_master(tmp_redux, out)
    for ia, gain in {1: 0.157, 2: 0.160, 3: 0.158, 4: 0.158}.items():
        assert hdr["BIASRN%d" % ia] == pytest.approx(
            d * gain / np.sqrt(2.0), abs=1e-3)
