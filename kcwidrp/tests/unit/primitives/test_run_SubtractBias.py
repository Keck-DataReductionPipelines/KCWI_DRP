"""SubtractBias run through ``run_primitive`` (task T.4, detector level).

Specification (class docstring): the master bias created by MakeMasterBias is
found, subtracted from the frame, the bias read noise BIASRNn is copied into
the frame header, and BIASSUB / MBFILE record the operation.  The master is
found through the proc table: a row of TYPE MBIAS, the frame's camera and
the frame's CCD configuration (CCDCFG); its file is
``redux/mbias_<CCDCFG>.fits`` in the working directory.

The master here is written the way MakeMasterBias writes it
(``kcwi_fits_writer`` to ``redux/`` plus a proc-table row) from synthetic
data.  The data are chosen exactly representable in float32 because the
writer stores float32.
"""
import numpy as np
import pytest

from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer
from kcwidrp.primitives.SubtractBias import SubtractBias


def master_bias(synthetic_ccddata, tmp_redux, *, camera="BLUE",
                ampmode="ALL", binning=(2, 2), value=None, rns=None,
                add_to_proctab=True, **header):
    """Write redux/mbias_<CCDCFG>.fits and register it in the proc table."""
    m = synthetic_ccddata(camera=camera, ampmode=ampmode, binning=binning,
                          FRAMENO=7, OFNAME="mbias_src.fits", IMTYPE="MBIAS",
                          **header)
    shape = m.data.shape
    rr, cc = np.mgrid[0:shape[0], 0:shape[1]]
    m.data = (value if value is not None else 10.0 + 0.5 * cc + rr)
    m.data = np.asarray(m.data, dtype=np.float64) + np.zeros(shape)
    for ia, rn in (rns or {}).items():
        m.header["BIASRN%d" % ia] = rn
    m.header["HISTORY"] = "synthetic master"
    name = "mbias_%s.fits" % m.header["CCDCFG"]
    kcwi_fits_writer(m, output_file=name, output_dir=str(tmp_redux.redux))
    if add_to_proctab:
        tmp_redux.proctab.update_proctab(frame=m, suffix="mbias",
                                         newtype="MBIAS",
                                         filename="kb260101_00007.fits")
    return m, name


def test_master_found_in_proctab_and_subtracted(run_primitive,
                                                synthetic_ccddata,
                                                tmp_redux):
    mbias, name = master_bias(synthetic_ccddata, tmp_redux,
                              rns={1: 2.5, 2: 2.75, 3: 3.5, 4: 3.25})
    ccd = synthetic_ccddata(camera="BLUE", ampmode="ALL", binning=(2, 2))
    ccd.data = np.full(ccd.data.shape, 1000.0)
    out = run_primitive(SubtractBias, ccd)
    rr, cc = np.mgrid[0:ccd.data.shape[0], 0:ccd.data.shape[1]]
    # source: analytic @ cf88fb8 - master = 10 + 0.5 col + row, frame = 1000,
    # so the difference is 990 - 0.5 col - row (all exact in binary floats)
    assert np.array_equal(out.ccddata.data, 990.0 - 0.5 * cc - rr)
    hdr = out.ccddata.header
    assert hdr["BIASSUB"] is True                  # source: docstring
    assert hdr["MBFILE"] == name                   # source: docstring
    # source: docstring - "transfer bias read noise": the master's BIASRNn
    for ia, rn in {1: 2.5, 2: 2.75, 3: 3.5, 4: 3.25}.items():
        assert hdr["BIASRN%d" % ia] == rn


def test_red_master_found_by_red_amp_numbers(run_primitive,
                                             synthetic_ccddata, tmp_redux):
    # source: analytic @ cf88fb8 - RED amps are 0-3, and the BIASRN0..3
    # cards are transferred (not BIASRN1..4).
    master_bias(synthetic_ccddata, tmp_redux, camera="RED",
                ampmode="L2U2L1U1", value=5.0,
                rns={0: 1.5, 1: 2.5, 2: 3.5, 3: 4.5})
    ccd = synthetic_ccddata(camera="RED", ampmode="L2U2L1U1")
    ccd.data = np.full(ccd.data.shape, 105.0)
    out = run_primitive(SubtractBias, ccd)
    assert np.array_equal(out.ccddata.data,
                          np.full(ccd.data.shape, 100.0))  # source: analytic
    for ia, rn in {0: 1.5, 1: 2.5, 2: 3.5, 3: 4.5}.items():
        assert out.ccddata.header["BIASRN%d" % ia] == rn   # source: docstring


def test_no_master_in_proctab_leaves_frame_alone(run_primitive,
                                                 synthetic_ccddata):
    # source: characterization @ cf88fb8 - with no MBIAS row the primitive
    # does not fail: the frame is unchanged and BIASSUB is False.
    ccd = synthetic_ccddata()
    before = ccd.data.copy()
    mbfile = ccd.header.get("MBFILE")
    out = run_primitive(SubtractBias, ccd)
    assert np.array_equal(out.ccddata.data, before)
    assert out.ccddata.header["BIASSUB"] is False
    assert out.ccddata.header.get("MBFILE") == mbfile   # not rewritten


def test_master_with_other_ccd_configuration_is_not_used(run_primitive,
                                                         synthetic_ccddata,
                                                         tmp_redux):
    # source: analytic @ cf88fb8 - docstring of search: biases must share the
    # CCD configuration.  A 1x1 master is not applied to a 2x2 frame.
    master_bias(synthetic_ccddata, tmp_redux, binning=(1, 1), value=5.0)
    ccd = synthetic_ccddata(binning=(2, 2))
    before = ccd.data.copy()
    out = run_primitive(SubtractBias, ccd)
    assert np.array_equal(out.ccddata.data, before)
    assert out.ccddata.header["BIASSUB"] is False


@pytest.mark.xfail(strict=True, reason="B17: a proc table whose only rows "
                   "are the other camera's gives IndexError in "
                   "search_proctab (kcwi_proctab.py:154) instead of "
                   "'no master found'")
def test_other_camera_master_is_not_used(run_primitive, synthetic_ccddata,
                                         tmp_redux):
    # source: analytic @ cf88fb8 - a RED master never matches a BLUE frame.
    master_bias(synthetic_ccddata, tmp_redux, camera="RED",
                ampmode="L2U2L1U1", value=5.0)
    ccd = synthetic_ccddata(camera="BLUE")
    before = ccd.data.copy()
    out = run_primitive(SubtractBias, ccd)
    assert np.array_equal(out.ccddata.data, before)
    assert out.ccddata.header["BIASSUB"] is False


def test_proctab_row_without_file_raises(run_primitive, synthetic_ccddata,
                                         tmp_redux):
    # source: characterization @ cf88fb8 - the proc table lists a master but
    # redux/ has no file: kcwi_fits_reader raises FileNotFoundError, which
    # propagates (the framework would log it and stop the frame).
    _, name = master_bias(synthetic_ccddata, tmp_redux)
    (tmp_redux.redux / name).unlink()
    with pytest.raises(FileNotFoundError):
        run_primitive(SubtractBias, synthetic_ccddata())


@pytest.mark.xfail(strict=True, reason="B17: search_proctab filters by type "
                   "with a mask built from the full table (kcwi_proctab.py:"
                   "154), so a proc table holding both cameras raises "
                   "IndexError")
def test_master_found_when_proctab_holds_both_cameras(run_primitive,
                                                      synthetic_ccddata,
                                                      tmp_redux):
    master_bias(synthetic_ccddata, tmp_redux, camera="BLUE", value=5.0)
    master_bias(synthetic_ccddata, tmp_redux, camera="RED",
                ampmode="L2U2L1U1", value=7.0)
    ccd = synthetic_ccddata(camera="BLUE")
    ccd.data = np.full(ccd.data.shape, 105.0)
    out = run_primitive(SubtractBias, ccd)
    # source: analytic @ cf88fb8 - the BLUE master (5) is subtracted
    assert np.array_equal(out.ccddata.data, np.full(ccd.data.shape, 100.0))
