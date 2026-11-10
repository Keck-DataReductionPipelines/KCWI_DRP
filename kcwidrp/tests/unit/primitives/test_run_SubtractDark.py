"""SubtractDark run through ``run_primitive`` (task T.4, detector level).

Specification (class docstring): the master dark is looked up in the proc
table (TYPE MDARK, same camera and CCD configuration; of several, the one
nearest in MJD), scaled by the ratio of exposure times (frame TTIME / dark
TTIME) and subtracted; DARKSUB, MDFILE and DARKSCL record it.  Without a
master dark the frame is left alone and DARKSUB is False.

The master file name is the proc-table file name with ``_mdark`` appended
(``redux/<stem>_mdark.fits``).  Masters are written as MakeMasterDark would
write them: ``kcwi_fits_writer`` to ``redux/`` plus a proc-table row.  Values
are exactly representable in float32 (the writer's storage type).
"""
import numpy as np
import pytest

from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer
from kcwidrp.primitives.SubtractDark import SubtractDark


def master_dark(synthetic_ccddata, tmp_redux, *, stem, value, ttime=600.0,
                mjd=61041.0, frameno=7, binning=(2, 2), drop_ttime=False):
    m = synthetic_ccddata(binning=binning, FRAMENO=frameno,
                          OFNAME=stem + ".fits", IMTYPE="MDARK",
                          TTIME=ttime, MJD=mjd)
    m.data = np.full(m.data.shape, value)
    if drop_ttime:
        del m.header["TTIME"]
    m.header["HISTORY"] = "synthetic master dark"
    name = stem + "_mdark.fits"
    kcwi_fits_writer(m, output_file=name, output_dir=str(tmp_redux.redux))
    if drop_ttime:
        m.header["TTIME"] = ttime      # the proc table needs a value
    tmp_redux.proctab.update_proctab(frame=m, suffix="mdark",
                                     newtype="MDARK", filename=stem + ".fits")
    return name


def test_master_dark_scaled_by_exposure_time_and_subtracted(
        run_primitive, synthetic_ccddata, tmp_redux):
    name = master_dark(synthetic_ccddata, tmp_redux, stem="kb260101_00007",
                       value=8.0, ttime=600.0)
    ccd = synthetic_ccddata(TTIME=300.0)
    ccd.data = np.full(ccd.data.shape, 100.0)
    out = run_primitive(SubtractDark, ccd)
    # source: analytic @ cf88fb8 - scale 300/600 = 0.5; 100 - 8*0.5 = 96
    assert np.array_equal(out.ccddata.data, np.full(ccd.data.shape, 96.0))
    hdr = out.ccddata.header
    assert hdr["DARKSUB"] is True                  # source: docstring
    assert hdr["MDFILE"] == name                   # source: docstring
    assert hdr["DARKSCL"] == 0.5                   # source: analytic


def test_longer_exposure_scales_dark_up(run_primitive, synthetic_ccddata,
                                        tmp_redux):
    master_dark(synthetic_ccddata, tmp_redux, stem="kb260101_00007",
                value=2.0, ttime=100.0)
    ccd = synthetic_ccddata(TTIME=400.0)
    ccd.data = np.full(ccd.data.shape, 100.0)
    out = run_primitive(SubtractDark, ccd)
    # source: analytic @ cf88fb8 - scale 400/100 = 4; 100 - 2*4 = 92
    assert np.array_equal(out.ccddata.data, np.full(ccd.data.shape, 92.0))
    assert out.ccddata.header["DARKSCL"] == 4.0    # source: analytic


def test_nearest_master_in_mjd_is_used(run_primitive, synthetic_ccddata,
                                       tmp_redux):
    # source: analytic @ cf88fb8 - docstring of the proctab search: with
    # several masters the nearest in MJD wins.  Frame MJD 61041.0; masters at
    # 61040.9 (value 1) and 61042.5 (value 3) -> the first.
    near = master_dark(synthetic_ccddata, tmp_redux, stem="kb260101_00007",
                       value=1.0, ttime=100.0, mjd=61040.9, frameno=7)
    master_dark(synthetic_ccddata, tmp_redux, stem="kb260101_00008",
                value=3.0, ttime=100.0, mjd=61042.5, frameno=8)
    ccd = synthetic_ccddata(TTIME=100.0, MJD=61041.0)
    ccd.data = np.full(ccd.data.shape, 10.0)
    out = run_primitive(SubtractDark, ccd)
    assert out.ccddata.header["MDFILE"] == near
    assert np.array_equal(out.ccddata.data, np.full(ccd.data.shape, 9.0))


def test_master_without_exposure_time_is_not_scaled(run_primitive,
                                                    synthetic_ccddata,
                                                    tmp_redux):
    # source: characterization @ cf88fb8 - if either header lacks TTIME the
    # dark is subtracted unscaled (factor 1.0) with a warning.
    master_dark(synthetic_ccddata, tmp_redux, stem="kb260101_00007",
                value=8.0, drop_ttime=True)
    ccd = synthetic_ccddata(TTIME=300.0)
    ccd.data = np.full(ccd.data.shape, 100.0)
    out = run_primitive(SubtractDark, ccd)
    assert np.array_equal(out.ccddata.data, np.full(ccd.data.shape, 92.0))
    assert out.ccddata.header["DARKSCL"] == 1.0


def test_no_master_dark_leaves_frame_alone(run_primitive, synthetic_ccddata):
    # source: characterization @ cf88fb8 - no MDARK row: frame unchanged,
    # DARKSUB False, and (docstring) the reduction continues.
    ccd = synthetic_ccddata()
    before = ccd.data.copy()
    mdfile = ccd.header.get("MDFILE")
    out = run_primitive(SubtractDark, ccd)
    assert np.array_equal(out.ccddata.data, before)
    assert out.ccddata.header["DARKSUB"] is False
    assert out.ccddata.header.get("MDFILE") == mdfile     # not rewritten


def test_master_dark_of_other_binning_not_used(run_primitive,
                                               synthetic_ccddata, tmp_redux):
    # source: analytic @ cf88fb8 - master darks must share the CCD
    # configuration, so a 1x1 master is not applied to a 2x2 frame.
    master_dark(synthetic_ccddata, tmp_redux, stem="kb260101_00007",
                value=8.0, binning=(1, 1))
    ccd = synthetic_ccddata(binning=(2, 2))
    before = ccd.data.copy()
    out = run_primitive(SubtractDark, ccd)
    assert np.array_equal(out.ccddata.data, before)
    assert out.ccddata.header["DARKSUB"] is False
