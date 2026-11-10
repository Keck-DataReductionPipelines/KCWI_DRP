"""SubtractSky run through ``run_primitive`` (task T.8, suite quality check).

Written in a fresh session from the class docstring and the plan, after a
deliberate breakage (the sky scale inverted) went unnoticed by the suite.

Specification (class docstring): the master sky made by MakeMasterSky is
subtracted from the observation.  Which model is used is decided by the
``kcwi.sky`` file if it names one for this frame, otherwise it is the
proc-table row of TYPE SKY with the frame's STATEID (the nearest in MJD).
If the exposure times of the observation and the model differ, the model is
scaled before it is subtracted, and the scale value is recorded in SKYSCL
(physically: sky counts grow in proportion to the exposure time, so the
model is multiplied by t_obs / t_sky).  The master used is recorded in
SKYMAST.  The MakeCube docstring lists a ``noskysub`` extension: the frame
as it was before the subtraction.

The master is written the way MakeMasterSky writes it (``kcwi_fits_writer``
to ``redux/<stem>_sky.fits`` plus a proc-table row of type SKY).  The data
are exact binary fractions, so float comparisons are exact (the writer
stores float32).
"""
import numpy as np
import pytest
from astropy.io import fits

from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer
from kcwidrp.primitives.SubtractSky import SubtractSky

HEAD = "35e41e7"


def master_sky(synthetic_ccddata, tmp_redux, *, stem="kb260101_00007",
               xposure=100.0, mjd=61041.0, level=5.0):
    """Write redux/<stem>_sky.fits and register it as TYPE SKY."""
    m = synthetic_ccddata(FRAMENO=7, OFNAME=stem + ".fits", XPOSURE=xposure,
                          MJD=mjd)
    cc = np.mgrid[0:m.data.shape[0], 0:m.data.shape[1]][1]
    m.data = level + 0.25 * cc          # exact in float32
    m.header["HISTORY"] = "synthetic master sky"
    name = stem + "_sky.fits"
    kcwi_fits_writer(m, output_file=name, output_dir=str(tmp_redux.redux))
    tmp_redux.proctab.update_proctab(frame=m, suffix="sky", newtype="SKY",
                                     filename=stem + ".fits")
    return m, name


def science(synthetic_ccddata, xposure, value=1000.0):
    ccd = synthetic_ccddata(FRAMENO=1, OFNAME="kb260101_00001.fits",
                            XPOSURE=xposure)
    ccd.data = np.full(ccd.data.shape, value)
    return ccd


@pytest.mark.parametrize("t_obs, t_sky, scale", [
    (300.0, 100.0, 3.0),     # science longer than the sky model: scale up
    (100.0, 400.0, 0.25),    # science shorter: scale down
    (100.0, 100.0, 1.0),
])
def test_sky_scaled_by_ratio_of_exposure_times_and_subtracted(
        run_primitive, synthetic_ccddata, tmp_redux, t_obs, t_sky, scale):
    _, name = master_sky(synthetic_ccddata, tmp_redux, xposure=t_sky)
    ccd = science(synthetic_ccddata, t_obs)
    cc = np.mgrid[0:ccd.data.shape[0], 0:ccd.data.shape[1]][1]
    out = run_primitive(SubtractSky, ccd).ccddata
    # source: analytic - out = obs - (t_obs / t_sky) * sky, with
    # sky = 5 + 0.25 col; the scales and values are exact binary fractions
    assert np.array_equal(out.data, 1000.0 - scale * (5.0 + 0.25 * cc))
    hdr = out.header
    assert hdr["SKYCOR"] is True                     # source: docstring
    assert hdr["SKYSCL"] == scale                    # source: docstring
    assert hdr["SKYMAST"] == name                    # source: docstring


def test_frame_before_subtraction_is_kept_as_noskysub(run_primitive,
                                                      synthetic_ccddata,
                                                      tmp_redux):
    # source: characterization @ 35e41e7 (MakeCube's docstring lists a
    # noskysub extension): the frame as it was before the subtraction
    master_sky(synthetic_ccddata, tmp_redux)
    ccd = science(synthetic_ccddata, 100.0)
    out = run_primitive(SubtractSky, ccd).ccddata
    assert np.array_equal(out.noskysub, np.full(ccd.data.shape, 1000.0))


def test_nearest_master_in_mjd_is_used(run_primitive, synthetic_ccddata,
                                       tmp_redux):
    # source: analytic - docstring: nearest model in MJD of the same state
    far, _ = master_sky(synthetic_ccddata, tmp_redux, stem="kb260101_00007",
                        mjd=61040.0, level=50.0)
    near, near_name = master_sky(synthetic_ccddata, tmp_redux,
                                 stem="kb260101_00008", mjd=61041.1,
                                 level=5.0)
    ccd = science(synthetic_ccddata, 100.0)         # MJD 61041.2083
    out = run_primitive(SubtractSky, ccd).ccddata
    assert out.header["SKYMAST"] == near_name


@pytest.mark.parametrize("bad", [0.0, -1.0])
def test_invalid_exposure_time_gives_unit_scale(run_primitive,
                                                synthetic_ccddata,
                                                tmp_redux, bad):
    # source: characterization @ 35e41e7 - with an exposure time <= 0 the
    # code logs a warning and subtracts the unscaled model
    master_sky(synthetic_ccddata, tmp_redux, xposure=100.0)
    ccd = science(synthetic_ccddata, bad)
    out = run_primitive(SubtractSky, ccd).ccddata
    assert out.header["SKYSCL"] == 1.0


def test_kcwi_sky_file_selects_the_model(run_primitive, synthetic_ccddata,
                                         tmp_redux):
    # source: analytic - docstring: the kcwi.sky file, if it names a sky for
    # this science file, decides which model is used ("<raw sci file> <raw
    # sky file>"; the model is <sky stem>_sky.fits)
    master_sky(synthetic_ccddata, tmp_redux, stem="kb260101_00007",
               level=50.0)
    _, chosen = master_sky(synthetic_ccddata, tmp_redux,
                           stem="kb260101_00008", mjd=61000.0, level=5.0)
    (tmp_redux.root / "kb260101_00008.fits").write_bytes(b"")  # must exist
    ccd = science(synthetic_ccddata, 100.0)
    name = ccd.header["OFNAME"]
    (tmp_redux.root / "kcwi.sky").write_text(
        "# science sky\n%s kb260101_00008.fits\n" % name)
    out = run_primitive(SubtractSky, ccd)
    assert out.ccddata.header["SKYMAST"] == chosen


def test_skipsky_leaves_frame_alone(run_primitive, synthetic_ccddata,
                                    tmp_redux):
    # source: analytic - docstring/config: skipsky = True turns sky
    # subtraction off; SKYCOR records that nothing was subtracted
    master_sky(synthetic_ccddata, tmp_redux)
    ccd = science(synthetic_ccddata, 100.0)
    out = run_primitive(SubtractSky, ccd,
                        config_overrides={"skipsky": True})
    # source: characterization @ 35e41e7 - a failed precondition makes
    # BasePrimitive.apply return None and the frame is left as it was
    assert out is None
    assert np.array_equal(ccd.data, np.full(ccd.data.shape, 1000.0))
    assert ccd.header["SKYCOR"] is False


def test_no_master_leaves_frame_alone(run_primitive, synthetic_ccddata):
    # source: analytic - precondition: with no sky model to subtract the
    # frame is passed on unchanged
    ccd = science(synthetic_ccddata, 100.0)
    before = ccd.header.copy()
    out = run_primitive(SubtractSky, ccd)
    # source: characterization @ 35e41e7 - failed precondition: apply
    # returns None and the frame is left as it was
    assert out is None
    assert np.array_equal(ccd.data, np.full(ccd.data.shape, 1000.0))
    assert ccd.header.get("SKYSCL") == before.get("SKYSCL")
    assert ccd.header.get("SKYMAST") == before.get("SKYMAST")


def test_output_intk_written_and_in_proctab(run_primitive, synthetic_ccddata,
                                            tmp_redux):
    # source: docstring - "Writes out a *_intk.fits file and adds an entry
    # to the proc file"
    master_sky(synthetic_ccddata, tmp_redux)
    ccd = science(synthetic_ccddata, 100.0)
    run_primitive(SubtractSky, ccd)
    path = tmp_redux.redux / "kb260101_00001_intk.fits"
    assert path.exists()
    with fits.open(path) as hdul:
        assert hdul[0].header["SKYCOR"] is True
    assert "intk" in [str(s) for s in tmp_redux.proctab.proctab["SUFF"]]
