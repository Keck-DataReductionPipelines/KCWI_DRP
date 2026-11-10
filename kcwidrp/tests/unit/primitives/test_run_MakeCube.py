"""MakeCube run through ``run_primitive`` (task T.8, suite quality check).

Written in a fresh session from the class docstring, after a deliberate
breakage (the 180-degree rotation of RED cubes removed) went unnoticed by
the suite.

Specification (class docstring):

* The 2D image is cut into 24 slices (columns ``xl0[i]:xl1[i]`` of the
  image) and each slice is warped with its transform ``tform[i]`` into an
  (ysize, xsize) image; the cube is assembled with the slice as the last
  axis, so ``cube[:, :, i]`` is slice i, for the intensity, uncertainty,
  mask and flag planes.  (Axis 0 is the image row, i.e. wavelength; axis 1
  the position along the slice; axis 2 the slice.)  Here every transform is
  the identity, so each warped slice must equal the cut-out of the image.
* "Rotates RED channel data by 180 degrees to align with BLUE channel
  data": for a RED frame the (along-slice, slice) plane of every cube is
  rotated by 180 degrees, ``out[w, i, j] = in[w, xsize-1-i, 23-j]``;
  BLUE data are not rotated.
* "Creates WCS keywords that reflect the new geometry of the cube": CRVAL3
  and CD3_3 from the geometry (``wave0out``, ``dwout``), CRPIX3 = 1, the
  celestial axes RA---TAN / DEC--TAN with CRVAL1/2 the pointing, and the CD
  matrix of a plane with pixel scales ``cdelt1 = -slscl`` (across slices)
  and ``cdelt2 = pxscl`` (along slices) rotated by ``crota``: [CD1_1 CD1_2;
  CD2_1 CD2_2] = [cdelt1 cos, -cdelt2 sin; cdelt1 sin, cdelt2 cos]
  (Calabretta & Greisen 2002, eq. 188).

The multiprocessing pool is replaced by a serial stand-in (a real pool
would start worker processes that re-import the package in every test);
``make_cube_helper`` itself runs unchanged.  The geometry file is a pickle
as SolveGeom writes it, with the keys MakeCube reads.
"""
import math
import pickle

import numpy as np
import pytest
from astropy.io import fits
from skimage import transform as tf

from kcwidrp.primitives import MakeCube as make_cube_module
from kcwidrp.primitives.MakeCube import MakeCube

HEAD = "35e41e7"
NSL, NX, NY = 24, 4, 12            # slices, pixels along a slice, rows
SLSCL, PXSCL = 0.00037718, 8.096e-05      # deg per pixel (Small slicer, 2x2)
WAVE0, DW = 3500.0, 1.25
ARC = "kb260101_00002"
OBJ = "kb260101_00001"


class SerialPool:
    """Stand-in for multiprocessing.Pool: same map/close interface."""

    def map(self, func, iterable):
        return [func(x) for x in iterable]

    def close(self):
        pass


@pytest.fixture(autouse=True)
def serial_pool(monkeypatch):
    monkeypatch.setattr(make_cube_module, "Pool", SerialPool)


def geometry():
    return dict(
        tform=[tf.AffineTransform() for _ in range(NSL)],   # identity
        xl0=[NX * i for i in range(NSL)], xl1=[NX * (i + 1)
                                                for i in range(NSL)],
        xsize=NX, ysize=NY, slscl=SLSCL, pxscl=PXSCL, barsep=1.0, bar0=2.0,
        waveall0=3400.0, waveall1=3700.0, wavegood0=3450.0,
        wavegood1=3650.0, wavemid=3550.0, wavensall0=3400.0,
        wavensall1=3700.0, wavensgood0=3450.0, wavensgood1=3650.0,
        wavensmid=3550.0, avwvsig=0.1, sdwvsig=0.01, cbarsno=3,
        cbarsfl="kb260101_00003.fits", arcno=2, arcfl=ARC + ".fits",
        wave0out=WAVE0, dwout=DW, dich_frac=1.0)


def planes():
    """Distinct values everywhere, exact in float32."""
    rr, cc = np.mgrid[0:NY, 0:NSL * NX]
    data = 100.0 * rr + cc                       # up to 1195
    unc = 0.5 + 0.25 * rr + 0.125 * (cc % 7)
    mask = ((rr * 7 + cc * 3) % 5 == 0)          # bool pattern
    flags = ((rr * 5 + cc) % 4 * 32).astype(np.uint8)   # 0, 32, 64, 96
    return data, unc, mask, flags


def make_frame(synthetic_ccddata, tmp_redux, camera, **header):
    ampmode = None
    ccd = synthetic_ccddata(
        camera=camera, ampmode=ampmode, IMTYPE="OBJECT", FRAMENO=1,
        OFNAME=OBJ + ".fits", RA="10:00:00.0", DEC="+02:00:00.0",
        ROTPOSN=0.0, ROTREFAN=0.0, **header)
    data, unc, mask, flags = planes()
    ccd.data = data.copy()
    ccd.uncertainty.array = unc.copy()
    ccd.mask = mask.copy()
    ccd.flags = flags.copy()
    # the master arc whose geometry is used (TYPE MARC, same state)
    arc = synthetic_ccddata(camera=camera, IMTYPE="MARC", FRAMENO=2,
                            OFNAME=ARC + ".fits")
    tmp_redux.proctab.update_proctab(frame=arc, suffix="marc",
                                     newtype="MARC", filename=ARC + ".fits")
    with open(tmp_redux.redux / (ARC + "_geom.pkl"), "wb") as f:
        pickle.dump(geometry(), f)
    return ccd


def run_cube(run_primitive, ccd, tmp_redux, **overrides):
    out = run_primitive(MakeCube, ccd,
                        config_overrides=dict({"ROTOFF": 0.0}, **overrides))
    with fits.open(tmp_redux.redux / (OBJ + "_icube.fits")) as hdul:
        cube = {n: hdul[n].data.copy() for n in ("PRIMARY", "UNCERT", "MASK",
                                                 "FLAGS")}
        hdr = hdul[0].header.copy()
    return out, cube, hdr


def expected_blue_cube():
    data, unc, mask, flags = planes()
    out = {}
    for name, plane in (("PRIMARY", data), ("UNCERT", unc),
                        ("MASK", mask.astype(np.uint8)), ("FLAGS", flags)):
        cube = np.zeros((NY, NX, NSL), dtype=np.float64)
        for i in range(NSL):
            cube[:, :, i] = plane[:, NX * i:NX * (i + 1)]
        out[name] = cube
    return out


def test_blue_slices_are_placed_unrotated(run_primitive, synthetic_ccddata,
                                          tmp_redux):
    ccd = make_frame(synthetic_ccddata, tmp_redux, "BLUE")
    _, cube, _ = run_cube(run_primitive, ccd, tmp_redux)
    # source: analytic - identity transforms: cube[:, :, i] is the cut-out
    # columns of slice i.  Tolerance: a cubic-spline warp at integer
    # coordinates reproduces the pixel to ~1e-12 relative (values <= 1195);
    # 1e-6 absolute also covers the writer's float32 storage (all inputs
    # are exact in float32, so only the spline round-off remains).
    expected = expected_blue_cube()
    np.testing.assert_allclose(cube["PRIMARY"], expected["PRIMARY"],
                               rtol=0, atol=1e-6)
    np.testing.assert_allclose(cube["UNCERT"], expected["UNCERT"],
                               rtol=0, atol=1e-6)
    # mask and flags: order-0 warp, exact
    assert np.array_equal(cube["MASK"], expected["MASK"])
    assert np.array_equal(cube["FLAGS"], expected["FLAGS"])


def test_red_cube_is_rotated_by_180_degrees_in_slice_plane(
        run_primitive, synthetic_ccddata, tmp_redux):
    ccd = make_frame(synthetic_ccddata, tmp_redux, "RED")
    _, cube, _ = run_cube(run_primitive, ccd, tmp_redux)
    expected = expected_blue_cube()
    for name in ("PRIMARY", "UNCERT", "MASK", "FLAGS"):
        # source: analytic - 180 degree rotation of the (along-slice, slice)
        # plane: out[w, i, j] = in[w, NX-1-i, NSL-1-j]; same tolerance as
        # the BLUE test for the float planes, exact for mask and flags
        rotated = expected[name][:, ::-1, ::-1]
        if name in ("PRIMARY", "UNCERT"):
            np.testing.assert_allclose(cube[name], rotated, rtol=0,
                                       atol=1e-6)
        else:
            assert np.array_equal(cube[name], rotated)


def test_wavelength_axis_and_geometry_keywords(run_primitive,
                                               synthetic_ccddata,
                                               tmp_redux):
    ccd = make_frame(synthetic_ccddata, tmp_redux, "BLUE")
    out, _, hdr = run_cube(run_primitive, ccd, tmp_redux)
    # source: docstring/analytic - values come straight from the geometry
    assert hdr["GEOMCOR"] is True
    assert hdr["CRPIX3"] == 1.0
    assert hdr["CRVAL3"] == WAVE0
    assert hdr["CD3_3"] == DW
    assert hdr["WAVMID"] == 3550.0 and hdr["WAVGOOD0"] == 3450.0
    assert hdr["PXSCL"] == PXSCL and hdr["SLSCL"] == SLSCL
    assert hdr["ARCFL"] == ARC + ".fits"
    assert hdr["GEOMFL"] == ARC + "_geom.pkl"
    assert (hdr["CTYPE1"], hdr["CTYPE2"]) == ("RA---TAN", "DEC--TAN")


def test_pointing_sets_reference_value_and_unrotated_cd_matrix(
        run_primitive, synthetic_ccddata, tmp_redux):
    ccd = make_frame(synthetic_ccddata, tmp_redux, "BLUE")
    _, _, hdr = run_cube(run_primitive, ccd, tmp_redux)
    # source: analytic - RA 10:00:00 = 150 deg, Dec +02:00:00 = 2 deg
    assert hdr["CRVAL1"] == pytest.approx(150.0, abs=1e-9)
    assert hdr["CRVAL2"] == pytest.approx(2.0, abs=1e-9)
    # source: published (Calabretta & Greisen 2002, eq. 188) - rotation 0:
    # diagonal CD with cdelt1 = -SLSCL (RA increases to the left) and
    # cdelt2 = PXSCL
    assert hdr["CD1_1"] == pytest.approx(-SLSCL, rel=1e-12)
    assert hdr["CD2_2"] == pytest.approx(PXSCL, rel=1e-12)
    assert hdr["CD1_2"] == pytest.approx(0.0, abs=1e-15)
    assert hdr["CD2_1"] == pytest.approx(0.0, abs=1e-15)


def test_rotated_cd_matrix(run_primitive, synthetic_ccddata, tmp_redux):
    ccd = make_frame(synthetic_ccddata, tmp_redux, "BLUE")
    ccd.header["ROTPOSN"] = 25.0
    ccd.header["ROTREFAN"] = 5.0
    _, _, hdr = run_cube(run_primitive, ccd, tmp_redux, ROTOFF=10.0)
    # source: characterization @ 35e41e7 - the rotation angle is
    # crota = -(ROTPOSN + ROTREFAN + ROTOFF) = -40 deg (the sign convention
    # of the instrument's rotator is not in the docstring); the CD matrix
    # for that angle is published (Calabretta & Greisen 2002, eq. 188):
    # [cdelt1 cos, -cdelt2 sin; cdelt1 sin, cdelt2 cos]
    rho = math.radians(-40.0)
    c1, c2 = -SLSCL, PXSCL
    assert hdr["CD1_1"] == pytest.approx(c1 * math.cos(rho), rel=1e-12)
    assert hdr["CD1_2"] == pytest.approx(-c2 * math.sin(rho), rel=1e-12)
    assert hdr["CD2_1"] == pytest.approx(c1 * math.sin(rho), rel=1e-12)
    assert hdr["CD2_2"] == pytest.approx(c2 * math.cos(rho), rel=1e-12)
    assert hdr["IFUPA"] == 30.0                  # ROTPOSN + ROTREFAN


def test_missing_master_arc_makes_no_cube(run_primitive, synthetic_ccddata,
                                          tmp_redux, caplog):
    # source: characterization @ 35e41e7 - with no MARC row the primitive
    # logs an error, sets GEOMCOR False and writes no cube (see B34)
    ccd = synthetic_ccddata(IMTYPE="OBJECT", FRAMENO=1,
                            OFNAME=OBJ + ".fits")
    run_primitive(MakeCube, ccd)
    assert ccd.header["GEOMCOR"] is False
    assert not (tmp_redux.redux / (OBJ + "_icube.fits")).exists()
    assert "No reference geometry" in caplog.text
