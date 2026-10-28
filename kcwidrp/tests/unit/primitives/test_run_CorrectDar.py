"""CorrectDar run through ``run_primitive`` (task T.4, numeric case).

(The T.3 file ``test_CorrectDar.py`` covers the module function atm_disper
against the published Ciddor (1996) / Peck & Reeder (1972) formula.)

Specification (class docstring, plan T.4 and M.4): each wavelength plane of
the cube is moved by the atmospheric dispersion between the plane's
wavelength and the reference wavelength WAVMID, ``atm_disper(WAVMID, wl,
AIRMASS)`` arcsec (Mauna Kea defaults: 10 C, 61100 Pa, 50 % RH, 400 ppm), so
that a point source ends up at the same spatial position at every wavelength.
The dispersion is along the parallactic direction; on the IFU it is projected
with the angle ``IFUPA - PARANG``: ``sin`` of it along the slice-to-slice
axis x (cube axis 2, SLSCL deg/px) and ``cos`` along the slice axis y (cube
axis 1, PXSCL deg/px).  The cube is padded by DARPADX / DARPADY pixels on
each side first, so the source's input position (y0, x0) becomes (y0 +
DARPADY, x0 + DARPADX).  Plane j has the wavelength CRVAL3 + j * CD3_3
(FITS WCS paper I with CRPIX3 = 1, as MakeCube writes it).

What is measured: the centroid of a Gaussian point source in every output
plane.  Expected: input position + padding + (dispersion / pixel scale)
projected as above.  The dispersion in arcsec comes from ``atm_disper``
(pinned to the published formula by T.3); everything else -- wavelength of
each plane, projection, pixel scales, the shift of the light -- is computed
here from the header, independently of the primitive's arithmetic.

Reuse (M.4): ``dar_reference_case()`` returns the inputs and the expected
per-plane source positions; the warp path's DAR folded into ``W`` must
place the source at the same positions (within the stated tolerance).
"""
import math

import numpy as np
import pytest
from astropy.io import fits

from kcwidrp.primitives.CorrectDar import CorrectDar, atm_disper
from kcwidrp.primitives.kcwi_file_primitives import KCCDData

HEAD = "c0fed5c"

# Inputs of the reference case (all in the cube header)
NW, NY, NX = 21, 40, 12            # planes, pixels along slice, slices
CRVAL3, CD3_3 = 3500.0, 100.0      # Angstrom; planes 3500 ... 5500
WAVMID = 4500.0                    # DAR reference wavelength
AIRMASS = 1.5
Y0, X0 = 20.0, 6.0                 # source position in the input cube (px)
SIG_Y, SIG_X = 1.6, 1.0            # Gaussian sigma (px)
PXSCL = 4.048e-05 * 2              # deg / px along slice (2x2 binning)
SLSCL = 3.7718e-04                 # deg / px across slices (Small slicer)
WAVMID_INDEX = int(round((WAVMID - CRVAL3) / CD3_3))   # plane 10


def dar_reference_case(projection_deg=0.0, parang=30.0):
    """Inputs and expected output source positions for the DAR reference.

    Returns a dict: ``header`` (cards to put on the cube), ``waves`` (plane
    wavelengths), ``dy_px`` and ``dx_px`` (expected shift of the source in
    each plane, output pixels, before padding), ``y0``/``x0`` (input
    position).  ``projection_deg`` is IFUPA - PARANG.
    """
    # source: published @ c0fed5c - FITS WCS paper I (Greisen & Calabretta
    # 2002): 1-based pixel p has wavelength CRVAL3 + (p - CRPIX3) * CD3_3; with
    # CRPIX3 = 1 the 0-based plane j is at CRVAL3 + j * CD3_3.
    waves = CRVAL3 + CD3_3 * np.arange(NW)
    # source: published @ c0fed5c - atm_disper is the Ciddor (1996) dispersion
    # (pinned by test_CorrectDar.py); the projection onto the IFU axes is the
    # one the class docstring and plan describe: sin(IFUPA - PARANG) across
    # slices (x), cos along the slice (y).
    ang = math.radians(projection_deg)
    y_scale = PXSCL * 3600.0          # arcsec / px
    x_scale = SLSCL * 3600.0
    disp = np.array([atm_disper(WAVMID, w, AIRMASS) for w in waves])
    header = {"CRVAL3": CRVAL3, "CD3_3": CD3_3, "CRPIX3": 1.0,
              "WAVMID": WAVMID, "WAVGOOD0": waves[0], "WAVGOOD1": waves[-1],
              "AIRMASS": AIRMASS, "IFUPA": parang + projection_deg,
              "PARANG": parang, "PXSCL": PXSCL, "SLSCL": SLSCL,
              "GEOMCOR": True}
    return dict(header=header, waves=waves, y0=Y0, x0=X0,
                dy_px=disp * math.cos(ang) / y_scale,
                dx_px=disp * math.sin(ang) / x_scale)


def make_cube(synthetic_ccddata, case):
    """A KCCDData cube holding one Gaussian point source in every plane."""
    frame = synthetic_ccddata(**case["header"])
    yy, xx = np.mgrid[0:NY, 0:NX]
    img = np.exp(-0.5 * (((yy - case["y0"]) / SIG_Y) ** 2 +
                         ((xx - case["x0"]) / SIG_X) ** 2))
    cube = np.repeat(1000.0 * img[np.newaxis, :, :], NW, axis=0)
    ccd = KCCDData(cube, meta=frame.header, unit="electron",
                   uncertainty=np.ones(cube.shape))
    ccd.mask = np.zeros(cube.shape, dtype=np.uint8)
    ccd.flags = np.zeros(cube.shape, dtype=np.uint8)
    return ccd


def centroids(cube):
    """First-moment (y, x) centroid of each plane (0-based pixels)."""
    yy, xx = np.mgrid[0:cube.shape[1], 0:cube.shape[2]]
    tot = cube.sum(axis=(1, 2))
    return ((cube * yy).sum(axis=(1, 2)) / tot,
            (cube * xx).sum(axis=(1, 2)) / tot)


def run(run_primitive, synthetic_ccddata, case):
    cube = make_cube(synthetic_ccddata, case)
    out = run_primitive(CorrectDar, cube)
    return out.ccddata


@pytest.mark.parametrize("projection_deg", [0.0, 90.0, 30.0, -135.0])
def test_point_source_is_moved_by_the_dispersion(run_primitive,
                                                 synthetic_ccddata,
                                                 projection_deg):
    case = dar_reference_case(projection_deg)
    out = run(run_primitive, synthetic_ccddata, case)
    pady, padx = out.header["DARPADY"], out.header["DARPADX"]
    cy, cx = centroids(np.asarray(out.data, dtype=np.float64))
    # source: independent @ c0fed5c - see dar_reference_case
    exp_y = case["y0"] + pady + case["dy_px"]
    exp_x = case["x0"] + padx + case["dx_px"]
    # Tolerance: the first moment of a cubic-spline-shifted, well sampled
    # Gaussian (sigma >= 1 px, far from the edges) moves by the shift up to
    # the spline's interpolation error plus float32 storage; measured max
    # 1.8e-6 px at c0fed5c.  1e-4 px is ~50x that and >= 4000x below the
    # shifts tested (-2.7 ... +1.9 px along y, up to 0.58 px across slices).
    np.testing.assert_allclose(cy, exp_y, rtol=0, atol=1e-4)
    np.testing.assert_allclose(cx, exp_x, rtol=0, atol=1e-4)


def test_shift_is_large_enough_to_be_a_test():
    # source: analytic @ c0fed5c - guards the reference case: across 3500-5500
    # A at airmass 1.5 the dispersion is several pixels along the slice, so a
    # missing or sign-flipped correction cannot pass the 0.01 px tolerance.
    case = dar_reference_case(0.0)
    assert case["dy_px"][0] < -1.0 and case["dy_px"][-1] > 1.0
    assert case["dy_px"][WAVMID_INDEX] == 0.0


def test_reference_plane_is_not_moved(run_primitive, synthetic_ccddata):
    # source: analytic @ c0fed5c - atm_disper(w, w) = 0, so the plane at WAVMID
    # keeps the source at its input position plus the padding (spline shift by
    # 0 is the identity, up to float32 storage).
    case = dar_reference_case(30.0)
    out = run(run_primitive, synthetic_ccddata, case)
    cy, cx = centroids(np.asarray(out.data, dtype=np.float64))
    # Tolerance: float32 storage of the plane, ~1e-6 px.
    assert cy[WAVMID_INDEX] == pytest.approx(
        Y0 + out.header["DARPADY"], abs=1e-5)
    assert cx[WAVMID_INDEX] == pytest.approx(
        X0 + out.header["DARPADX"], abs=1e-5)


def test_header_and_padding(run_primitive, synthetic_ccddata, tmp_redux):
    case = dar_reference_case(30.0, parang=12.0)
    out = run(run_primitive, synthetic_ccddata, case)
    hdr = out.header
    # source: analytic @ c0fed5c - docstring
    assert hdr["DARCOR"] is True
    # source: analytic @ c0fed5c - DARANG records IFUPA - PARANG = 30 deg
    assert hdr["DARANG"] == pytest.approx(30.0, abs=1e-12)
    # source: analytic @ c0fed5c - docstring
    assert hdr["DAREFWL"] == WAVMID
    # source: characterization @ c0fed5c - padding is int(4"/scale) for a
    # BL frame (2" for H, 3" for M gratings): 4/0.29146 -> 13 px along the
    # slice, 4/1.3578 -> 2 px across slices.
    assert (hdr["DARPADY"], hdr["DARPADX"]) == (13, 2)
    assert out.data.shape == (NW, NY + 2 * 13, NX + 2 * 2)
    # source: analytic @ c0fed5c - the reference pixel moves with the padding
    assert hdr["CRPIX1"] == 1.0 + 2 and hdr["CRPIX2"] == 1.0 + 13
    # source: analytic @ c0fed5c - docstring - padded pixels carry flag 128.
    # In the reference plane (shift 0) the padding is exactly the border and
    # the data region keeps its input flags (0).
    ref = out.flags[WAVMID_INDEX]
    inner = np.zeros(ref.shape, dtype=bool)
    inner[13:13 + NY, 2:2 + NX] = True
    assert np.all(ref[inner] == 0) and np.all(ref[~inner] == 128)
    assert (tmp_redux.redux / "kb260101_00001_icubed.fits").exists()
    with fits.open(tmp_redux.redux / "kb260101_00001_icubed.fits") as hdul:
        assert hdul[0].header["DARCOR"] is True


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=(
    "B30: CorrectDar shifts the flags (and mask) cube with the same "
    "cubic-spline scipy.ndimage.shift as the data (CorrectDar.py:252-255): "
    "flag values become arbitrary bytes (e.g. 34, 94, 129, 141) at the "
    "padding edge, and padding pixels shifted in from outside the array get "
    "the fill value 0, i.e. are unflagged although they hold no data"))
def test_shifted_planes_keep_padding_flagged(run_primitive,
                                             synthetic_ccddata):
    # source: analytic @ c0fed5c - docstring - "Sets flags in padded region of
    # output cube to a value of 128."  With all input flags 0, every output
    # flag is 0 or 128, and the outer 5 rows at both ends of the slice axis
    # (padding is 13 px, the largest shift is < 4 px) hold no data, so they are
    # 128 in every plane.
    case = dar_reference_case(30.0)
    out = run(run_primitive, synthetic_ccddata, case)
    assert set(np.unique(out.flags)) <= {0, 128}
    assert np.all(out.flags[:, :5, :] == 128)
    assert np.all(out.flags[:, -5:, :] == 128)


def test_not_geometry_corrected_is_skipped(run_primitive, synthetic_ccddata):
    # source: analytic @ c0fed5c - docstring of _pre_condition - only
    # geometrically corrected cubes are DAR corrected; otherwise DARCOR is
    # False and data unchanged.
    case = dar_reference_case(0.0)
    case["header"]["GEOMCOR"] = False
    cube = make_cube(synthetic_ccddata, case)
    before = cube.data.copy()
    run_primitive(CorrectDar, cube)
    assert np.array_equal(cube.data, before)
    assert cube.header["DARCOR"] is False
