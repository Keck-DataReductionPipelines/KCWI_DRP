"""RemoveCosmicRays run through ``run_primitive`` (task T.4, detector level).

Specification (class docstring): cosmic rays are detected with astroscrappy,
the cleaned image replaces the data, the flags array gets the value 4 at the
cosmic-ray pixels, NCRCLEAN is the number of such pixels and CRCLEAN True.
Frames shorter than CRR_MINEXPTIME (60 s in kcwi.cfg) are not cleaned.

The input is a smooth, noise-free ramp (1000 + 0.5 col + 0.25 row) with
cosmic-ray hits added at known pixels: single-pixel spikes are the textbook
case a Laplacian detector must find, and a smooth ramp has no other
candidate.  The exact flagged set is asserted, and that no other pixel
changes.
"""
import numpy as np
import pytest

from kcwidrp.primitives.RemoveCosmicRays import RemoveCosmicRays


def ramp_frame(synthetic_ccddata, ttime=300.0, hits=(), **header):
    ccd = synthetic_ccddata(TTIME=ttime, TELAPSE=ttime, **header)
    ny, nx = 64, 64
    rr, cc = np.mgrid[0:ny, 0:nx]
    ccd.data = 1000.0 + 0.5 * cc + 0.25 * rr
    ccd.flags = np.zeros((ny, nx), dtype=np.uint8)
    for r, c, amp in hits:
        ccd.data[r, c] += amp
    return ccd


def test_injected_cosmic_rays_are_flagged_and_cleaned(run_primitive,
                                                      synthetic_ccddata):
    hits = [(20, 30, 8000.0), (45, 12, 5000.0), (10, 50, 20000.0)]
    ccd = ramp_frame(synthetic_ccddata, hits=hits)
    clean_truth = ramp_frame(synthetic_ccddata).data     # the ramp, no hits
    out = run_primitive(RemoveCosmicRays, ccd)
    expected_flags = np.zeros((64, 64), dtype=np.uint8)
    for r, c, _ in hits:
        expected_flags[r, c] = 4
    # source: independent @ cf88fb8 - exactly the three injected pixels are
    # flagged with value 4 (docstring)
    assert np.array_equal(out.ccddata.flags, expected_flags)
    hdr = out.ccddata.header
    assert hdr["NCRCLEAN"] == 3                      # source: independent
    assert hdr["CRCLEAN"] is True                    # source: docstring
    # source: independent @ cf88fb8 - every pixel that is not a hit is
    # bit-for-bit what it was
    mask = expected_flags == 0
    assert np.array_equal(out.ccddata.data[mask], clean_truth[mask])
    # source: analytic @ cf88fb8 - a hit is replaced by the mean of its
    # unmasked neighbours ('meanmask'); on a plane the mean of a set of
    # neighbours placed symmetrically about the pixel equals the pixel's
    # value, so the ramp is recovered.
    # Tolerance: float64 sums of ~25 values near 1000: ~1e-12 absolute.
    for r, c, _ in hits:
        np.testing.assert_allclose(out.ccddata.data[r, c], clean_truth[r, c],
                                   rtol=0, atol=1e-9)


def test_clean_smooth_image_is_untouched(run_primitive, synthetic_ccddata):
    # source: independent @ cf88fb8 - a ramp without cosmic rays has nothing
    # to flag or change
    ccd = ramp_frame(synthetic_ccddata)
    before = ccd.data.copy()
    out = run_primitive(RemoveCosmicRays, ccd)
    assert not out.ccddata.flags.any()
    assert np.array_equal(out.ccddata.data, before)
    assert out.ccddata.header["NCRCLEAN"] == 0


def test_existing_flags_kept_and_added_to(run_primitive, synthetic_ccddata):
    # source: characterization @ cf88fb8 - value 4 is added to the pixel's
    # existing flags (saturation 8 + cosmic ray 4 = 12); others are kept.
    ccd = ramp_frame(synthetic_ccddata, hits=[(20, 30, 8000.0)])
    ccd.flags[20, 30] = 8
    ccd.flags[5, 5] = 2
    out = run_primitive(RemoveCosmicRays, ccd)
    assert out.ccddata.flags[20, 30] == 12
    assert out.ccddata.flags[5, 5] == 2


def test_short_exposure_is_not_cleaned(run_primitive, synthetic_ccddata):
    # source: analytic @ cf88fb8 - docstring: below CRR_MINEXPTIME (60 s) no
    # cleaning: a 30 s frame keeps its hit, flags stay clear, CRCLEAN False.
    ccd = ramp_frame(synthetic_ccddata, ttime=30.0,
                     hits=[(20, 30, 8000.0)])
    before = ccd.data.copy()
    out = run_primitive(RemoveCosmicRays, ccd)
    assert np.array_equal(out.ccddata.data, before)
    assert not out.ccddata.flags.any()
    assert out.ccddata.header["CRCLEAN"] is False    # source: docstring
    assert out.ccddata.header["NCRCLEAN"] == 0       # source: docstring


def test_exposure_time_threshold_edge(run_primitive, synthetic_ccddata):
    # source: analytic @ cf88fb8 - "exposure time below which no CR cleaning
    # is done": exactly 60 s is not below, so it is cleaned.
    ccd = ramp_frame(synthetic_ccddata, ttime=60.0,
                     hits=[(20, 30, 8000.0)])
    out = run_primitive(RemoveCosmicRays, ccd)
    assert out.ccddata.flags[20, 30] == 4
    assert out.ccddata.header["CRCLEAN"] is True


def test_shuffled_exposure_uses_shuffle_time(run_primitive,
                                             synthetic_ccddata):
    # source: characterization @ cf88fb8 - the exposure time compared with
    # the minimum is max(TELAPSE, NSHFUP * TTIME): TELAPSE 30 s but 3
    # shuffles of 30 s make 90 s >= 60 s, so the frame is cleaned.
    ccd = ramp_frame(synthetic_ccddata, ttime=30.0, NSHFUP=3,
                     hits=[(20, 30, 8000.0)])
    out = run_primitive(RemoveCosmicRays, ccd)
    assert out.ccddata.header["CRCLEAN"] is True
    assert out.ccddata.flags[20, 30] == 4
