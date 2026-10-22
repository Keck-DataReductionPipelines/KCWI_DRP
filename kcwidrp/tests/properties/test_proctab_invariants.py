"""Property: a frame added to the proc table is found by search_proctab (T.5).

The proc table records every processed frame; later primitives look up
their calibrations with ``search_proctab(frame, target_type, target_group,
nearest)``.  Whatever else the table holds, a frame that was added with
``update_proctab(frame, newtype=T)`` must be among the rows returned when
the same frame searches for type T: it has, by construction, the same
camera, CCD configuration (CCDCFG), exposure time, group and instrument
state as itself, which are the keys the documented matching rules use
(biases and master darks: CCDCFG; raw darks: CCDCFG and TTIME, and the group
when given; objects: the group; other types: STATEID).  With ``nearest`` it
is the row at MJD offset 0, so it is the one returned (MJDs are distinct).

Tables are random sets of 1-8 frames of random type, group, state, CCD
configuration (binning), exposure time and MJD, all of one camera.  Bug B17
(search_proctab raises IndexError when the table holds both cameras) is
covered separately by a strict xfail with mixed cameras.  T.3 checks fixed
hand-built tables.
"""
import pytest
from hypothesis import Phase, given, settings
from hypothesis import strategies as st

from kcwidrp.core.kcwi_proctab import Proctab
from kcwidrp.tests.properties._profiles import scaled

TYPES = ["MBIAS", "BIAS", "DARK", "MDARK", "OBJECT", "MARC", "MFLAT",
         "CONTBARS"]
GROUPS = ["kb260101_00001", "kb260101_00010", "kb260101_00020"]
STATES = ["6968443881c6cc2a85772b8f", "otherstate1", "otherstate2"]


@st.composite
def frame_specs(draw, cameras):
    """1-8 frame descriptions with distinct FRAMENO and MJD."""
    n = draw(st.integers(1, 8))
    framenos = draw(st.lists(st.integers(1, 9999), min_size=n, max_size=n,
                             unique=True))
    mjds = draw(st.lists(st.integers(0, 10**6), min_size=n, max_size=n,
                         unique=True))
    specs = []
    for fno, mjd in zip(framenos, mjds):
        specs.append(dict(
            camera=draw(st.sampled_from(cameras)),
            binning=draw(st.sampled_from([(1, 1), (2, 2)])),
            type=draw(st.sampled_from(TYPES)),
            FRAMENO=fno,
            MJD=61000.0 + mjd * 1e-4,
            GROUPID=draw(st.sampled_from(GROUPS)),
            STATEID=draw(st.sampled_from(STATES)),
            TTIME=float(draw(st.sampled_from([0, 1, 30, 300]))),
            OFNAME="f%04d.fits" % fno,
        ))
    return specs


def _build(synthetic_ccddata, specs):
    pt = Proctab()
    pt.new_proctab()
    frames = []
    for s in specs:
        hdr = {k: v for k, v in s.items()
               if k not in ("camera", "binning", "type")}
        f = synthetic_ccddata(camera=s["camera"], binning=s["binning"], **hdr)
        pt.update_proctab(f, suffix="raw", newtype=s["type"],
                          filename=s["OFNAME"])
        frames.append(f)
    return pt, frames


def _check_found(pt, frame, spec, nearest):
    tab = pt.search_proctab(frame, target_type=spec["type"],
                            target_group=frame.header["GROUPID"],
                            nearest=nearest)
    names = list(tab["filename"])
    # source: analytic - the frame matches itself on every documented key
    assert spec["OFNAME"] in names
    if nearest:
        # source: analytic - MJD offset 0 is the unique minimum
        assert names == [spec["OFNAME"]]


@given(specs=frame_specs(cameras=["BLUE"]) | frame_specs(cameras=["RED"]),
       pick=st.integers(0, 7), nearest=st.booleans())
@scaled(0.6)
def test_added_frame_is_found_single_camera(synthetic_ccddata, specs, pick,
                                            nearest):
    pt, frames = _build(synthetic_ccddata, specs)
    i = pick % len(specs)
    _check_found(pt, frames[i], specs[i], nearest)


@st.composite
def mixed_camera_specs(draw):
    specs = draw(frame_specs(cameras=["BLUE", "RED"]).filter(
        lambda s: len({x["camera"] for x in s}) == 2))
    return specs


@pytest.mark.xfail(strict=True, reason="B17: search_proctab raises "
                   "IndexError when the table holds frames of both cameras "
                   "(kcwi_proctab.py:154 filters the camera-filtered table "
                   "with a mask built from the full table)")
@given(specs=mixed_camera_specs(), pick=st.integers(0, 7),
       nearest=st.booleans())
# Every mixed table triggers B17, so shrinking would only spend time (about
# 100 s) minimising a known failure: generate phase only.
@settings(phases=[Phase.generate], max_examples=5)
def test_added_frame_is_found_mixed_cameras(synthetic_ccddata, specs, pick,
                                            nearest):
    pt, frames = _build(synthetic_ccddata, specs)
    i = pick % len(specs)
    _check_found(pt, frames[i], specs[i], nearest)
