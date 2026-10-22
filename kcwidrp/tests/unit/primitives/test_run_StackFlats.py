"""StackFlats run through ``run_primitive`` (task T.4, detector level).

(The T.3 file ``test_StackFlats.py`` covers the module function inv_median.)

Specification (class docstring and ``ccdproc`` semantics): the flats of one
state (proc-table rows of TYPE ``want_type`` sharing the frame's STATEID),
once ``min_files`` are present and no stack of ``stack_type`` exists yet, are
each scaled to unit median, combined with 'average' and an upper 2-sigma
clip, and the result is scaled back up by the median of the flats' medians.
The stack is written as ``redux/<first flat stem>_<stack_type lower>.fits``,
entered in the proc table as ``stack_type``, with the mask of a flat; the
header records NSTACK, STCKMETH, STCKSIGU and STACKFn.

Inputs are the ``*_intd.fits`` files SubtractScatteredLight leaves in
``redux/`` (written here with ``kcwi_fits_writer`` and a proc-table row).
"""
import numpy as np
import pytest
from astropy.io import fits

from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer
from kcwidrp.primitives.StackFlats import StackFlats

SHAPE = (16, 16)
EXTRA = {"want_type": "FLATLAMP", "stack_type": "SFLAT", "new_type": "MFLAT"}


def write_flats(synthetic_ccddata, tmp_redux, datas, masks=None,
                imtype="FLATLAMP"):
    """Write redux/<stem>_intd.fits for each array and register it."""
    for i, data in enumerate(datas):
        stem = "kb260101_%05d" % (i + 1)
        f = synthetic_ccddata(IMTYPE=imtype, FRAMENO=i + 1,
                              OFNAME=stem + ".fits", TTIME=60.0)
        f.data = np.asarray(data, dtype=np.float64)
        f.uncertainty = None
        f.mask = None if masks is None else masks[i]
        f.flags = None
        f.header["HISTORY"] = "synthetic intd"
        kcwi_fits_writer(f, output_file=stem + ".fits",
                         output_dir=str(tmp_redux.redux), suffix="intd")
        tmp_redux.proctab.update_proctab(frame=f, suffix="intd",
                                         filename=stem + ".fits")


def stack(run_primitive, synthetic_ccddata, min_files, **extra):
    ccd = synthetic_ccddata(IMTYPE="FLATLAMP", FRAMENO=99,
                            OFNAME="kb260101_00099.fits")
    ex = dict(EXTRA, min_files=min_files)
    ex.update(extra)
    return run_primitive(StackFlats, ccd, action_extra=ex)


def read(tmp_redux, name):
    with fits.open(tmp_redux.redux / name) as hdul:
        unc = hdul["UNCERT"].data.copy() if "UNCERT" in hdul else None
        mask = hdul["MASK"].data.copy() if "MASK" in hdul else None
        return hdul[0].data.copy(), unc, mask, hdul[0].header.copy()


def test_identical_flats_stack_to_the_flat(run_primitive, synthetic_ccddata,
                                           tmp_redux):
    rr, cc = np.mgrid[0:16, 0:16]
    flat = 1024.0 * (1.0 + cc / 16.0)           # exact in float32
    write_flats(synthetic_ccddata, tmp_redux, [flat] * 5)
    out = stack(run_primitive, synthetic_ccddata, 5)
    assert out is not None
    data, unc, _, hdr = read(tmp_redux, "kb260101_00001_sflat.fits")
    # source: analytic @ cf88fb8 - five identical flats: each is scaled to
    # unit median, the average of identical frames is that frame, and the
    # scale-back by the median of medians restores the flat.
    # Tolerance: the scale by 1/median and back is two float roundings in
    # float64 then float32 storage: relative 1e-6 is >> 6e-8.
    np.testing.assert_allclose(data, flat, rtol=1e-6)
    # source: analytic @ cf88fb8 - identical inputs have zero spread.
    # Tolerance: 1/median scaling leaves rounding differences ~1e-16 between
    # nominally identical frames (seen 5e-17); 1e-9 is far below any real
    # spread (flats scatter at ~1e-2 of unit median).
    assert np.abs(unc).max() < 1e-9
    assert hdr["NSTACK"] == 5                  # source: docstring
    assert hdr["STCKMETH"] == "average"        # source: docstring
    assert hdr["STCKSIGU"] == 2.0              # source: docstring
    assert hdr["IMTYPE"] == "SFLAT"            # source: docstring
    for i in range(5):
        assert hdr["STACKF%d" % (i + 1)] == "kb260101_%05d_intd.fits" % (
            i + 1)                             # source: docstring


def test_output_name_and_proctab_entry(run_primitive, synthetic_ccddata,
                                       tmp_redux):
    write_flats(synthetic_ccddata, tmp_redux, [np.full(SHAPE, 1000.0)] * 3)
    stack(run_primitive, synthetic_ccddata, 3)
    # source: analytic @ cf88fb8 - docstring: FLATLAMP -> *_sflat.fits, named
    # from the first flat, with an SFLAT row in the proc table
    assert (tmp_redux.redux / "kb260101_00001_sflat.fits").is_file()
    tab = tmp_redux.proctab.proctab
    rows = tab[tab["TYPE"] == "SFLAT"]
    assert len(rows) == 1
    assert rows["SUFF"][0] == "sflat"
    assert rows["filename"][0] == "kb260101_00001.fits"
    assert tmp_redux.proctab_file.is_file()


def test_too_few_flats_does_nothing(run_primitive, synthetic_ccddata,
                                    tmp_redux):
    # source: analytic @ cf88fb8 - docstring: needs min_files flats; with 2
    # of 3 the precondition fails and apply() returns None, writing nothing
    write_flats(synthetic_ccddata, tmp_redux, [np.full(SHAPE, 1000.0)] * 2)
    assert stack(run_primitive, synthetic_ccddata, 3) is None
    assert not list(tmp_redux.redux.glob("*_sflat.fits"))
    assert not (tmp_redux.proctab.proctab["TYPE"] == "SFLAT").any()


def test_existing_stack_blocks_a_second_one(run_primitive, synthetic_ccddata,
                                            tmp_redux):
    # source: analytic @ cf88fb8 - precondition: no stack of this type may
    # exist yet ("already have ..., expecting 0"), so a second run is a no-op
    write_flats(synthetic_ccddata, tmp_redux, [np.full(SHAPE, 1000.0)] * 3)
    assert stack(run_primitive, synthetic_ccddata, 3) is not None
    (tmp_redux.redux / "kb260101_00001_sflat.fits").unlink()
    assert stack(run_primitive, synthetic_ccddata, 3) is None
    assert not list(tmp_redux.redux.glob("*_sflat.fits"))


def test_other_state_flats_are_not_counted(run_primitive, synthetic_ccddata,
                                           tmp_redux):
    # source: analytic @ cf88fb8 - flats are selected by STATEID: 3 flats of
    # the frame's state plus 2 of another state, min_files 5 -> not enough.
    write_flats(synthetic_ccddata, tmp_redux, [np.full(SHAPE, 1000.0)] * 3)
    for j in (4, 5):
        f = synthetic_ccddata(IMTYPE="FLATLAMP", FRAMENO=j,
                              OFNAME="kb260101_%05d.fits" % j,
                              STATEID="otherstate0000000000000")
        f.data = np.full(SHAPE, 1000.0)
        tmp_redux.proctab.update_proctab(frame=f, suffix="intd",
                                         filename="kb260101_%05d.fits" % j)
    assert stack(run_primitive, synthetic_ccddata, 5) is None


def test_mask_of_a_flat_is_carried_to_the_stack(run_primitive,
                                                synthetic_ccddata,
                                                tmp_redux):
    # source: analytic @ cf88fb8 - docstring: the bad-pixel mask (the same in
    # every flat) becomes the stack's mask.
    mask = np.zeros(SHAPE, dtype=bool)
    mask[3, 4] = True
    mask[9, 9] = True
    write_flats(synthetic_ccddata, tmp_redux, [np.full(SHAPE, 1000.0)] * 3,
                masks=[mask] * 3)
    stack(run_primitive, synthetic_ccddata, 3)
    _, _, m, _ = read(tmp_redux, "kb260101_00001_sflat.fits")
    assert np.array_equal(m.astype(bool), mask)


# -- flats of different brightness ---------------------------------------------

LEVELS = [1.0, 2.0, 3.0, 4.0, 5.0]           # relative brightness a_i
DELTAS = [0.0, 0.01, 0.02, 0.03, 0.04]       # per-flat bump at three pixels
BUMP = [(2, 2), (7, 11), (12, 5)]            # < half the pixels: median 1000


def bright_flats():
    flats = []
    for a, d in zip(LEVELS, DELTAS):
        f = np.full(SHAPE, 1000.0 * a)
        for r, c in BUMP:
            f[r, c] = 1000.0 * a * (1.0 + d)
        flats.append(f)
    return flats


@pytest.mark.xfail(strict=True, reason="B26: StackFlats.py:92-95 builds "
                   "every flat's path from stackf[0], so all N inputs are the "
                   "first flat and the other flats are never read")
def test_flats_of_different_brightness_are_scaled_and_averaged(
        run_primitive, synthetic_ccddata, tmp_redux):
    write_flats(synthetic_ccddata, tmp_redux, bright_flats())
    stack(run_primitive, synthetic_ccddata, 5)
    data, _, _, _ = read(tmp_redux, "kb260101_00001_sflat.fits")
    # source: analytic @ cf88fb8 - flat i is a_i * 1000 * (1 + d_i at the
    # three bump pixels).  Each has median a_i * 1000 (the bumps are fewer
    # than half the pixels), so after unit-median scaling away from the
    # bumps every flat is 1; at a bump the five values are 1 + d_i, whose
    # average is 1.02 (none is clipped: median 1.02, MAD 0.01 -> sigma
    # 0.0148, limit 1.0497 > 1.04).  Scale-back by the median of the medians
    # 3000 gives 3000 away from the bumps and 3060 at them.
    # Tolerance: float32 storage and the scale round trip, relative 1e-6.
    expected = np.full(SHAPE, 3000.0)
    for r, c in BUMP:
        expected[r, c] = 3060.0
    np.testing.assert_allclose(data, expected, rtol=1e-6)


@pytest.mark.xfail(strict=True, reason="B26: all N inputs are the first "
                   "flat (StackFlats.py:92-95), so the uncertainty is 0; and "
                   "B27: StackFlats.py:115 rescales the stacked data "
                   "by median(fmeds) but not its uncertainty, which stays in "
                   "unit-median units")
def test_stack_uncertainty_is_in_the_units_of_the_stack(
        run_primitive, synthetic_ccddata, tmp_redux):
    write_flats(synthetic_ccddata, tmp_redux, bright_flats())
    stack(run_primitive, synthetic_ccddata, 5)
    _, unc, _, _ = read(tmp_redux, "kb260101_00001_sflat.fits")
    # source: analytic @ cf88fb8 - at a bump the five unit-median values are
    # 1 + d_i with population variance 2e-4 (d = 0, .01, ..., .04 about their
    # mean .02: (4+1+0+1+4) 1e-4 / 5), so the standard error of their
    # average is sqrt(2e-4 / 5) = 6.3246e-3 in unit-median units, and
    # 3000 * 6.3246e-3 = 18.974 in the units of the stack (scaled by the
    # same 3000 as the data).  Elsewhere the five values are equal: 0.
    # Tolerance: float32 storage and the scale round trip, relative 1e-5.
    expected = np.zeros(SHAPE)
    for r, c in BUMP:
        expected[r, c] = 3000.0 * np.sqrt(2e-4 / 5.0)
    np.testing.assert_allclose(unc, expected, rtol=1e-5, atol=1e-9)
