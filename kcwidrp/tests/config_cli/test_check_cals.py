"""Characterization tests for kcwidrp.scripts.check_cals.check_cal_type
(task T.3), driven by a hand-built Proctab.

check_cal_type(proctab, ccd_frame, setup_frame, targ_type, minimum, logger)
returns "PASSED" when the proc table holds at least ``minimum`` rows of type
``targ_type`` matching ``ccd_frame``'s configuration, otherwise
"FAILED: found <n>, needed <m>" (both numbers left-justified in width 3).
"""
import logging

import pytest

from kcwidrp.core.kcwi_proctab import Proctab
from kcwidrp.scripts.check_cals import check_cal_type

LOG = logging.getLogger("test_check_cals")


@pytest.fixture
def table(synthetic_ccddata):
    """Proc table with 3 BIAS (2x2), 1 BIAS (1x1), 2 CONTBARS of the
    default state and 1 CONTBARS of another state."""
    pt = Proctab()
    pt.new_proctab()

    def add(fno, typ, **kw):
        f = synthetic_ccddata(FRAMENO=fno, **kw)
        pt.update_proctab(f, suffix="raw", newtype=typ,
                          filename="f%d.fits" % fno)

    for fno in (1, 2, 3):
        add(fno, "BIAS")
    add(4, "BIAS", binning=(1, 1))
    add(5, "CONTBARS")
    add(6, "CONTBARS")
    add(7, "CONTBARS", STATEID="otherstate")
    return pt


def _setup(frame):
    # a proc-table row stand-in; check_cal_type only reads setup_frame['CID']
    return {"CID": frame.header["STATEID"]}


def test_enough_biases_passes(table, synthetic_ccddata):
    # source: analytic - 3 BIAS rows have this frame's CCDCFG (binning 2x2)
    f = synthetic_ccddata()
    assert check_cal_type(table, f, _setup(f), "BIAS", 3, LOG) == "PASSED"
    assert check_cal_type(table, f, _setup(f), "BIAS", 1, LOG) == "PASSED"


def test_too_few_biases_fails_with_counts(table, synthetic_ccddata):
    # source: analytic - 3 found, 5 needed; the numbers are formatted '<3'
    f = synthetic_ccddata()
    assert check_cal_type(table, f, _setup(f), "BIAS", 5, LOG) == \
        "FAILED: found 3  , needed 5  "


def test_other_binning_counted_separately(table, synthetic_ccddata):
    # source: analytic - a 1x1 frame matches only the single 1x1 bias
    f = synthetic_ccddata(binning=(1, 1))
    assert check_cal_type(table, f, _setup(f), "BIAS", 1, LOG) == "PASSED"
    assert check_cal_type(table, f, _setup(f), "BIAS", 2, LOG) == \
        "FAILED: found 1  , needed 2  "


def test_contbars_matched_on_state_id(table, synthetic_ccddata):
    # source: analytic - CONTBARS rows are matched on STATEID: 2 for the
    # default state, 1 for 'otherstate'
    f = synthetic_ccddata()
    g = synthetic_ccddata(STATEID="otherstate")
    assert check_cal_type(table, f, _setup(f), "CONTBARS", 2, LOG) == "PASSED"
    assert check_cal_type(table, f, _setup(f), "CONTBARS", 3, LOG) == \
        "FAILED: found 2  , needed 3  "
    assert check_cal_type(table, g, _setup(g), "CONTBARS", 1, LOG) == "PASSED"
    assert check_cal_type(table, g, _setup(g), "CONTBARS", 2, LOG) == \
        "FAILED: found 1  , needed 2  "


def test_missing_type_fails_with_zero_found(table, synthetic_ccddata):
    # source: analytic - no ARCLAMP rows in the table
    f = synthetic_ccddata()
    assert check_cal_type(table, f, _setup(f), "ARCLAMP", 1, LOG) == \
        "FAILED: found 0  , needed 1  "


def test_zero_minimum_always_passes(table, synthetic_ccddata):
    # source: analytic - 0 found >= 0 needed
    f = synthetic_ccddata()
    assert check_cal_type(table, f, _setup(f), "ARCLAMP", 0, LOG) == "PASSED"
