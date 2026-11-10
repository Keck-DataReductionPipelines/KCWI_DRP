"""Characterization tests for kcwidrp.core.kcwi_proctab.Proctab (task T.3).

Source of expected values: hand-built tables.  Every frame is built from the
``synthetic_ccddata`` fixture with explicit header overrides; the expected
table contents are the override values written in the test (source:
analytic - the row is the header values that were put in).
"""
import numpy as np
import pytest

from kcwidrp.core.kcwi_proctab import Proctab


@pytest.fixture
def make_frame(synthetic_ccddata):
    def make(camera="BLUE", **kw):
        return synthetic_ccddata(camera=camera, **kw)
    return make


def _new():
    pt = Proctab()
    pt.new_proctab()
    return pt


def _add(pt, frame, newtype, fname, suffix="mbias"):
    pt.update_proctab(frame, suffix=suffix, newtype=newtype, filename=fname)


def test_new_proctab_columns_and_empty():
    # source: analytic - column names/order are those new_proctab declares in
    # the KCWI proc-table definition (19 columns).
    pt = _new()
    assert len(pt.proctab) == 0
    assert pt.proctab.colnames == [
        'FRAMENO', 'CID', 'DID', 'TYPE', 'GRPID', 'TTIME', 'CAM', 'IFU',
        'GRAT', 'GANG', 'CWAVE', 'BIN', 'FILT', 'MJD', 'STAGE', 'SUFF',
        'OFNAME', 'TARGNAME', 'filename']
    # string columns must be object so that long names are not truncated
    assert pt.proctab['filename'].dtype == object
    assert pt.proctab['CID'].dtype == object


def test_update_adds_row_with_header_values(make_frame):
    # source: analytic - each cell is the header value given in the overrides
    # (CCDCFG of the fixture is "2211000": CCDSUM 22, CCDMODE 1, GAINMUL 10,
    # AMPMNUM 00).
    pt = _new()
    f = make_frame(MJD=100.5, FRAMENO=7, TTIME=12.5, GROUPID="grpA",
                   OFNAME="kb260101_00007.fits", TARGNAME="Hz 44")
    _add(pt, f, "MBIAS", "kb_7.fits")
    assert len(pt.proctab) == 1
    row = pt.proctab[0]
    assert row['FRAMENO'] == 7
    assert row['CID'] == "6968443881c6cc2a85772b8f"   # STATEID default
    assert row['DID'] == 2211000
    assert row['TYPE'] == "MBIAS"                      # newtype
    assert row['GRPID'] == "grpA"
    assert row['TTIME'] == 12.5
    assert row['CAM'] == "BLUE"
    assert row['IFU'] == "Small"
    assert row['GRAT'] == "BL"
    assert row['BIN'] == "2,2"
    assert row['MJD'] == 100.5
    assert row['STAGE'] == 1                           # suffix 'mbias'
    assert row['SUFF'] == "mbias"
    assert row['OFNAME'] == "kb260101_00007.fits"
    assert row['TARGNAME'] == "Hz44"                   # spaces removed
    assert row['filename'] == "kb_7.fits"


@pytest.mark.parametrize("suffix, stage", [
    ("raw", 9),       # not in the stage table
    ("mbias", 1), ("intd", 2), ("mcbars", 3), ("mflat", 4), ("intk", 5),
    ("icube", 6), ("icubed", 7), ("icubes", 8), ("somethingelse", 9),
])
def test_stage_from_suffix(make_frame, suffix, stage):
    # source: characterization @ b0e2550 - the suffix->stage numbers are what
    # update_proctab did at this commit (no external specification exists).
    # 'raw' is lower-case and so is not the table key 'RAW'.
    pt = _new()
    _add(pt, make_frame(), "OBJECT", "x.fits", suffix=suffix)
    assert pt.proctab['STAGE'][0] == stage


def test_stateid_zero_and_empty_group_become_none(make_frame):
    # source: characterization @ b0e2550 - STATEID '0' and a blank GROUPID
    # were stored as 'NONE'.
    pt = _new()
    f = make_frame(STATEID="0", GROUPID="   ")
    _add(pt, f, "BIAS", "a.fits")
    assert pt.proctab['CID'][0] == "NONE"
    assert pt.proctab['GRPID'][0] == "NONE"


def test_target_name_falls_back_to_object(make_frame):
    # source: analytic - empty TARGNAME -> OBJECT keyword, spaces removed
    pt = _new()
    _add(pt, make_frame(TARGNAME="", OBJECT="My Object"), "OBJECT", "a.fits")
    assert pt.proctab['TARGNAME'][0] == "MyObject"


def test_red_camera_row_uses_red_keywords(make_frame):
    # source: analytic - RED rows take RGRATNAM/RGRANGLE/RCWAVE and no filter
    # (the fixture's RED header carries RGRATNAM 'RM1', RGRANGLE 123.0774384,
    # RCWAVE 7199.9848633).
    pt = _new()
    _add(pt, make_frame("RED"), "OBJECT", "r.fits")
    row = pt.proctab[0]
    assert row['CAM'] == "RED"
    assert row['GRAT'] == "RM1"
    assert row['GANG'] == 123.0774384
    assert row['CWAVE'] == 7199.9848633


def test_same_key_replaces_earlier_row_and_rows_sorted_by_frameno(make_frame):
    # source: analytic - rows are unique on (CID, FRAMENO, TYPE), keeping the
    # last added, and the table is sorted on FRAMENO.
    pt = _new()
    _add(pt, make_frame(FRAMENO=5), "MBIAS", "five_a.fits")
    _add(pt, make_frame(FRAMENO=2), "MBIAS", "two.fits")
    _add(pt, make_frame(FRAMENO=5), "MBIAS", "five_b.fits")
    assert list(pt.proctab['FRAMENO']) == [2, 5]
    assert list(pt.proctab['filename']) == ["two.fits", "five_b.fits"]


def _bias_table(make_frame):
    pt = _new()
    for fno, mjd, grp, fname in [(1, 100.1, "g1", "b1.fits"),
                                 (2, 100.5, "g1", "b2.fits"),
                                 (3, 101.0, "g2", "b3.fits")]:
        _add(pt, make_frame(FRAMENO=fno, MJD=mjd, GROUPID=grp), "MBIAS", fname)
    _add(pt, make_frame(FRAMENO=4, MJD=100.2, TTIME=30.0, GROUPID="g1"),
         "DARK", "d4.fits")
    _add(pt, make_frame(FRAMENO=5, MJD=100.3, TTIME=60.0, GROUPID="g1"),
         "DARK", "d5.fits")
    _add(pt, make_frame(FRAMENO=6, MJD=100.4, GROUPID="g1"), "MDARK",
         "md6.fits")
    _add(pt, make_frame(FRAMENO=7, MJD=100.6, GROUPID="g1"), "OBJECT",
         "o7.fits")
    _add(pt, make_frame(FRAMENO=8, MJD=100.7, GROUPID="g2"), "OBJECT",
         "o8.fits")
    _add(pt, make_frame(FRAMENO=9, MJD=100.8, STATEID="otherstate"), "MARC",
         "arc9.fits")
    return pt


def test_search_bias_returns_all_matching_ccdcfg(make_frame):
    # source: analytic - three MBIAS rows were added with the frame's CCDCFG
    pt = _bias_table(make_frame)
    tab = pt.search_proctab(make_frame(MJD=100.45), target_type="MBIAS")
    assert list(tab['filename']) == ["b1.fits", "b2.fits", "b3.fits"]


def test_search_bias_other_ccdcfg_not_matched(make_frame):
    # source: analytic - a frame with different binning has CCDCFG '1111000'
    pt = _bias_table(make_frame)
    other = make_frame(binning=(1, 1))
    tab = pt.search_proctab(other, target_type="MBIAS")
    assert len(tab) == 0


def test_search_nearest_picks_smallest_mjd_offset(make_frame):
    # source: analytic - |MJD - 100.45| is 0.35, 0.05, 0.55 for b1, b2, b3,
    # so b2 is nearest.
    pt = _bias_table(make_frame)
    tab = pt.search_proctab(make_frame(MJD=100.45), target_type="MBIAS",
                            nearest=True)
    assert list(tab['filename']) == ["b2.fits"]


def test_search_nearest_below_and_above_range(make_frame):
    # source: analytic - 99.0 is nearest to 100.1 (b1); 105.0 to 101.0 (b3)
    pt = _bias_table(make_frame)
    lo = pt.search_proctab(make_frame(MJD=99.0), target_type="MBIAS",
                           nearest=True)
    hi = pt.search_proctab(make_frame(MJD=105.0), target_type="MBIAS",
                           nearest=True)
    assert list(lo['filename']) == ["b1.fits"]
    assert list(hi['filename']) == ["b3.fits"]


def test_search_nearest_false_does_not_reduce(make_frame):
    # source: analytic - without nearest, all three rows come back
    pt = _bias_table(make_frame)
    tab = pt.search_proctab(make_frame(MJD=100.45), target_type="MBIAS",
                            nearest=False)
    assert len(tab) == 3


def test_search_dark_needs_same_ttime(make_frame):
    # source: analytic - raw DARK rows are matched on CCDCFG and TTIME
    pt = _bias_table(make_frame)
    tab = pt.search_proctab(make_frame(TTIME=30.0), target_type="DARK")
    assert list(tab['filename']) == ["d4.fits"]
    tab = pt.search_proctab(make_frame(TTIME=45.0), target_type="DARK")
    assert len(tab) == 0


def test_search_dark_with_group(make_frame):
    # source: analytic - target_group further restricts DARK rows on GRPID
    pt = _bias_table(make_frame)
    tab = pt.search_proctab(make_frame(TTIME=30.0), target_type="DARK",
                            target_group="g1")
    assert list(tab['filename']) == ["d4.fits"]
    tab = pt.search_proctab(make_frame(TTIME=30.0), target_type="DARK",
                            target_group="g2")
    assert len(tab) == 0


def test_search_mdark_ignores_ttime(make_frame):
    # source: analytic - MDARK rows only need the same CCDCFG (scaled later)
    pt = _bias_table(make_frame)
    tab = pt.search_proctab(make_frame(TTIME=999.0), target_type="MDARK")
    assert list(tab['filename']) == ["md6.fits"]


def test_search_object_by_group(make_frame):
    # source: analytic - OBJECT rows are matched on GRPID only
    pt = _bias_table(make_frame)
    f = make_frame()
    assert list(pt.search_proctab(f, target_type="OBJECT",
                                  target_group="g1")['filename']) == ["o7.fits"]
    assert list(pt.search_proctab(f, target_type="OBJECT",
                                  target_group="g2")['filename']) == ["o8.fits"]


def test_search_other_types_match_stateid(make_frame):
    # source: analytic - any other type (here MARC) is matched on STATEID
    pt = _bias_table(make_frame)
    hit = pt.search_proctab(make_frame(STATEID="otherstate"),
                            target_type="MARC")
    miss = pt.search_proctab(make_frame(), target_type="MARC")
    assert list(hit['filename']) == ["arc9.fits"]
    assert len(miss) == 0


def test_search_without_target_type_returns_none(make_frame):
    # source: analytic - no target type: warning and None
    pt = _bias_table(make_frame)
    assert pt.search_proctab(make_frame()) is None


def test_search_with_no_table_returns_none(make_frame):
    # source: analytic - a Proctab that was never created returns None
    pt = Proctab()
    assert pt.search_proctab(make_frame(), target_type="MBIAS") is None


@pytest.mark.xfail(strict=True, reason="B17: search_proctab raises "
                   "IndexError when the table holds rows of both cameras "
                   "(kcwi_proctab.py:154 indexes the camera-filtered table "
                   "with a mask built from the full table)")
def test_search_with_both_cameras_in_table(make_frame):
    # source: analytic - a BLUE frame searching a table that holds one RED and
    # two BLUE biases of the same configuration must find the two BLUE ones.
    pt = _new()
    _add(pt, make_frame("BLUE", FRAMENO=1, MJD=100.1), "MBIAS", "b1.fits")
    _add(pt, make_frame("RED", FRAMENO=2, MJD=100.2), "MBIAS", "r2.fits")
    _add(pt, make_frame("BLUE", FRAMENO=3, MJD=100.5), "MBIAS", "b3.fits")
    q = make_frame("BLUE", MJD=100.45)
    tab = pt.search_proctab(q, target_type="MBIAS")
    assert list(tab['filename']) == ["b1.fits", "b3.fits"]


def test_write_and_read_back(make_frame, tmp_path):
    # source: analytic - what was put in must come back: integers, strings
    # and floats compare exactly (the fixed-width text keeps full precision).
    pt = _bias_table(make_frame)
    path = tmp_path / "kcwi.proc"
    pt.write_proctab(tfil=str(path))
    assert path.exists()
    rd = Proctab()
    rd.read_proctab(tfil=str(path))
    assert rd.proctab.colnames == pt.proctab.colnames
    assert len(rd.proctab) == len(pt.proctab)
    for name in ('FRAMENO', 'DID', 'STAGE', 'TTIME', 'GANG', 'CWAVE', 'MJD'):
        assert np.array_equal(np.asarray(rd.proctab[name]),
                              np.asarray(pt.proctab[name])), name
    for name in ('CID', 'TYPE', 'GRPID', 'CAM', 'IFU', 'GRAT', 'BIN', 'FILT',
                 'SUFF', 'OFNAME', 'TARGNAME', 'filename'):
        assert list(rd.proctab[name]) == list(pt.proctab[name]), name


def test_write_twice_overwrites(make_frame, tmp_path):
    # source: analytic - a second write to the same file replaces it
    pt = _bias_table(make_frame)
    path = str(tmp_path / "kcwi.proc")
    pt.write_proctab(tfil=path)
    _add(pt, make_frame(FRAMENO=20, MJD=102.0), "MBIAS", "b20.fits")
    pt.write_proctab(tfil=path)
    rd = Proctab()
    rd.read_proctab(tfil=path)
    assert len(rd.proctab) == len(pt.proctab)
    assert "b20.fits" in list(rd.proctab['filename'])


def test_read_missing_file_gives_empty_table(tmp_path):
    # source: analytic - a missing file yields a new, empty proc table
    pt = Proctab()
    pt.read_proctab(tfil=str(tmp_path / "nope.proc"))
    assert len(pt.proctab) == 0
    assert 'FRAMENO' in pt.proctab.colnames


def test_read_sets_display_formats(make_frame, tmp_path):
    # source: analytic - read_proctab formats GANG '7.2f', CWAVE '8.2f',
    # MJD '15.6f' (display only; values are unchanged, see round-trip test)
    pt = _bias_table(make_frame)
    path = str(tmp_path / "kcwi.proc")
    pt.write_proctab(tfil=path)
    rd = Proctab()
    rd.read_proctab(tfil=path)
    assert rd.proctab['GANG'].format == '7.2f'
    assert rd.proctab['CWAVE'].format == '8.2f'
    assert rd.proctab['MJD'].format == '15.6f'


def test_write_with_no_table_writes_nothing(tmp_path):
    # source: analytic - an uncreated proctab writes no file
    pt = Proctab()
    path = tmp_path / "kcwi.proc"
    pt.write_proctab(tfil=str(path))
    assert not path.exists()
