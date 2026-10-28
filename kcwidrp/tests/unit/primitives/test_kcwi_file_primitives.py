"""Characterization tests for the stand-alone functions of
kcwidrp.primitives.kcwi_file_primitives (task T.3): parse_imsec, strip_fname,
get_master_name, get_unique_CCD_master_name, get_unique_STATEID_master_name,
fix_header, plotlabel, kcwi_fits_writer <-> kcwi_fits_reader, write_table <->
read_table.

Expected values are hand-built (inputs written in the test), round trips, or
labelled characterizations of tables that have no external specification.
"""
import logging
from types import SimpleNamespace

import numpy as np
import pytest
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.table import Table
from astropy.time import Time

from kcwidrp.primitives.kcwi_file_primitives import (
    KCCDData, fix_header, get_master_name, get_unique_CCD_master_name,
    get_unique_STATEID_master_name, kcwi_fits_reader, kcwi_fits_writer,
    parse_imsec, plotlabel, read_table, strip_fname, write_table)


# ------------------------------------------------------------- parse_imsec
@pytest.mark.parametrize("section, sec, stride", [
    # source: analytic - 1-based inclusive IRAF section [x0:x1,y0:y1] becomes
    # zero-based (y0, y1, x0, x1) with inclusive end, strides (ystride, xstride)
    ("[1:10,1:20]", (0, 19, 0, 9), (1, 1)),
    ("[5:12,3:8]", (2, 7, 4, 11), (1, 1)),
    # a third element is the stride along that axis
    ("[1:10:2,1:20:3]", (0, 19, 0, 9), (3, 2)),
    # descending axes: limits are swapped and the stride is made negative
    ("[10:1,20:1]", (0, 19, 0, 9), (-1, -1)),
    ("[10:1:2,1:20]", (0, 19, 0, 9), (1, -2)),
    ("[1:10,20:1:4]", (0, 19, 0, 9), (-4, 1)),
    # the lower limit is clamped so that no index is negative
    ("[0:10,0:5]", (0, 4, 0, 9), (1, 1)),
])
def test_parse_imsec(section, sec, stride):
    assert parse_imsec(section) == (sec, stride)


def test_parse_imsec_section_slices_expected_pixels():
    # source: analytic - section [3:6,2:4] selects columns 3..6 and rows 2..4
    # (1-based, inclusive): numpy slice [1:4, 2:6]
    sec, _ = parse_imsec("[3:6,2:4]")
    img = np.arange(100).reshape(10, 10)
    cut = img[sec[0]:sec[1] + 1, sec[2]:sec[3] + 1]
    assert cut.shape == (3, 4)
    assert np.array_equal(cut, img[1:4, 2:6])


# ------------------------------------------------------------- name helpers
@pytest.mark.parametrize("name, stem", [
    # source: analytic - pathlib stem: strip directory and last extension
    ("kb260101_00001.fits", "kb260101_00001"),
    ("/data/redux/kb260101_00001_intf.fits", "kb260101_00001_intf"),
    ("relative/dir/kr1.fits", "kr1"),
    ("a.b.fits", "a.b"),
    ("noextension", "noextension"),
])
def test_strip_fname(name, stem):
    assert strip_fname(name) == stem


def test_strip_fname_accepts_path_object(tmp_path):
    # source: analytic - a pathlib.Path is accepted
    assert strip_fname(tmp_path / "kb1.fits") == "kb1"


@pytest.mark.parametrize("empty", ["", None])
def test_strip_fname_empty_returned_unchanged(empty):
    # source: analytic - documented: "or filename if error occurs"
    assert strip_fname(empty) == empty


def test_get_master_name():
    # source: analytic - name = stem(filename) + '_' + lower(type) + '.fits',
    # directory dropped
    tab = {"filename": ["/x/y/kb260101_00005.fits", "kb260101_00006_int.fits"]}
    assert get_master_name(tab, "MBIAS") == "kb260101_00005_mbias.fits"
    assert get_master_name(tab, "MFLAT", loc=1) == \
        "kb260101_00006_int_mflat.fits"


def test_get_unique_ccd_master_name(synthetic_ccddata):
    # source: analytic - lower(type) + '_' + CCDCFG + '.fits'; the fixture's
    # CCDCFG is CCDSUM '22' + CCDMODE '1' + GAINMUL '10' + AMPMNUM '00'
    ccd = synthetic_ccddata()
    assert get_unique_CCD_master_name(ccd) == "mbias_2211000.fits"
    assert get_unique_CCD_master_name(ccd, "MDARK") == "mdark_2211000.fits"


def test_get_unique_stateid_master_name(synthetic_ccddata):
    # source: analytic - lower(type) + '_' + STATEID + '.fits'
    ccd = synthetic_ccddata(STATEID="abc123")
    assert get_unique_STATEID_master_name(ccd, "MFLAT") == "mflat_abc123.fits"


# --------------------------------------------------------------- fix_header
def _blue(synthetic_ccddata, gainmul, ampids=(0, 1, 2, 3), **kw):
    ccd = synthetic_ccddata(camera="BLUE", GAINMUL=gainmul, **kw)
    for i, a in enumerate(ampids, start=1):
        ccd.header["AMPID%d" % i] = a
    for i in range(1, 5):                      # remove the fixture's gains
        ccd.header["GAIN%d" % i] = -1.0
    return ccd


@pytest.mark.parametrize("gainmul, gains", [
    # source: characterization @ b0e2550 - blue_amp_gain table (slow readout);
    # GAINMUL 10 equals the GAINn of the real frame kb260115_00011 recorded in
    # conftest._GAINS (0.157, 0.160, 0.158, 0.158)
    (1, [1.570, 1.600, 1.610, 1.600]),
    (2, [0.785, 0.800, 0.805, 0.800]),
    (5, [0.314, 0.320, 0.325, 0.319]),
    (10, [0.157, 0.160, 0.158, 0.158]),
])
def test_fix_header_blue_sets_gain_per_amp(synthetic_ccddata, gainmul, gains):
    ccd = _blue(synthetic_ccddata, gainmul)
    fix_header(ccd)
    assert [ccd.header["GAIN%d" % i] for i in range(1, 5)] == gains


def test_fix_header_blue_gain_follows_ampid(synthetic_ccddata):
    # source: characterization @ b0e2550 - GAINn is looked up with AMPIDn as
    # the key: amp 1 reports AMPID 3, so gets the AMPID-3 gain of GAINMUL 5.
    ccd = _blue(synthetic_ccddata, 5, ampids=(3, 2, 1, 0))
    fix_header(ccd)
    assert [ccd.header["GAIN%d" % i] for i in range(1, 5)] == \
        [0.319, 0.325, 0.320, 0.314]


RED_AMPS = [
    # (AMPMODE, expected AMPMNUM, expected {GAINk: value}) with k the red amp
    # index L1=0, L2=1, U1=2, U2=3
    # source: characterization @ b0e2550 - AMPMODE -> AMPMNUM table and the
    # red_amp_gain values (slow readout, high gain)
    ("L2U2L1U1", 0, {0: 1.54, 1: 1.551, 2: 1.61, 3: 1.526}),
    ("L2U2", 1, {1: 1.551, 3: 1.526}),
    ("L1U1", 2, {0: 1.54, 2: 1.61}),
    ("L2L1", 3, {0: 1.54, 1: 1.551}),
    ("U2U1", 4, {2: 1.61, 3: 1.526}),
    ("L2", 5, {1: 1.551}),
    ("U2", 6, {3: 1.526}),
    ("L1", 7, {0: 1.54}),
    ("U1", 8, {2: 1.61}),
    ("XX", 0, {}),
]


@pytest.mark.parametrize("ampmode, ampnum, gains", RED_AMPS)
def test_fix_header_red(synthetic_ccddata, ampmode, ampnum, gains):
    ccd = synthetic_ccddata(camera="RED", AMPMODE=ampmode, TAPLINES=3,
                            CDSSPEED=1)
    for k in range(4):
        ccd.header["GAIN%d" % k] = -1.0
    fix_header(ccd)
    assert ccd.header["AMPMNUM"] == ampnum
    assert ccd.header["NVIDINP"] == 3          # copied from TAPLINES
    assert ccd.header["GAINMUL"] == 1
    assert ccd.header["CCDMODE"] == 1          # copied from CDSSPEED
    for k in range(4):
        assert ccd.header["GAIN%d" % k] == gains.get(k, -1.0)


def test_fix_header_red_ccdmode_defaults_to_zero(synthetic_ccddata):
    # source: characterization @ b0e2550 - no CDSSPEED keyword -> CCDMODE 0
    ccd = synthetic_ccddata(camera="RED", ampmode="L2U2L1U1")
    del ccd.header["CDSSPEED"]
    fix_header(ccd)
    assert ccd.header["CCDMODE"] == 0


def test_fix_header_red_ait_frame_gets_mjd_and_targname(synthetic_ccddata):
    # source: independent - astropy.time for the MJD of DATE-END; frames
    # without TELESCOP (Caltech AIT) get TARGNAME from OBJECT and MJD from
    # DATE-END.
    ccd = synthetic_ccddata(camera="RED", OBJECT="AIT test",
                            **{"DATE-END": "2021-03-04T06:07:08"})
    del ccd.header["TELESCOP"]
    fix_header(ccd)
    assert ccd.header["TARGNAME"] == "AIT test"
    # Tolerance: MJD ~6e4 days in float64 has resolution 7e-12 d; the code
    # sums hour/24 + minute/1440 + second/86400, accurate to ~1e-11 d; 1e-9 d
    # (86 microseconds).
    assert ccd.header["MJD"] == pytest.approx(
        Time("2021-03-04T06:07:08").mjd, rel=0, abs=1e-9)


def test_fix_header_unknown_camera_prints_error(synthetic_ccddata, capsys):
    # source: characterization @ b0e2550 - an unknown CAMERA prints an error
    # and leaves the header unchanged
    ccd = synthetic_ccddata()
    ccd.header["CAMERA"] = "GREEN"
    before = dict(ccd.header)
    fix_header(ccd)
    assert "illegal CAMERA keyword: GREEN" in capsys.readouterr().out
    assert dict(ccd.header) == before


# ---------------------------------------------------------------- plotlabel
def _args(camera, cwave=4599.83, **kw):
    hdr = fits.Header({"FRAMENO": 12, "CAMERA": camera})
    d = dict(ccddata=SimpleNamespace(header=hdr), illum="Dome",
             ifuname="Small", grating="BL", cwave=cwave, filter="None")
    d.update(kw)
    return SimpleNamespace(**d)


def test_plotlabel_blue_includes_filter():
    # source: analytic - hand-assembled from the documented fields: frame
    # number, illumination, slicer, grating, central wavelength (truncated
    # to int by %d of int()), filter for the blue channel
    assert plotlabel(_args("BLUE")) == \
        "[ Img # 12 (Dome) Slicer: Small, Grating: BL, CWave: 4599, " \
        "Filter: None ] "


def test_plotlabel_red_has_no_filter():
    # source: analytic - the red channel has no filter field
    assert plotlabel(_args("RED", cwave=7200.4, grating="RM1")) == \
        "[ Img # 12 (Dome) Slicer: Small, Grating: RM1, CWave: 7200] "


# ----------------------------------------------- kcwi_fits_writer / reader
@pytest.fixture
def frame(synthetic_ccddata):
    """A blue frame whose header is acceptable to the writer."""
    ccd = synthetic_ccddata(AMPID1=0, AMPID2=1, AMPID3=2, AMPID4=3)
    # The fixture stores NAXIS1/NAXIS2 as ordinary cards at the end of the
    # header, which astropy refuses to write; real headers come from the file.
    del ccd.header["NAXIS1"], ccd.header["NAXIS2"]
    ccd.header.add_history("a history line")
    ccd.mask[3, 4] = True
    ccd.flags[5, 6] = 128
    ccd.flags[1, 1] = 4
    return ccd


def test_writer_reader_roundtrip_all_planes(frame, tmp_path):
    # source: analytic - what was written must be read back.  Data and
    # uncertainty are written as float32, so the expected values are the
    # float32 rounding of the originals, read back as float64 (data) or float32
    # (uncertainty); mask is boolean; flags are uint8 and exact.
    data = frame.data.copy()
    unc = frame.uncertainty.array.copy()
    noskysub = data * 2.0
    frame.noskysub = noskysub
    kcwi_fits_writer(frame, output_file="rt.fits", output_dir=str(tmp_path))
    back, table = kcwi_fits_reader(str(tmp_path / "rt.fits"))

    assert isinstance(back, KCCDData)
    assert back.data.dtype == np.float64
    assert np.array_equal(back.data, data.astype(np.float32))
    assert isinstance(back.uncertainty, StdDevUncertainty)
    assert np.array_equal(back.uncertainty.array, unc.astype(np.float32))
    assert back.mask.dtype == bool
    assert np.array_equal(back.mask, frame.mask)
    assert back.flags.dtype == np.uint8
    assert np.array_equal(back.flags, frame.flags)
    assert back.flags[5, 6] == 128 and back.flags[1, 1] == 4
    assert np.array_equal(back.noskysub, noskysub.astype(np.float32))
    assert table is None


def test_writer_reader_roundtrip_header(frame, tmp_path):
    # source: analytic - every card present before writing has the same value
    # after the round trip; the reader sets the unit from BUNIT.
    before = {k: frame.header[k] for k in frame.header
              if k not in ("HISTORY", "COMMENT")}
    kcwi_fits_writer(frame, output_file="rt.fits", output_dir=str(tmp_path))
    back, _ = kcwi_fits_reader(str(tmp_path / "rt.fits"))
    for key, value in before.items():
        assert back.header[key] == value, key
    assert back.unit == "adu"
    assert back.uncertainty.unit == "adu"
    assert "a history line" in [str(h) for h in back.header["HISTORY"]]


def test_writer_adds_version_history_once(frame, tmp_path):
    # source: analytic - 'kcwidrp version=...' is added when absent and not
    # repeated by a second write of the same object
    kcwi_fits_writer(frame, output_file="a.fits", output_dir=str(tmp_path))
    kcwi_fits_writer(frame, output_file="b.fits", output_dir=str(tmp_path))
    versions = [str(h) for h in frame.header["HISTORY"]
                if "kcwidrp version" in str(h)]
    assert len(versions) == 1


def test_writer_converts_float64_to_float32(frame, tmp_path):
    # source: documented in the writer docstring ("Converts float64 data to
    # float32"): the file holds float32 for data and uncertainty
    kcwi_fits_writer(frame, output_file="f.fits", output_dir=str(tmp_path))
    with fits.open(tmp_path / "f.fits") as hdul:
        assert hdul["PRIMARY"].data.dtype == np.dtype(">f4")
        assert hdul["UNCERT"].data.dtype == np.dtype(">f4")
        assert hdul["FLAGS"].data.dtype == np.dtype("uint8")
        assert hdul["MASK"].data.dtype == np.dtype("uint8")


def test_writer_suffix_and_directory_handling(frame, tmp_path):
    # source: analytic - output_file's directory part is dropped; suffix is
    # inserted before the extension as '_<suffix>'
    kcwi_fits_writer(frame, output_file="/elsewhere/kb1.fits",
                     output_dir=str(tmp_path), suffix="icube")
    assert (tmp_path / "kb1_icube.fits").exists()


def test_writer_overwrites_existing_file(frame, tmp_path):
    # source: analytic - writeto(overwrite=True): a second write succeeds
    kcwi_fits_writer(frame, output_file="o.fits", output_dir=str(tmp_path))
    frame.data[:] = 1.0
    kcwi_fits_writer(frame, output_file="o.fits", output_dir=str(tmp_path))
    back, _ = kcwi_fits_reader(str(tmp_path / "o.fits"))
    assert np.all(back.data == 1.0)


def test_writer_does_not_write_the_table(frame, tmp_path):
    # source: documented - the writer docstring says "Currently fits tables
    # are not written out"; the table argument is ignored.
    tab = fits.BinTableHDU.from_columns(
        [fits.Column(name="t", format="D", array=np.arange(3.0))],
        name="Exposure Events")
    kcwi_fits_writer(frame, table=tab, output_file="t.fits",
                     output_dir=str(tmp_path))
    back, table = kcwi_fits_reader(str(tmp_path / "t.fits"))
    assert table is None


@pytest.mark.xfail(strict=True, reason="B20: kcwi_fits_writer raises "
                   "KeyError for a header without any HISTORY card "
                   "(kcwi_file_primitives.py:1197 ccddata.header['HISTORY'])")
def test_writer_accepts_header_without_history(frame, tmp_path):
    # source: analytic - a frame whose header has no HISTORY card yet is
    # still writable: the writer adds the version history itself.
    del frame.header["HISTORY"]
    kcwi_fits_writer(frame, output_file="nh.fits", output_dir=str(tmp_path))
    assert (tmp_path / "nh.fits").exists()


def test_reader_returns_exposure_events_table(frame, tmp_path):
    # source: analytic - a hand-built file with an 'Exposure Events' table
    # extension gives back that table as the second return value.
    kcwi_fits_writer(frame, output_file="base.fits", output_dir=str(tmp_path))
    with fits.open(tmp_path / "base.fits") as hdul:
        events = fits.BinTableHDU.from_columns(
            [fits.Column(name="t", format="D", array=np.array([1.0, 2.5]))],
            name="Exposure Events")
        out = fits.HDUList(list(hdul) + [events])
        out.writeto(tmp_path / "with_events.fits")
    back, table = kcwi_fits_reader(str(tmp_path / "with_events.fits"))
    assert table is not None
    assert list(table.data["t"]) == [1.0, 2.5]
    assert back.data.shape == frame.data.shape


def test_reader_builds_ccdcfg_when_missing(frame, tmp_path):
    # source: analytic - CCDCFG = CCDSUM without space + CCDMODE (1 digit) +
    # GAINMUL (2 digits) + AMPMNUM (2 digits): '22' '1' '10' '00'
    del frame.header["CCDCFG"]
    kcwi_fits_writer(frame, output_file="c.fits", output_dir=str(tmp_path))
    back, _ = kcwi_fits_reader(str(tmp_path / "c.fits"))
    assert back.header["CCDCFG"] == "2211000"


def test_reader_missing_file_raises(tmp_path, capsys):
    # source: documented - FileNotFoundError for a missing file
    with pytest.raises(FileNotFoundError):
        kcwi_fits_reader(str(tmp_path / "nope.fits"))


def test_reader_applies_fix_header(frame, tmp_path):
    # source: analytic - the reader runs fix_header, so a stale GAIN1 written
    # in the file is replaced by the table value for AMPID1 = 0, GAINMUL = 10
    frame.header["GAIN1"] = -5.0
    kcwi_fits_writer(frame, output_file="g.fits", output_dir=str(tmp_path))
    back, _ = kcwi_fits_reader(str(tmp_path / "g.fits"))
    assert back.header["GAIN1"] == 0.157   # characterization @ b0e2550


# ------------------------------------------------- write_table / read_table
def test_write_table_read_table_roundtrip(tmp_path):
    # source: analytic - columns, names, values, comment and keywords come
    # back as written.
    cols = [np.arange(5, dtype=np.int32), np.linspace(0.0, 1.0, 5),
            np.array(["a", "bb", "ccc", "d", "e"])]
    write_table(output_dir=str(tmp_path), table=cols,
                names=["i", "x", "s"], comment="hello table",
                keywords={"MYKEY": 42}, output_name="t.fits")
    tab = read_table(input_dir=str(tmp_path), file_name="t.fits")
    assert isinstance(tab, Table)
    assert tab.colnames == ["i", "x", "s"]
    assert np.array_equal(tab["i"], cols[0])
    assert tab["i"].dtype.kind == "i" and tab["i"].dtype.itemsize == 4  # FITS is big-endian
    assert np.array_equal(tab["x"], cols[1])
    assert list(tab["s"]) == ["a", "bb", "ccc", "d", "e"]
    assert tab.meta["MYKEY"] == 42
    assert tab.meta["comments"] == ["hello table"]   # astropy's name for COMMENT


def test_write_table_does_not_clobber_by_default(tmp_path):
    # source: analytic - documented: clobber=False leaves an existing file
    write_table(str(tmp_path), [np.arange(3)], ["a"], output_name="t.fits")
    write_table(str(tmp_path), [np.arange(10, 13)], ["a"],
                output_name="t.fits")
    tab = read_table(str(tmp_path), "t.fits")
    assert list(tab["a"]) == [0, 1, 2]


def test_write_table_clobber_true_replaces(tmp_path):
    # source: analytic - documented: clobber=True overwrites
    write_table(str(tmp_path), [np.arange(3)], ["a"], output_name="t.fits")
    write_table(str(tmp_path), [np.arange(10, 13)], ["a"],
                output_name="t.fits", clobber=True)
    tab = read_table(str(tmp_path), "t.fits")
    assert list(tab["a"]) == [10, 11, 12]


def test_read_table_missing_returns_none(tmp_path):
    # source: documented - "None if unsuccessful"
    assert read_table(str(tmp_path), "nope.fits") is None
