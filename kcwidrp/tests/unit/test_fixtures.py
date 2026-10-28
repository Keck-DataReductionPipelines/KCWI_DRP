"""Tests of the shared fixtures in conftest.py (infrastructure tests, T.1)."""
import numpy as np
import pytest
from astropy.utils import iers
from keckdrpframework.primitives.base_primitive import BasePrimitive

from kcwidrp.primitives.SubtractOverscan import SubtractOverscan
from kcwidrp.primitives.kcwi_file_primitives import (
    KCCDData, kcwi_fits_reader, kcwi_fits_writer, parse_imsec,
    red_amp_gain)
from kcwidrp.tests.conftest import _amp_sections


class _ReturnArgs(BasePrimitive):
    """Hands back the args the harness built through ingest_file.map_ccd."""

    def _perform(self):
        return self.action.args


def test_iers_auto_download_is_off():
    # source: analytic -- D18: conftest switches it off for the session.
    assert iers.conf.auto_download is False


def test_kcwi_config_defaults(kcwi_config):
    config = kcwi_config()
    # source: analytic -- the fixture's documented defaults.
    assert config.plot_level == 0
    assert config.plot_pause == 0
    # source: characterization @ b685577 -- values in configs/kcwi.cfg.
    assert config.output_directory == "redux"
    assert config.NBARS == 120
    assert config.enable_bokeh is True


def test_kcwi_config_overrides(kcwi_config):
    config = kcwi_config(plot_level=3, NBARS=7, new_key="x")
    assert config.plot_level == 3
    assert config.plot_pause == 0
    assert config.NBARS == 7
    assert config.new_key == "x"
    # overrides do not leak into the next call
    assert kcwi_config().plot_level == 0


def test_synthetic_ccddata_blue_defaults(synthetic_ccddata):
    ccd = synthetic_ccddata()
    assert isinstance(ccd, KCCDData)
    # source: analytic -- default frame: two amps of (29 pre-scan + 32 data +
    # 97 overscan) = 158 columns; 2 x 64 data rows + 8 gap rows.
    assert ccd.data.shape == (136, 316)
    assert ccd.data.dtype == np.float64
    assert ccd.uncertainty.array.shape == ccd.data.shape
    assert ccd.uncertainty.array.dtype == np.float64
    assert ccd.mask.shape == ccd.data.shape
    # source: characterization @ b685577 -- astropy 6.1.7 CCDData stores the
    # mask as bool even when uint8 is assigned.
    assert ccd.mask.dtype == np.bool_
    assert not ccd.mask.any()
    assert ccd.flags.shape == ccd.data.shape
    assert ccd.flags.dtype == np.uint8
    assert ccd.noskysub is None
    hdr = ccd.header
    assert hdr["CAMERA"] == "BLUE"
    assert hdr["BINNING"] == "2,2"
    assert hdr["CCDSUM"] == "2 2"
    assert hdr["NVIDINP"] == 4
    # source: analytic -- CCDCFG = CCDSUM(no blank) + CCDMODE(1) + GAINMUL(2)
    # + AMPMNUM(2), as kcwi_fits_reader builds it; CCDMODE 1 and GAINMUL 10
    # as in kb260115_00011, whose CCDCFG is "1111000" at 1x1 binning.
    assert hdr["CCDCFG"] == "2211000"
    assert hdr["NAXIS1"] == 316 and hdr["NAXIS2"] == 136


def test_synthetic_ccddata_has_every_scanned_keyword(synthetic_ccddata,
                                                     header_keyword_info):
    info = header_keyword_info
    hdr = synthetic_ccddata().header
    expected = info.scanned - info.structural
    # source: analytic -- every keyword in the AST scan of conftest.py is
    # present in the default header.
    assert sorted(expected - set(hdr.keys())) == []
    # and every default is a keyword that was scanned or is derived
    # (BARTANG/RARTANG, per-amp keys are added outside the defaults table).
    assert set(info.defaults) <= info.scanned


@pytest.mark.parametrize("camera, ampmode, binning, amps", [
    ("BLUE", "ALL", (2, 2), [1, 2, 3, 4]),
    ("BLUE", "TUP", (2, 2), [1, 2]),
    ("BLUE", "ALL", (1, 1), [1, 2, 3, 4]),
    ("RED", "L2U2L1U1", (2, 2), [0, 1, 2, 3]),
    ("RED", "L2U2", (2, 2), [1, 3]),
    ("RED", "L2U2L1U1", (1, 1), [0, 1, 2, 3]),
])
def test_amp_sections_are_disjoint_and_inside_the_frame(
        synthetic_ccddata, camera, ampmode, binning, amps):
    ccd = synthetic_ccddata(camera=camera, ampmode=ampmode, binning=binning)
    ny, nx = ccd.data.shape
    covered = np.zeros(ccd.data.shape, dtype=int)
    for ia in amps:
        for kind in ("BSEC", "DSEC"):
            (y0, y1, x0, x1), _ = parse_imsec(ccd.header["%s%d" % (kind, ia)])
            assert 0 <= y0 <= y1 < ny and 0 <= x0 <= x1 < nx
            covered[y0:y1 + 1, x0:x1 + 1] += 1
    # source: analytic -- the data and overscan sections of all amps never
    # overlap (each pixel is in at most one section).
    assert covered.max() == 1
    # source: analytic -- the amps present are exactly those with sections.
    present = [k for k in range(5) if "DSEC%d" % k in ccd.header]
    assert present == amps


# Section keywords copied from real frames of 2026-01-15 (raw headers; the
# fixture reproduces them when given the real data size).
_REAL_SECTIONS = [
    # (id, camera, ampmode, binning, nd, nyd, prescan, overscan, ygap,
    #  naxis (nx, ny), {keyword: section})
    ("kb260115_00010", "BLUE", "ALL", (1, 1), 2048, 2056, 54, 198, 288,
     (4600, 4400),
     {"ASEC1": "[1:2048,1:2056]", "BSEC1": "[2103:2300,1:2056]",
      "CSEC1": "[1:2048,1:2056]", "DSEC1": "[55:2102,1:2056]",
      "BSEC2": "[2498:2301,1:2056]", "CSEC2": "[4096:2049,1:2056]",
      "DSEC2": "[4546:2499,1:2056]",
      "BSEC3": "[2103:2300,4400:2345]", "CSEC3": "[1:2048,4112:2057]",
      "DSEC3": "[55:2102,4400:2345]", "CSEC4": "[4096:2049,4112:2057]",
      "DSEC4": "[4546:2499,4400:2345]"}),
    ("kb260115_00001", "BLUE", "TUP", (2, 2), 1024, 2056, 29, 97, 144,
     (2300, 2200),
     {"ASEC1": "[1:2048,1:4112]", "BSEC1": "[1054:1150,1:2056]",
      "CSEC1": "[1:2048,4112:1]", "DSEC1": "[30:1053,1:2056]",
      "BSEC2": "[1247:1151,1:2056]", "CSEC2": "[4096:2049,4112:1]",
      "DSEC2": "[2271:1248,1:2056]"}),
    ("kr260115_00010", "RED", "L2U2L1U1", (1, 1), 2057, 2064, 7, 53, 120,
     (4234, 4248),
     {"ASEC0": "[1:2057,2064:1:-1]", "BSEC0": "[2065:2117,4248:2185:-1]",
      "CSEC0": "[1:2057,2065:4128:-1]", "DSEC0": "[8:2064,2185:4248]",
      "TSEC0": "[1:2117,4248:2125:-1]",
      "ASEC1": "[1:2057,1:2064]", "BSEC1": "[2065:2117,1:2064]",
      "CSEC1": "[1:2057,1:2064]", "DSEC1": "[8:2064,1:2064]",
      "TSEC1": "[1:2117,1:2124]",
      "ASEC2": "[2057:1:-1,2064:1:-1]", "BSEC2": "[2170:2118:-1,4248:2185:-1]",
      "CSEC2": "[2058:4114:-1,2065:4128:-1]", "DSEC2": "[2171:4227,2185:4248]",
      "TSEC2": "[4234:2118:-1,4248:2125:-1]",
      "ASEC3": "[2057:1:-1,1:2064]", "BSEC3": "[2170:2118:-1,1:2064]",
      "CSEC3": "[2058:4114:-1,1:2064]", "DSEC3": "[2171:4227,1:2064]",
      "TSEC3": "[4234:2118:-1,1:2124]"}),
    ("kr260115_00001", "RED", "L2U2", (2, 2), 1028, 2064, 4, 72, 60,
     (2208, 2124),
     {"ASEC1": "[2:2057:2,1:4128:2]", "BSEC1": "[1033:1104,1:2064]",
      "CSEC1": "[2:2057:2,1:4128:2]", "DSEC1": "[5:1032,1:2064]",
      "TSEC1": "[1:1104,1:2124]",
      "ASEC3": "[2057:2:-2,1:4128:2]", "BSEC3": "[1176:1105:-1,1:2064]",
      "CSEC3": "[2058:4113:-2,1:4128:2]", "DSEC3": "[1177:2204,1:2064]",
      "TSEC3": "[2208:1105:-1,1:2124]"}),
]


@pytest.mark.parametrize("case", _REAL_SECTIONS, ids=[c[0] for c in
                                                       _REAL_SECTIONS])
def test_amp_sections_reproduce_the_real_frames(case):
    (_, camera, ampmode, binning, nd, nyd, prescan, overscan, ygap,
     (nx, ny), expected) = case
    got_ny, got_nx, sections = _amp_sections(
        camera, ampmode, binning, nd, nyd, prescan, overscan, ygap)
    # source: the raw headers of the named frame (read with astropy,
    # 2026-01-15); the pre-scan, overscan and gap widths are those frames'.
    assert (got_nx, got_ny) == (nx, ny)
    for key, value in expected.items():
        assert sections[key] == value, key


def test_blue_has_one_based_amps_and_ampids(synthetic_ccddata):
    hdr = synthetic_ccddata().header
    # source: kb260115_00011 (post-ingest) -- ALL: AMPID1..4 = 0..3.
    assert [hdr["AMPID%d" % i] for i in range(1, 5)] == [0, 1, 2, 3]
    assert "ASEC0" not in hdr and "GAIN0" not in hdr
    tup = synthetic_ccddata(ampmode="TUP").header
    # source: kb260115_00003 (post-ingest) -- TUP: two amps, AMPID 2 and 3,
    # GAIN 0.158 each; AMPMNUM 9, CCDMODE 0.
    assert (tup["AMPID1"], tup["AMPID2"]) == (2, 3)
    assert (tup["GAIN1"], tup["GAIN2"]) == (0.158, 0.158)
    assert "ASEC3" not in tup and tup["NVIDINP"] == 2
    assert tup["AMPMNUM"] == 9 and tup["CCDMODE"] == 0


def test_red_gain_keys_are_numbered_from_zero(synthetic_ccddata):
    hdr = synthetic_ccddata(camera="RED").header
    # source: kr260115_00013 post-ingest GAIN0..3 = red_amp_gain (L1, L2, U1,
    # U2) = 1.54, 1.551, 1.61, 1.526 (also the table in kcwi_file_primitives).
    assert [hdr["GAIN%d" % i] for i in range(4)] == [1.54, 1.551, 1.61, 1.526]
    assert [hdr["GAIN%d" % i] for i in range(4)] == [
        red_amp_gain[a] for a in ("L1", "L2", "U1", "U2")]
    assert "GAIN4" not in hdr and "ASEC4" not in hdr
    assert "AMPID1" not in hdr and "ASEC0" in hdr
    two = synthetic_ccddata(camera="RED", ampmode="L2U2").header
    # source: kr260115_00003 post-ingest -- L2U2: GAIN1 1.551, GAIN3 1.526.
    assert sorted(k for k in two if k.startswith("GAIN") and k[4:].isdigit()
                  ) == ["GAIN1", "GAIN3"]
    assert (two["GAIN1"], two["GAIN3"]) == (1.551, 1.526)
    assert two["NVIDINP"] == 2 and two["AMPMNUM"] == 1


def test_fix_header_leaves_a_default_frame_unchanged(synthetic_ccddata):
    from kcwidrp.primitives.kcwi_file_primitives import fix_header
    for camera in ("BLUE", "RED"):
        ccd = synthetic_ccddata(camera=camera)
        before = {k: ccd.header[k] for k in ccd.header if k != "HISTORY"}
        fix_header(ccd)
        # source: analytic -- the fixture holds post-ingest values, so
        # fix_header (which only recomputes GAINn etc.) is idempotent on it.
        assert {k: ccd.header[k] for k in ccd.header
                if k != "HISTORY"} == before


def test_unsupported_ampmode_raises(synthetic_ccddata):
    with pytest.raises(ValueError, match="unsupported camera/ampmode"):
        synthetic_ccddata(camera="BLUE", ampmode="L1")


@pytest.mark.parametrize("camera, n_amps, aoff", [
    ("BLUE", 4, 1), ("RED", 4, 0)])
def test_default_frames_build_args_through_map_ccd(
        run_primitive, synthetic_ccddata, camera, n_amps, aoff):
    run_primitive(_ReturnArgs, synthetic_ccddata(camera=camera))
    args = run_primitive.last.action.args
    bsec, dsec, tsec, direc, amps, off = args.map_ccd
    # source: analytic -- four amps; BLUE numbers them from 1, RED from 0
    # (ingest_file.map_ccd).
    assert len(amps) == n_amps
    assert amps == list(range(aoff, aoff + n_amps))
    assert off == aoff
    # source: analytic -- the section lists hold a slot per amp number
    assert len(bsec) == len(dsec) == len(tsec) == len(direc) == 4


def test_red_two_amp_frame_builds_args(run_primitive, synthetic_ccddata):
    run_primitive(_ReturnArgs, synthetic_ccddata(camera="RED",
                                                 ampmode="L2U2"))
    _, _, _, _, amps, off = run_primitive.last.action.args.map_ccd
    # source: analytic -- L2U2 = red_amp_dict entries 1 and 3.
    assert amps == [1, 3] and off == 0


@pytest.mark.parametrize("camera, ampmode", [
    ("BLUE", "ALL"), ("BLUE", "TUP"), ("RED", "L2U2L1U1"), ("RED", "L2U2")])
def test_default_frame_round_trips_through_writer_and_reader(
        synthetic_ccddata, tmp_path, camera, ampmode):
    ccd = synthetic_ccddata(camera=camera, ampmode=ampmode)
    ccd.flags[3, 4] = 7
    ccd.mask[5, 6] = True
    # Real raw/ingested headers have no HISTORY card (checked on
    # kb260115_00011 and kr260115_00013), and the writer reads
    # header["HISTORY"] without a default (known bug, strict xfail in
    # unit/primitives/test_kcwi_file_primitives.py).  Every primitive appends
    # a HISTORY line, so frames reaching the writer have one.
    ccd.header.add_history("synthetic frame")
    # kcwi_fits_writer stores float64 as float32; compare with that rounding.
    data = ccd.data.astype(np.float32)
    unc = ccd.uncertainty.array.astype(np.float32)
    mask = ccd.mask.copy()
    flags = ccd.flags.copy()
    kcwi_fits_writer(ccd, output_file="rt.fits", output_dir=str(tmp_path))
    back, table = kcwi_fits_reader(str(tmp_path / "rt.fits"))
    # source: analytic -- a write/read round trip preserves the planes (data
    # and uncertainty at float32 precision, as the writer converts them).
    assert table is None
    assert np.array_equal(back.data.astype(np.float32), data)
    assert np.array_equal(back.uncertainty.array.astype(np.float32), unc)
    assert np.array_equal(back.mask, mask)
    assert np.array_equal(back.flags, flags)
    assert back.header["NAXIS1"] == data.shape[1]
    assert back.header["NAXIS2"] == data.shape[0]
    assert back.header["CCDCFG"] == ccd.header["CCDCFG"]


@pytest.mark.parametrize("camera, ampmode, binning", [
    ("BLUE", "ALL", (2, 2)), ("BLUE", "ALL", (1, 1)), ("BLUE", "TUP", (2, 2)),
    ("RED", "L2U2L1U1", (2, 2)), ("RED", "L2U2L1U1", (1, 1)),
    ("RED", "L2U2", (2, 2))])
def test_subtract_overscan_runs_on_default_frames(
        run_primitive, synthetic_ccddata, camera, ampmode, binning):
    ccd = synthetic_ccddata(camera=camera, ampmode=ampmode, binning=binning)
    out = run_primitive(SubtractOverscan, ccd)
    # source: analytic -- the default overscan is wider than minoscanpix
    # (75 BLUE, 20 RED in kcwi.cfg), so every amp is fitted and OSCANSUB is
    # True.
    assert out.ccddata.header["OSCANSUB"] is True


def test_subtract_overscan_removes_the_overscan_level(run_primitive,
                                                      synthetic_ccddata):
    ccd = synthetic_ccddata(camera="RED", ampmode="L2U2")
    (y0, y1, x0, x1), _ = parse_imsec(ccd.header["BSEC1"])
    ccd.data[y0:y1 + 1, x0:x1 + 1] = 500.0
    (y0, y1, x0, x1), _ = parse_imsec(ccd.header["DSEC1"])
    ccd.data[y0:y1 + 1, x0:x1 + 1] = 520.0
    out = run_primitive(SubtractOverscan, ccd)
    # source: analytic -- constant overscan 500 is subtracted from the data
    # section of amp 1 (520 -> 20); the overscan value is recorded.
    assert np.all(out.ccddata.data[y0:y1 + 1, x0:x1 + 1] == 20.0)
    assert out.ccddata.header["OSCNVAL1"] == 500


def test_synthetic_ccddata_is_deterministic(synthetic_ccddata):
    assert np.array_equal(synthetic_ccddata().data, synthetic_ccddata().data)
    assert not np.array_equal(synthetic_ccddata(seed=1).data,
                              synthetic_ccddata(seed=2).data)


def test_tmp_redux(tmp_redux, tmp_path):
    assert tmp_redux.root == tmp_path
    assert tmp_redux.redux.is_dir()
    assert tmp_redux.redux.parent == tmp_path
    # source: analytic -- a new Proctab has the 19 columns of new_proctab()
    # and no rows.
    table = tmp_redux.proctab.proctab
    assert len(table) == 0
    assert len(table.colnames) == 19
    assert tmp_redux.proctab_file.parent == tmp_redux.redux
