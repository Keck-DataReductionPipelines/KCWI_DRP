"""Tests of the shared fixtures in conftest.py (infrastructure tests, T.1)."""
import numpy as np
from astropy.utils import iers

from kcwidrp.primitives.kcwi_file_primitives import KCCDData, parse_imsec


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
    # source: analytic -- default shape in the fixture
    assert ccd.data.shape == (128, 144)
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
    # + AMPMNUM(2), as kcwi_fits_reader builds it.
    assert hdr["CCDCFG"] == "2200100"
    assert hdr["NAXIS1"] == 144 and hdr["NAXIS2"] == 128


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


def test_synthetic_ccddata_amp_sections_tile_the_frame(synthetic_ccddata):
    ccd = synthetic_ccddata()
    covered = np.zeros(ccd.data.shape, dtype=int)
    for ia in range(1, 5):
        (y0, y1, x0, x1), _ = parse_imsec(ccd.header["ASEC%d" % ia])
        covered[y0:y1 + 1, x0:x1 + 1] += 1
    # source: analytic -- the four amp sections partition the raw frame.
    assert np.array_equal(covered, np.ones_like(covered))


def test_synthetic_ccddata_overrides(synthetic_ccddata):
    ccd = synthetic_ccddata(camera="RED", binning=(1, 2), ampmode="L1",
                            shape=(40, 50), OBJECT="Foo", EXPTIME=12.5)
    assert ccd.data.shape == (40, 50)
    assert ccd.header["CAMERA"] == "RED"
    assert ccd.header["BINNING"] == "1,2"
    assert ccd.header["AMPMODE"] == "L1"
    assert ccd.header["NVIDINP"] == 1
    assert ccd.header["OBJECT"] == "Foo"
    assert ccd.header["EXPTIME"] == 12.5
    assert "RARTANG" in ccd.header and "BARTANG" not in ccd.header


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
