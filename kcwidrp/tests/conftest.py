"""Shared fixtures for the KCWI DRP test suite (task T.1; rules in docs/testing.md).

Fixtures
--------
kcwi_config(**overrides)
    Factory: the parsed ``configs/kcwi.cfg``, loaded the way
    ``scripts/reduce_kcwi.py`` loads it, with plot_level = 0 and
    plot_pause = 0 by default.
synthetic_ccddata(camera="BLUE", binning=(2, 2), ampmode=None, shape=None,
                  **header_overrides)
    Factory: a KCCDData with data, uncertainty, mask, flags and a header that
    carries the keywords primitives read.
tmp_redux
    A temporary output directory with a ``redux/`` subdirectory and an empty
    Proctab.

Session-wide: astropy IERS auto-download is switched off (D18).
"""
import os
from types import SimpleNamespace

import numpy as np
import pytest
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from astropy.utils import iers
from keckdrpframework.config.framework_config import ConfigClass

from kcwidrp.core.kcwi_pkg_resources import get_resource_path
from kcwidrp.core.kcwi_proctab import Proctab
from kcwidrp.primitives.kcwi_file_primitives import KCCDData

# D18: no network access for IERS tables, for every test in the session.
iers.conf.auto_download = False

# ---------------------------------------------------------------------------
# Header keyword scan (T.1, 2026-10-22, commit b685577 + T.1 working tree).
#
# Produced by an AST scan of kcwidrp/primitives/*.py and kcwidrp/core/*.py for
# string-literal reads:  header['K'], hdr['K'], <x>.header['K'],
# <x>.header.get('K'), and <x>.get_keyword('K') (the last is how
# ingest_file reads most keywords; it was added to the scan because the plain
# header[...] scan misses IFUNUM, NUMOPEN, BNASNAM, LMPnSTAT, ...).
# Format: KEYWORD count-of-reads.  Keywords built with % or format()
# (ASECn, BSECn, CSECn, DSECn, TSECn, ATSECn, GAINn, BIASRNn, OSCNRNn,
# OSCNVALn, STACKFn, and BARTANG/RARTANG, selected by camera) are not
# literals; they are handled explicitly in synthetic_ccddata.
# ---------------------------------------------------------------------------
HEADER_KEYWORD_SCAN = """
AIRMASS 2
AMPMNUM 2
AMPMODE 6
ARCFL 4
ARCNO 3
AUTMSKTY 1
AUTOMASK 1
AVEXCOR 1
AVWVSIG 3
BAR0 2
BARSEP 2
BCWAVE 1
BFILTNAM 2
BGRANGLE 1
BGRATNAM 1
BINNING 6
BNASNAM 1
BPFILE 1
BUNIT 3
CAMERA 15
CBARSFL 3
CBARSNO 3
CCDCFG 9
CCDMODE 4
CCDSUM 1
CD1_1 3
CD1_2 3
CD2_1 3
CD2_2 3
CD3_3 9
CDELT1 4
CDELT2 2
CDSSPEED 1
CNAME1 6
CNAME2 5
CNAME3 3
COMMENT 2
CONTPOS 1
CONTWID 1
CRPIX1 8
CRPIX2 7
CRPIX3 5
CRVAL1 7
CRVAL2 5
CRVAL3 9
CTYPE1 7
CTYPE2 5
CTYPE3 6
CUNIT1 6
CUNIT2 5
CUNIT3 3
DARANG 1
DARCOR 4
DAREFWL 1
DARKSCL 1
DARPADX 1
DARPADY 2
DATE-BEG 1
DATE-END 1
DEC 1
DECBASE 1
DICHFRAC 4
EQUINOX 2
EXPTIME 1
EXTCOR 1
FLIMAGIN 1
FLSPECTR 1
FRAMENO 25
GAINMUL 4
GEOMCOR 4
GEOMFL 2
GROUPID 5
HISTORY 47
IFUNAM 3
IFUNUM 2
IFUPA 2
IFUROFF 1
IMTYPE 39
INVFORD 1
INVFW0 1
INVFW1 1
INVSENS 1
INVSL0 1
INVSL1 1
INVSLMX 1
INVSLY 1
INVSW0 1
INVSW1 1
INVSY0 1
INVSY1 1
INVSZ0 1
INVSZ1 1
KOAID 2
LATPOLE 3
LMP0NAM 1
LMP0SHST 1
LMP0STAT 1
LMP1NAM 1
LMP1SHST 1
LMP1STAT 1
LMP3STAT 2
LONPOLE 3
MASTFLAT 1
MBFILE 1
MDFILE 1
MFFILE 3
MJD 6
MSFILE 1
MSIMNO 1
NASSCL 2
NAXIS 1
NAXIS1 2
NAXIS2 4
NAXIS3 1
NBPCLEAN 1
NCRCLEAN 2
NSATFLAG 1
NSHFDN 1
NSHFUP 2
NSTACK 7
NUMOPEN 1
NVIDINP 7
OBJECT 6
OBJROW0 1
OBJROW1 1
OFNAME 4
PARANG 1
PONAME 1
POSMAPF 2
PXSCL 4
RA 1
RABASE 1
RADESYS 3
RCWAVE 2
RGRANGLE 1
RGRATNAM 1
RNASNAM 1
ROTPOSN 1
ROTREFAN 1
SDWVSIG 3
SHUFROWS 2
SKYCOR 4
SKYIMAGE 1
SKYMAST 2
SKYMODEL 1
SKYMSK 2
SKYMSKF 1
SKYROW0 1
SKYROW1 1
SKYSCL 1
SLIMAPF 2
SLSCL 4
STATEID 6
STATENAM 1
STCKMETH 6
STCKSIGU 6
SUBROW0 1
SUBROW1 1
TAPLINES 1
TARGDEC 2
TARGNAME 7
TARGRA 2
TELAPSE 1
TELESCOP 1
TTIME 9
TUCKALPH 1
VCORR 4
VCORRTYP 3
WAVALL0 8
WAVALL1 8
WAVGOOD0 9
WAVGOOD1 9
WAVMAPF 2
WAVMID 7
WCSDIM 4
WCSNAME 2
XPOSURE 5
history 5
"""

# Reads that are not instrument keywords (FITS structure, free text).
STRUCTURAL_KEYWORDS = {"COMMENT", "HISTORY", "history", "NAXIS", "NAXIS1",
                       "NAXIS2", "NAXIS3"}


def scanned_keywords():
    """Return the set of keywords in HEADER_KEYWORD_SCAN."""
    return {line.split()[0] for line in HEADER_KEYWORD_SCAN.strip().splitlines()}


# Plausible values for a KCWI BLUE, 2x2 binned, four-amplifier OBJECT frame.
# Values marked "?" in the trailing comment are guesses (not read from real
# data); see _UNSURE below.  None of them is a calibrated number: they exist so
# that header reads do not raise KeyError.
_DEFAULTS = {
    "AIRMASS": 1.1,
    "AMPMNUM": 0,              # fix_header(): AMPMODE not in its table -> 0
    "ARCFL": "kb260101_00010.fits",          # ?
    "ARCNO": 10,                             # ?
    "AUTMSKTY": "none",                      # ?
    "AUTOMASK": False,                       # ?
    "AVEXCOR": 1.0,                          # ?
    "AVWVSIG": 0.5,                          # ? Angstrom
    "BAR0": 0,                               # ?
    "BARSEP": 16.0,                          # ? pixels
    "BCWAVE": 4500.0,
    "BFILTNAM": "KBlue",
    "BGRANGLE": 61.0,                        # ? degrees
    "BGRATNAM": "BM",
    "BNASNAM": "Object",                     # 'Mask' in value -> nasmask()
    "BPFILE": "",                            # ?
    "BUNIT": "adu",
    "CBARSFL": "kb260101_00011.fits",        # ?
    "CBARSNO": 11,                           # ?
    "CCDMODE": 0,
    "CCDSUM": "2 2",
    "CD1_1": -0.0003,                        # ? deg (slicer scale ~0.37 arcsec)
    "CD1_2": 0.0,
    "CD2_1": 0.0,
    "CD2_2": 4.048e-05 * 2,                  # PIXSCALE * ybin (kcwi.cfg)
    "CD3_3": 1.0,                            # ? Angstrom/px, BM at 2x binning
    "CDELT1": 4.048e-05 * 2,
    "CDELT2": 3.7718e-04,                    # SLICESCALE (kcwi.cfg)
    "CDSSPEED": 0,                           # ?
    "CNAME1": "KCWI y",                      # ?
    "CNAME2": "KCWI x",                      # ?
    "CNAME3": "KCWI wave",                   # ?
    "CONTPOS": 0.0,                          # ?
    "CONTWID": 0.0,                          # ?
    "CRPIX1": 1.0,
    "CRPIX2": 1.0,
    "CRPIX3": 1.0,
    "CRVAL1": 0.0,
    "CRVAL2": 0.0,
    "CRVAL3": 3500.0,                        # ?
    "CTYPE1": "RA---TAN",
    "CTYPE2": "DEC--TAN",
    "CTYPE3": "AWAV",
    "CUNIT1": "deg",
    "CUNIT2": "deg",
    "CUNIT3": "Angstrom",
    "DARANG": 0.0,                           # ?
    "DARCOR": False,
    "DAREFWL": 4500.0,                       # ?
    "DARKSCL": 1.0,                          # ?
    "DARPADX": 0,                            # ?
    "DARPADY": 0,                            # ?
    "DATE-BEG": "2026-01-01T05:00:00.000",
    "DATE-END": "2026-01-01T05:00:01.000",
    "DEC": 2.0,
    "DECBASE": 2.0,                          # ?
    "DICHFRAC": 0.0,                         # ?
    "EQUINOX": 2000.0,
    "EXPTIME": 1.0,
    "EXTCOR": False,
    "FLIMAGIN": "off",
    "FLSPECTR": "off",
    "FRAMENO": 1,
    "GAINMUL": 1,                            # fix_header() sets 1
    "GEOMCOR": False,
    "GEOMFL": "",                            # ?
    "GROUPID": "kb260101_00001",             # ?
    "IFUNAM": "Large",
    "IFUNUM": 1,
    "IFUPA": 0.0,
    "IFUROFF": 0.0,                          # ?
    "IMTYPE": "OBJECT",
    "INVFORD": 3,                            # ?
    "INVFW0": 0.0,                           # ?
    "INVFW1": 0.0,                           # ?
    "INVSENS": 1.0,                          # ?
    "INVSL0": 0.0,                           # ?
    "INVSL1": 0.0,                           # ?
    "INVSLMX": 0.0,                          # ?
    "INVSLY": 0.0,                           # ?
    "INVSW0": 0.0,                           # ?
    "INVSW1": 0.0,                           # ?
    "INVSY0": 0,                             # ?
    "INVSY1": 0,                             # ?
    "INVSZ0": 0,                             # ?
    "INVSZ1": 0,                             # ?
    "KOAID": "KB.20260101.00001.00",         # ?
    "LATPOLE": 90.0,                         # ?
    "LMP0NAM": "FEAR",
    "LMP0SHST": 0,
    "LMP0STAT": 0,
    "LMP1NAM": "THAR",
    "LMP1SHST": 0,
    "LMP1STAT": 0,
    "LMP3STAT": 0,
    "LONPOLE": 180.0,                        # ?
    "MASTFLAT": "",                          # ?
    "MBFILE": "",                            # ?
    "MDFILE": "",                            # ?
    "MFFILE": "",                            # ?
    "MJD": 61041.2083,                       # = 2026-01-01T05:00 UT
    "MSFILE": "",                            # ?
    "MSIMNO": 0,                             # ?
    "NASSCL": False,                         # ?
    "NBPCLEAN": 0,
    "NCRCLEAN": 0,
    "NSATFLAG": 0,
    "NSHFDN": 0,
    "NSHFUP": 0,
    "NSTACK": 1,
    "NUMOPEN": 0,                            # ?
    "NVIDINP": 4,                            # overridden from ampmode
    "OBJECT": "Test",
    "OBJROW0": 0,                            # ?
    "OBJROW1": 0,                            # ?
    "OFNAME": "kb260101_00001.fits",
    "PARANG": 0.0,                           # ?
    "PONAME": "IFU",                         # ?
    "POSMAPF": "",                           # ?
    "PXSCL": 4.048e-05 * 2,                  # ? PIXSCALE * xbin
    "RA": 150.0,
    "RABASE": 150.0,                         # ?
    "RADESYS": "FK5",
    "RCWAVE": 7000.0,                        # ? presence means dichroic in
    "RGRANGLE": 0.0,                         # ?   BLUE headers (dich())
    "RGRATNAM": "RL",                        # ?
    "RNASNAM": "Object",
    "ROTPOSN": 0.0,                          # ?
    "ROTREFAN": 0.0,                         # ?
    "SDWVSIG": 0.1,                          # ?
    "SHUFROWS": 0,
    "SKYCOR": False,
    "SKYIMAGE": "",                          # ?
    "SKYMAST": "",                           # ?
    "SKYMODEL": "",                          # ?
    "SKYMSK": "",                            # ?
    "SKYMSKF": "",                           # ?
    "SKYROW0": 0,                            # ?
    "SKYROW1": 0,                            # ?
    "SKYSCL": 1.0,                           # ?
    "SLIMAPF": "",                           # ?
    "SLSCL": 3.7718e-04,                     # SLICESCALE (kcwi.cfg)
    "STATEID": "0123456789abcdef",           # ? hash-like string
    "STATENAM": "synthetic",                 # ?
    "STCKMETH": "average",                   # ?
    "STCKSIGU": 3.0,                         # ?
    "SUBROW0": 0,                            # ?
    "SUBROW1": 0,                            # ?
    "TAPLINES": 4,                           # overridden from ampmode
    "TARGDEC": 2.0,
    "TARGNAME": "Test",
    "TARGRA": 150.0,
    "TELAPSE": 1.0,
    "TELESCOP": "Keck II",
    "TTIME": 1.0,
    "TUCKALPH": 0.2,                         # TUKEYALPHA (kcwi.cfg)
    "VCORR": 0.0,
    "VCORRTYP": "none",                      # ?
    "WAVALL0": 3500.0,                       # ?
    "WAVALL1": 5500.0,                       # ?
    "WAVGOOD0": 3600.0,                      # ?
    "WAVGOOD1": 5400.0,                      # ?
    "WAVMAPF": "",                           # ?
    "WAVMID": 4500.0,
    "WCSDIM": 3,
    "WCSNAME": "KCWI",                       # ?
    "XPOSURE": 1.0,
}

# Values marked "?" in the comments above are guesses.

# Per-amplifier read noise and gain (e-/ADU, e-): plausible, not calibrated.
_GAINS = (1.00, 1.02, 1.04, 1.06)           # ?
_BIAS_RN = (3.0, 3.1, 3.2, 3.3)             # ?


def _sec(x0, x1, y0, y1):
    return "[%d:%d,%d:%d]" % (x0, x1, y0, y1)


def _amp_sections(nx, ny, n_amps, overscan, xbin=1, ybin=1):
    """Return dict of ASECn/BSECn/CSECn/DSECn (1-biased, IRAF style).

    Four amps: amps 1,2 lower left/right, 3,4 upper left/right (docstring of
    ingest_file.map_ccd).  Left amps read data then overscan; right amps read
    overscan then data (mirror image).  CSECn is in unbinned CCD pixels
    (map_ccd divides its span by the binning).  One amp: the whole frame, overscan on
    the right.  The exact on-chip layout of real frames is not reproduced.
    """
    out = {}
    if n_amps == 1:
        nd = nx - overscan
        out["ASEC1"] = _sec(1, nx, 1, ny)
        out["DSEC1"] = _sec(1, nd, 1, ny)
        out["BSEC1"] = _sec(nd + 1, nx, 1, ny)
        out["CSEC1"] = _sec(1, nd * xbin, 1, ny * ybin)
        return out
    bw, bh = nx // 2, ny // 2
    nd = bw - overscan
    for ia in range(1, 5):
        right = ia in (2, 4)
        upper = ia in (3, 4)
        y0, y1 = (bh + 1, ny) if upper else (1, bh)
        if right:
            bsec = _sec(bw + 1, bw + overscan, y0, y1)
            dsec = _sec(bw + overscan + 1, nx, y0, y1)
            cx = (2 * nd * xbin, nd * xbin + 1)   # descending: read backwards
        else:
            dsec = _sec(1, nd, y0, y1)
            bsec = _sec(nd + 1, bw, y0, y1)
            cx = (1, nd * xbin)
        cy = ((2 * bh * ybin, bh * ybin + 1) if upper else (1, bh * ybin))
        out["ASEC%d" % ia] = _sec(bw + 1 if right else 1, nx if right else bw,
                                  y0, y1)
        out["BSEC%d" % ia] = bsec
        out["DSEC%d" % ia] = dsec
        out["CSEC%d" % ia] = _sec(cx[0], cx[1], cy[0], cy[1])
    return out


@pytest.fixture(scope="session", autouse=True)
def _no_iers_download():
    """Re-assert D18 for the session (also set at import time above)."""
    iers.conf.auto_download = False
    yield


@pytest.fixture
def kcwi_config():
    """Factory for the parsed kcwi.cfg (see module docstring).

    The file is parsed with keckdrpframework's ConfigClass exactly as
    scripts/reduce_kcwi.py does:
    ``ConfigClass(get_resource_path('kcwidrp', 'configs/kcwi.cfg'),
    default_section='KCWI')``, and ``cwd`` is set as reduce_kcwi does.
    Overrides are written into the parsed properties, which is where
    ConfigClass.__getattr__ looks values up (keys are lower-case).
    """
    def make(**overrides):
        path = get_resource_path("kcwidrp", "configs/kcwi.cfg")
        config = ConfigClass(path, default_section="KCWI")
        config.cwd = os.getcwd()
        config.properties["plot_level"] = 0
        config.properties["plot_pause"] = 0
        for key, value in overrides.items():
            config.properties[key.lower()] = value
        return config
    return make


@pytest.fixture
def synthetic_ccddata():
    """Factory for a KCCDData frame with the header keywords primitives read."""
    def make(camera="BLUE", binning=(2, 2), ampmode=None, shape=None,
             overscan=8, seed=0, **header_overrides):
        camera = camera.upper()
        if ampmode is None:
            ampmode = "ALL" if camera == "BLUE" else "L2U2L1U1"
        n_amps = 4 if ampmode in ("ALL", "L2U2L1U1") else 1
        if shape is None:
            shape = (128, 144) if n_amps == 4 else (128, 72)
        ny, nx = shape
        hdr = fits.Header()
        values = dict(_DEFAULTS)
        values["CAMERA"] = camera
        values["BINNING"] = "%d,%d" % tuple(binning)
        values["CCDSUM"] = "%d %d" % tuple(binning)
        values["AMPMODE"] = ampmode
        values["NVIDINP"] = n_amps
        values["TAPLINES"] = n_amps
        values["PXSCL"] = 4.048e-05 * binning[0]
        values["CD2_2"] = 4.048e-05 * binning[1]
        values["CDELT1"] = 4.048e-05 * binning[1]
        values["BARTANG" if camera == "BLUE" else "RARTANG"] = 0.0   # ?
        values.update(_amp_sections(nx, ny, n_amps, overscan,
                                      binning[0], binning[1]))
        for ia in range(1, n_amps + 1):
            values["GAIN%d" % ia] = _GAINS[ia - 1]
            values["BIASRN%d" % ia] = _BIAS_RN[ia - 1]
        # CCDCFG as kcwi_fits_reader builds it (CCDSUM, CCDMODE, GAINMUL,
        # AMPMNUM).
        values["CCDCFG"] = ("%s%1d%02d%02d" % (
            values["CCDSUM"].replace(" ", ""), values["CCDMODE"],
            values["GAINMUL"], values["AMPMNUM"]))
        if camera == "RED":
            values.update({"RCWAVE": 7000.0, "RGRATNAM": "RL",
                           "RGRANGLE": 0.0, "RNASNAM": "Object"})
        values.update(header_overrides)
        for key, value in values.items():
            hdr[key] = value
        hdr["NAXIS1"] = nx
        hdr["NAXIS2"] = ny
        rng = np.random.default_rng(seed)
        data = 1000.0 + rng.normal(0.0, 3.0, size=shape)
        ccd = KCCDData(data, meta=hdr, unit="adu",
                       uncertainty=StdDevUncertainty(np.sqrt(np.abs(data))))
        # mask and flags are set as attributes, as CreateUncertaintyImage
        # does with uint8 zeros. (astropy 6.1.7's CCDData stores the mask as
        # bool whatever dtype is assigned; flags is a plain attribute.)
        ccd.mask = np.zeros(shape, dtype=np.uint8)
        ccd.flags = np.zeros(shape, dtype=np.uint8)
        return ccd
    return make


@pytest.fixture
def tmp_redux(tmp_path):
    """Temporary output directory with redux/ and an empty Proctab."""
    redux = tmp_path / "redux"
    redux.mkdir()
    proctab = Proctab()
    proctab.new_proctab()
    return SimpleNamespace(root=tmp_path, redux=redux, proctab=proctab,
                           proctab_file=redux / "kcwi.proc")


@pytest.fixture
def header_keyword_info():
    """The keyword scan, the structural exclusions and the default table."""
    return SimpleNamespace(scanned=scanned_keywords(),
                           structural=STRUCTURAL_KEYWORDS,
                           defaults=dict(_DEFAULTS))
