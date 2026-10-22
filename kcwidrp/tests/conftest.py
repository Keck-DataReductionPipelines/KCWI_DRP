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


# Header values for a KCWI BLUE, 2x2 binned, four-amplifier OBJECT frame.
# Trailing comments say where each value comes from:
#   "# kb260115_00011": copied from a real blue arc (2026-01-15, BL/RM1,
#       Small slicer), after kcwi_fits_reader()/fix_header() where ingest
#       changes the value.  The frame is 1x1 binned; binning-dependent values
#       are computed in synthetic_ccddata() instead.
#   "# object value: synthetic": object-frame keywords (pointing, target,
#       timing); 2026jan15 has no object frame, so these are placeholders.
#   "# DRP-written": keywords that primitives write; they never appear in
#       raw headers, and the value only has to be of the right type.
#   "# ?": still an unverified guess.
_DEFAULTS = {
    "AIRMASS": 1.1,                          # object value: synthetic (no object frame yet)
    "AMPMNUM": 0,              # fix_header(): AMPMODE not in its table -> 0
    "ARCFL": "kb260101_00010.fits",          # DRP-written: not in raw headers
    "ARCNO": 10,                             # DRP-written: not in raw headers
    "AUTMSKTY": "none",                      # DRP-written: not in raw headers
    "AUTOMASK": False,                       # DRP-written: not in raw headers
    "AVEXCOR": 1.0,                          # DRP-written: not in raw headers
    "AVWVSIG": 0.5,                          # DRP-written: not in raw headers (Angstrom)
    "BAR0": 0,                               # DRP-written: not in raw headers
    "BARSEP": 16.0,                          # DRP-written: not in raw headers (pixels)
    "BCWAVE": 4599.8271484,                  # kb260115_00011
    "BFILTNAM": 'None',                      # kb260115_00011
    "BGRANGLE": 28.7556705,                  # kb260115_00011
    "BGRATNAM": 'BL',                        # kb260115_00011
    "BNASNAM": 'Open',                       # kb260115_00011
    "BPFILE": "",                            # DRP-written: not in raw headers
    "BUNIT": "adu",
    "CBARSFL": "kb260101_00011.fits",        # DRP-written: not in raw headers
    "CBARSNO": 11,                           # DRP-written: not in raw headers
    "CCDMODE": 1,                            # kb260115_00011
    "CCDSUM": "2 2",
    "CD1_1": -0.0003,                        # DRP-written: not in raw headers (deg (slicer scale ~0.37 arcsec))
    "CD1_2": 0.0,
    "CD2_1": 0.0,
    "CD2_2": 4.048e-05 * 2,                  # PIXSCALE * ybin (kcwi.cfg)
    "CD3_3": 1.0,                            # DRP-written: not in raw headers (Angstrom/px, BM at 2x binning)
    "CDELT1": 4.048e-05 * 2,
    "CDELT2": 3.7718e-04,                    # SLICESCALE (kcwi.cfg)
    "CDSSPEED": 0,                           # raw RED only (kr260115_00013: 1)
    "CNAME1": "KCWI y",                      # DRP-written: not in raw headers
    "CNAME2": "KCWI x",                      # DRP-written: not in raw headers
    "CNAME3": "KCWI wave",                   # DRP-written: not in raw headers
    "CONTPOS": 0.0,                          # DRP-written: not in raw headers
    "CONTWID": 0.0,                          # DRP-written: not in raw headers
    "CRPIX1": 1.0,
    "CRPIX2": 1.0,
    "CRPIX3": 1.0,
    "CRVAL1": 0.0,
    "CRVAL2": 0.0,
    "CRVAL3": 3500.0,                        # DRP-written: not in raw headers
    "CTYPE1": "RA---TAN",
    "CTYPE2": "DEC--TAN",
    "CTYPE3": "AWAV",
    "CUNIT1": "deg",
    "CUNIT2": "deg",
    "CUNIT3": "Angstrom",
    "DARANG": 0.0,                           # DRP-written: not in raw headers
    "DARCOR": False,
    "DAREFWL": 4500.0,                       # DRP-written: not in raw headers
    "DARKSCL": 1.0,                          # DRP-written: not in raw headers
    "DARPADX": 0,                            # DRP-written: not in raw headers
    "DARPADY": 0,                            # DRP-written: not in raw headers
    "DATE-BEG": "2026-01-01T05:00:00.000",   # object value: synthetic (no object frame yet)
    "DATE-END": "2026-01-01T05:00:01.000",   # object value: synthetic (no object frame yet)
    "DEC": 2.0,                              # object value: synthetic (no object frame yet)
    "DECBASE": 2.0,                          # object value: synthetic (no object frame yet)
    "DICHFRAC": 0.0,                         # DRP-written: not in raw headers
    "EQUINOX": 2000.0,                       # kb260115_00011
    "EXPTIME": 1.0,                          # object value: synthetic (no object frame yet)
    "EXTCOR": False,
    "FLIMAGIN": "off",
    "FLSPECTR": "off",
    "FRAMENO": 1,                            # object value: synthetic (no object frame yet)
    "GAINMUL": 10,                           # kb260115_00011 (post-ingest)
    "GEOMCOR": False,
    "GEOMFL": "",                            # DRP-written: not in raw headers
    "GROUPID": "kb260101_00001",             # object value: synthetic (no object frame yet)
    "IFUNAM": 'Small',                       # kb260115_00011
    "IFUNUM": 3,                             # kb260115_00011
    "IFUPA": 0.0,                            # object value: synthetic (no object frame yet)
    "IFUROFF": 0.0,                          # DRP-written: not in raw headers
    "IMTYPE": "OBJECT",                      # object value: synthetic (no object frame yet)
    "INVFORD": 3,                            # DRP-written: not in raw headers
    "INVFW0": 0.0,                           # DRP-written: not in raw headers
    "INVFW1": 0.0,                           # DRP-written: not in raw headers
    "INVSENS": 1.0,                          # DRP-written: not in raw headers
    "INVSL0": 0.0,                           # DRP-written: not in raw headers
    "INVSL1": 0.0,                           # DRP-written: not in raw headers
    "INVSLMX": 0.0,                          # DRP-written: not in raw headers
    "INVSLY": 0.0,                           # DRP-written: not in raw headers
    "INVSW0": 0.0,                           # DRP-written: not in raw headers
    "INVSW1": 0.0,                           # DRP-written: not in raw headers
    "INVSY0": 0,                             # DRP-written: not in raw headers
    "INVSY1": 0,                             # DRP-written: not in raw headers
    "INVSZ0": 0,                             # DRP-written: not in raw headers
    "INVSZ1": 0,                             # DRP-written: not in raw headers
    "KOAID": "KB.20260101.00001.00",         # object value: synthetic (no object frame yet)
    "LATPOLE": 90.0,                         # DRP-written: not in raw headers
    "LMP0NAM": 'FeAr',                       # kb260115_00011
    "LMP0SHST": 0,
    "LMP0STAT": 0,
    "LMP1NAM": 'ThAr',                       # kb260115_00011
    "LMP1SHST": 0,
    "LMP1STAT": 0,
    "LMP3STAT": 0,
    "LONPOLE": 180.0,                        # DRP-written: not in raw headers
    "MASTFLAT": "",                          # DRP-written: not in raw headers
    "MBFILE": "",                            # DRP-written: not in raw headers
    "MDFILE": "",                            # DRP-written: not in raw headers
    "MFFILE": "",                            # DRP-written: not in raw headers
    "MJD": 61041.2083,                       # object value: synthetic (no object frame yet)
    "MSFILE": "",                            # DRP-written: not in raw headers
    "MSIMNO": 0,                             # DRP-written: not in raw headers
    "NASSCL": False,                         # DRP-written: not in raw headers
    "NBPCLEAN": 0,
    "NCRCLEAN": 0,
    "NSATFLAG": 0,
    "NSHFDN": 0,
    "NSHFUP": 0,
    "NSTACK": 1,
    "NUMOPEN": 1,                            # kb260115_00011
    "NVIDINP": 4,                            # overridden from ampmode
    "OBJECT": "Test",                        # object value: synthetic (no object frame yet)
    "OBJROW0": 0,                            # DRP-written: not in raw headers
    "OBJROW1": 0,                            # DRP-written: not in raw headers
    "OFNAME": "kb260101_00001.fits",         # object value: synthetic (no object frame yet)
    "PARANG": 0.0,                           # object value: synthetic (no object frame yet)
    "PONAME": 'REF-N',                       # kb260115_00011
    "POSMAPF": "",                           # DRP-written: not in raw headers
    "PXSCL": 4.048e-05 * 2,                  # DRP-written: not in raw headers (PIXSCALE * xbin)
    "RA": 150.0,                             # object value: synthetic (no object frame yet)
    "RABASE": 150.0,                         # object value: synthetic (no object frame yet)
    "RADESYS": "FK5",
    "RCWAVE": 7199.9848633,                  # kb260115_00011
    "RGRANGLE": 123.0774384,                 # kb260115_00011
    "RGRATNAM": 'RM1',                       # kb260115_00011
    "RNASNAM": 'Open',                       # kb260115_00011
    "ROTPOSN": 0.0,                          # object value: synthetic (no object frame yet)
    "ROTREFAN": 0.0,                         # object value: synthetic (no object frame yet)
    "SDWVSIG": 0.1,                          # DRP-written: not in raw headers
    "SHUFROWS": 1370,                        # kb260115_00011
    "SKYCOR": False,
    "SKYIMAGE": "",                          # DRP-written: not in raw headers
    "SKYMAST": "",                           # DRP-written: not in raw headers
    "SKYMODEL": "",                          # DRP-written: not in raw headers
    "SKYMSK": "",                            # DRP-written: not in raw headers
    "SKYMSKF": "",                           # DRP-written: not in raw headers
    "SKYROW0": 0,                            # DRP-written: not in raw headers
    "SKYROW1": 0,                            # DRP-written: not in raw headers
    "SKYSCL": 1.0,                           # DRP-written: not in raw headers
    "SLIMAPF": "",                           # DRP-written: not in raw headers
    "SLSCL": 3.7718e-04,                     # SLICESCALE (kcwi.cfg)
    "STATEID": '6968443881c6cc2a85772b8f',   # kb260115_00011
    "STATENAM": 'Small_BL_4600_RM1_7200',    # kb260115_00011
    "STCKMETH": "average",                   # DRP-written: not in raw headers
    "STCKSIGU": 3.0,                         # DRP-written: not in raw headers
    "SUBROW0": 0,                            # DRP-written: not in raw headers
    "SUBROW1": 0,                            # DRP-written: not in raw headers
    "TAPLINES": 4,                           # overridden from ampmode
    "TARGDEC": 2.0,                          # object value: synthetic (no object frame yet)
    "TARGNAME": "Test",                      # object value: synthetic (no object frame yet)
    "TARGRA": 150.0,                         # object value: synthetic (no object frame yet)
    "TELAPSE": 1.0,                          # object value: synthetic (no object frame yet)
    "TELESCOP": 'Keck II',                   # kb260115_00011
    "TTIME": 1.0,                            # object value: synthetic (no object frame yet)
    "TUCKALPH": 0.2,                         # TUKEYALPHA (kcwi.cfg)
    "VCORR": 0.0,
    "VCORRTYP": "none",                      # DRP-written: not in raw headers
    "WAVALL0": 3500.0,                       # DRP-written: not in raw headers
    "WAVALL1": 5500.0,                       # DRP-written: not in raw headers
    "WAVGOOD0": 3600.0,                      # DRP-written: not in raw headers
    "WAVGOOD1": 5400.0,                      # DRP-written: not in raw headers
    "WAVMAPF": "",                           # DRP-written: not in raw headers
    "WAVMID": 4500.0,
    "WCSDIM": 3,
    "WCSNAME": "KCWI",                       # DRP-written: not in raw headers
    "XPOSURE": 1.0,                          # object value: synthetic (no object frame yet)
}

# Values marked "?" in the comments above are guesses.

# Per-amplifier gain (e-/ADU) as kcwi_fits_reader()/fix_header() set it for
# kb260115_00011 (BLUE, GAINMUL 10, amps 1-4).  RED frames use amps numbered
# from 0 (GAIN0..GAIN3) with fix_header's red_amp_gain table; the fixture's RED
# variant does not model that yet (T.4).
_GAINS = (0.157, 0.160, 0.158, 0.158)       # kb260115_00011 (post-ingest)
# Read noise per amplifier (e-): DRP-written (not in raw headers), plausible.
_BIAS_RN = (3.0, 3.1, 3.2, 3.3)             # DRP-written


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
        values["BARTANG" if camera == "BLUE" else "RARTANG"] = (
            23.1307564 if camera == "BLUE" else 56.022934)  # kb/kr260115 arcs
        values.update(_amp_sections(nx, ny, n_amps, overscan,
                                      binning[0], binning[1]))
        for ia in range(1, n_amps + 1):
            values["GAIN%d" % ia] = _GAINS[ia - 1]
            values["BIASRN%d" % ia] = _BIAS_RN[ia - 1]
        if camera == "RED":
            # kr260115_00013 shares the instrument state of the blue arc.
            values.update({"RCWAVE": 7199.9848633, "RGRATNAM": 'RM1',
                           "RGRANGLE": 123.0774384, "RNASNAM": 'Open',
                           "GAINMUL": 1, "CCDMODE": 1})
        # CCDCFG as kcwi_fits_reader builds it (CCDSUM, CCDMODE, GAINMUL,
        # AMPMNUM).
        values["CCDCFG"] = ("%s%1d%02d%02d" % (
            values["CCDSUM"].replace(" ", ""), values["CCDMODE"],
            values["GAINMUL"], values["AMPMNUM"]))
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
