"""Characterization tests for kcwidrp.core.kcwi_get_std (task T.3).

Expected values: names of standards present in ``kcwidrp/data/stds`` (the
directory listing is the specification; the file names are lower-case with
whitespace removed, e.g. ``hz44.fits``, ``bd28d4211.fits``), and names that
are not in that directory.
"""
import logging
import os

import pytest
from astropy.io import fits

from kcwidrp.core.kcwi_get_std import kcwi_get_std, is_file_kcwi_std

STDS_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "data",
                        "stds")
LOG = logging.getLogger("test_get_std")


@pytest.mark.parametrize("target, name", [
    ("hz44", "hz44"),               # source: independent - data/stds/hz44.fits
    ("bd28d4211", "bd28d4211"),     # source: independent - bd28d4211.fits
    ("g191b2b", "g191b2b"),         # source: independent - g191b2b.fits
    ("feige110", "feige110"),       # source: independent - feige110.fits
])
def test_known_standard_found(target, name):
    stdfile, stdname = kcwi_get_std(target, logger=LOG)
    assert stdname == name
    assert os.path.basename(stdfile) == name + ".fits"
    assert os.path.exists(stdfile)
    # the file returned is the one in the data/stds directory
    assert os.path.samefile(stdfile, os.path.join(STDS_DIR, name + ".fits"))


@pytest.mark.parametrize("target, name", [
    ("HZ44", "hz44"),               # upper case
    ("Hz44", "hz44"),               # mixed case
    ("HZ 44", "hz44"),              # embedded space
    ("  hz44 ", "hz44"),            # leading and trailing space
    ("Feige 110", "feige110"),      # source: independent - feige110.fits
    ("BD28D4211", "bd28d4211"),
])
def test_case_and_whitespace_variants(target, name):
    # source: independent - data/stds file names are lower case, no spaces
    stdfile, stdname = kcwi_get_std(target, logger=LOG)
    assert stdname == name
    assert os.path.basename(stdfile) == name + ".fits"


@pytest.mark.parametrize("target", ["NotAStar", "", "hz4444", "vega"])
def test_non_standard_returns_none(target):
    # source: independent - none of these names has a file in data/stds
    assert kcwi_get_std(target, logger=LOG) == (None, None)


@pytest.mark.xfail(strict=True, reason="B18: kcwi_get_std(targname) "
                   "without a logger raises AttributeError although "
                   "logger=None is the default (kcwi_get_std.py:14,17)")
def test_get_std_without_logger():
    # source: analytic - the signature offers logger=None, so a call without
    # a logger must work and find hz44.
    stdfile, stdname = kcwi_get_std("hz44")
    assert stdname == "hz44"


def _write(tmp_path, name, imtype, targname):
    path = tmp_path / name
    hdr = fits.Header({"IMTYPE": imtype, "TARGNAME": targname})
    fits.PrimaryHDU(data=None, header=hdr).writeto(path)
    return str(path)


def test_is_std_true_for_object_frame_on_standard(tmp_path):
    # source: independent - OBJECT frame whose TARGNAME is in data/stds
    f = _write(tmp_path, "a.fits", "OBJECT", "HZ 44")
    assert is_file_kcwi_std(f, logger=LOG) is True


def test_is_std_false_for_object_frame_on_science_target(tmp_path):
    # source: independent - 'NGC1234' has no file in data/stds
    f = _write(tmp_path, "b.fits", "OBJECT", "NGC1234")
    assert is_file_kcwi_std(f, logger=LOG) is False


def test_is_std_false_for_non_object_frame(tmp_path):
    # source: analytic - only OBJECT frames are standards; an ARCLAMP frame
    # with a standard's name in TARGNAME is not
    f = _write(tmp_path, "c.fits", "ARCLAMP", "hz44")
    assert is_file_kcwi_std(f, logger=LOG) is False


def test_is_std_false_for_missing_file(tmp_path):
    # source: analytic - a missing file is reported as False, not raised
    assert is_file_kcwi_std(str(tmp_path / "missing.fits"), logger=LOG) is False


@pytest.mark.xfail(strict=True, reason="B18: is_file_kcwi_std on a "
                   "missing file without a logger raises AttributeError "
                   "(kcwi_get_std.py:32-33 calls logger.error on None)")
def test_is_std_missing_file_without_logger(tmp_path):
    # source: analytic - logger=None is the default, so this must return False
    assert is_file_kcwi_std(str(tmp_path / "missing.fits")) is False
