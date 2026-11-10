"""Characterization tests for kcwidrp.core.kcwi_correct_extin (task T.3).

Expected values come from the text table ``data/extin/snfext.dat`` (Buton et
al. 2013, A&A 549, A8; columns: wavelength in Angstrom, extinction in
mag/airmass), read here with numpy.loadtxt.  The code under test reads the
other file of the pair, ``snfext.fits``, so the two are independent copies of
the same table.  The test wavelengths lie on the 2 A nodes of the table, where
any interpolation must reproduce the tabulated value.

The flux ratio for extinction k (mag/airmass) at airmass X is the definition
of magnitudes: F_true / F_obs = 10**(0.4 * k * X).
"""
import logging
import os

import numpy as np
import pytest
from astropy.io import fits

from kcwidrp.core.kcwi_correct_extin import kcwi_correct_extin

DAT = os.path.join(os.path.dirname(__file__), "..", "..", "..", "data",
                   "extin", "snfext.dat")
LOG = logging.getLogger("test_extin")

# Tolerance: the .dat file has 7 significant digits and the FITS copy is
# float32 (relative 6e-8); the ratio's exponent is 0.4*k*X*ln10 <= ~1, so the
# ratio differs by < 1e-7 relative.  1e-6 leaves a factor 10.
RTOL = 1e-6

NW = 10
W0 = 4000.0   # a node of the table
DW = 2.0      # the table spacing


@pytest.fixture(scope="module")
def table():
    lam, ext = np.loadtxt(DAT, usecols=(0, 1), unpack=True)
    return lam, ext


def _expected_ratio(table, airmass):
    lam, ext = table
    idx = np.rint((W0 + DW * np.arange(NW) - lam[0]) / DW).astype(int)
    assert np.array_equal(lam[idx], W0 + DW * np.arange(NW))  # on nodes
    return 10.0 ** (0.4 * ext[idx] * airmass)


def _hdr(airmass):
    return fits.Header({"AIRMASS": airmass, "CD3_3": DW, "CRVAL3": W0})


@pytest.mark.parametrize("airmass", [1.0, 1.5, 2.0])
def test_vector_ratio_matches_table(table, airmass):
    # source: independent - snfext.dat (Buton et al. 2013) and the magnitude
    # definition 10**(0.4 k X)
    img = np.ones(NW)
    kcwi_correct_extin(img, _hdr(airmass), logger=LOG)
    np.testing.assert_allclose(img, _expected_ratio(table, airmass),
                               rtol=RTOL, atol=0)


def test_ratio_between_two_airmasses(table):
    # source: independent - the ratio of corrections at X=2.0 and X=1.0 is
    # 10**(0.4 k (2.0 - 1.0)) from snfext.dat
    a = np.ones(NW)
    b = np.ones(NW)
    kcwi_correct_extin(a, _hdr(1.0), logger=LOG)
    kcwi_correct_extin(b, _hdr(2.0), logger=LOG)
    lam, ext = table
    idx = np.rint((W0 + DW * np.arange(NW) - lam[0]) / DW).astype(int)
    np.testing.assert_allclose(b / a, 10.0 ** (0.4 * ext[idx] * 1.0),
                               rtol=RTOL, atol=0)


def test_airmass_one_correction_is_positive_boost(table):
    # source: analytic - extinction k > 0 mag/airmass in this range, so the
    # corrected flux is larger than the input at every wavelength
    img = np.ones(NW)
    kcwi_correct_extin(img, _hdr(1.2), logger=LOG)
    assert np.all(img > 1.0)


def test_cube_gets_same_ratio_in_every_spaxel(table):
    # source: independent - snfext.dat; the correction depends on wavelength
    # (axis 0) only, and scales the data multiplicatively
    cube = np.full((NW, 2, 3), 5.0)
    kcwi_correct_extin(cube, _hdr(1.3), logger=LOG)
    expect = 5.0 * _expected_ratio(table, 1.3)
    for iy in range(2):
        for ix in range(3):
            np.testing.assert_allclose(cube[:, iy, ix], expect, rtol=RTOL,
                                       atol=0)


def test_header_keywords_updated(table):
    # source: independent - AVEXCOR is the mean of the flux ratios from
    # snfext.dat; EXTCOR and HISTORY are set as the primitive's output contract
    hdr = _hdr(1.4)
    img = np.ones(NW)
    kcwi_correct_extin(img, hdr, logger=LOG)
    assert hdr["EXTCOR"] is True
    np.testing.assert_allclose(hdr["AVEXCOR"],
                               _expected_ratio(table, 1.4).mean(),
                               rtol=RTOL, atol=0)
    assert "kcwi_correct_extin" in [str(h) for h in hdr["HISTORY"]]


def test_works_without_logger(table, capsys):
    # source: analytic - a missing logger falls back to print()
    img = np.ones(NW)
    kcwi_correct_extin(img, _hdr(1.0))
    assert "Extinction corrected" in capsys.readouterr().out
    np.testing.assert_allclose(img, _expected_ratio(table, 1.0), rtol=RTOL,
                               atol=0)
