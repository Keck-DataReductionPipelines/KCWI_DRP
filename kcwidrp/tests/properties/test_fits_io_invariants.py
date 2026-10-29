"""Property: kcwi_fits_writer -> kcwi_fits_reader round-trips a frame (T.5).

What is written must be read back: the data, uncertainty, mask and flags
arrays (data containing NaN, mask with arbitrary masked pixels, flags with
any uint8 value) and every header card added to the frame.

The writer stores data and uncertainty as float32 (documented in its
docstring), so the arrays are drawn as float32 values (held in float64, as
in the DRP): the round trip is then exact.  NaN must stay NaN at the same
pixels (FITS IEEE floats carry NaN).

Header cards are random but valid FITS (FITS standard 4.0, section 4):
keyword names of up to 8 characters from A-Z, 0-9, '-' and '_'; values that
are logical, integer (64-bit), real (finite float64) or character strings
of printable ASCII (0x20-0x7E).  Keyword names start with 'HY' so that they
cannot collide with structural keywords (SIMPLE, NAXISn, ...), with the
keywords the reader rewrites (fix_header: GAINn, ...; CCDCFG; BUNIT) or
with the frame's own cards.  Trailing blanks in a string value are not
significant in FITS (standard 4.2.1.1), so the expected string has them
stripped, and real values are compared to the 13 significant digits a
20-character FITS value field holds.

The frame is built here, not from the synthetic_ccddata fixture, with the
few keywords the reader needs (CAMERA, GAINMUL, NVIDINP, AMPID1, CCDCFG) and
a HISTORY card (the writer requires one: bug B20, pinned in T.3).  T.3 checks
one fixed frame; here the arrays and cards are random.
"""
import numpy as np
import pytest
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra import numpy as hnp

from kcwidrp.primitives.kcwi_file_primitives import (KCCDData,
                                                     kcwi_fits_reader,
                                                     kcwi_fits_writer)
from kcwidrp.tests.properties._profiles import scaled

SHAPE = st.tuples(st.integers(1, 12), st.integers(1, 12))

f32 = st.floats(width=32, allow_nan=False, allow_infinity=False)


def _float32_arrays(shape, nan=True):
    elements = st.one_of(f32, st.just(np.nan)) if nan else \
        st.floats(min_value=0.0, width=32, allow_infinity=False)
    return hnp.arrays(np.float32, shape, elements=elements)


keyword = st.text(alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_",
                  min_size=0, max_size=6).map(lambda s: "HY" + s)
card_value = st.one_of(
    st.booleans(),
    st.integers(min_value=-2**63, max_value=2**63 - 1),
    st.floats(allow_nan=False, allow_infinity=False),
    st.text(alphabet=st.characters(min_codepoint=0x20, max_codepoint=0x7E),
            max_size=60),
)
cards = st.dictionaries(keyword, card_value, max_size=8)


@st.composite
def frames(draw):
    shape = draw(SHAPE)
    data = draw(_float32_arrays(shape))
    unc = draw(_float32_arrays(shape, nan=False))
    mask = draw(hnp.arrays(np.bool_, shape))
    flags = draw(hnp.arrays(np.uint8, shape))
    return data, unc, mask, flags, draw(cards)


def _header(extra):
    hdr = fits.Header()
    # Keywords the reader needs (fix_header for a BLUE frame; CCDCFG present
    # so the reader does not rebuild it).  Values as in the KCWI blue frame
    # kb260115_00011 (see the conftest synthetic_ccddata defaults).
    hdr["CAMERA"] = "BLUE"
    hdr["GAINMUL"] = 10
    hdr["NVIDINP"] = 1
    hdr["AMPID1"] = 0
    hdr["CCDCFG"] = "2210100"
    hdr.add_history("property test frame")
    for key, value in extra.items():
        hdr[key] = value
    return hdr


@given(frame=frames())
@scaled(0.6)
def test_writer_reader_round_trip(tmp_path_factory, frame):
    data, unc, mask, flags, extra = frame
    out_dir = tmp_path_factory.mktemp("io")
    ccd = KCCDData(data.astype(np.float64), meta=_header(extra), unit="adu",
                   uncertainty=StdDevUncertainty(unc.astype(np.float64)))
    ccd.mask = mask.copy()
    ccd.flags = flags.copy()
    kcwi_fits_writer(ccd, output_file="rt.fits", output_dir=str(out_dir))
    back, _ = kcwi_fits_reader(str(out_dir / "rt.fits"))

    # source: analytic - what was written is read back (module docstring)
    assert np.array_equal(back.data, data.astype(np.float64), equal_nan=True)
    assert np.array_equal(back.uncertainty.array, unc)
    assert np.array_equal(np.asarray(back.mask, dtype=bool), mask)
    assert back.flags.dtype == np.uint8
    assert np.array_equal(back.flags, flags)
    for key, value in extra.items():
        got = back.header[key]
        if isinstance(value, str):
            assert got == value.rstrip(" "), key
        elif isinstance(value, float):
            assert type(got) is float, key
            # Tolerance: a FITS fixed-format real value has 20 characters
            # (FITS 4.0, section 4.2.4: columns 11-30).  With a sign, the
            # point and an exponent like 'E-300' that leaves 13 significant
            # digits, and astropy truncates longer representations, so the
            # relative error is < 1e-12.  (astropy/FITS behaviour, not DRP.)
            assert got == pytest.approx(value, rel=1e-12, abs=0), key
        else:
            assert type(got) is type(value), key
            assert got == value, key
