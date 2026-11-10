"""Property: parse_imsec inverts IRAF image-section notation (task T.5).

IRAF image sections (IRAF "imsec" convention, as in the ASECn/BSECn/DSECn
keywords of KCWI headers) are written ``[x1:x2,y1:y2]`` or with a step,
``[x1:x2:s,y1:y2:t]``: 1-based, both limits inclusive, x (column) first.  A
first limit larger than the second means the axis is read in descending
order; the step is a positive magnitude.

parse_imsec documents its result as zero-based python slice limits
``(y0, y1, x0, x1)`` (inclusive upper limit) and strides ``(ystride,
xstride)``, a negative stride marking a descending axis.  ``_format`` below
is the inverse written from that convention, not from the code: format then
parse gives back the limits and strides, and parse then format gives back the
canonical section string (no explicit step when the step is 1).

parse_imsec is not exercised by T.3 on random input; T.3 pins 7 hand-picked
sections.
"""
from hypothesis import given
from hypothesis import strategies as st

from kcwidrp.primitives.kcwi_file_primitives import parse_imsec

# KCWI detectors are 4096 x 4096 plus overscan; any 1-based limit up to 5000
# is a valid section limit.
limit = st.integers(min_value=1, max_value=5000)
step = st.integers(min_value=1, max_value=16)


def _axis(lo0, hi0, stride):
    """IRAF text for one axis from zero-based limits and signed stride."""
    a, b = lo0 + 1, hi0 + 1                  # 1-based, inclusive
    if stride < 0:
        a, b = b, a                          # descending: first > second
    text = "%d:%d" % (a, b)
    if abs(stride) != 1:
        text += ":%d" % abs(stride)
    return text


def _format(sec, stride):
    y0, y1, x0, x1 = sec
    ystride, xstride = stride
    return "[%s,%s]" % (_axis(x0, x1, xstride), _axis(y0, y1, ystride))


@st.composite
def axis_spec(draw):
    """(zero-based lo, hi, signed stride) for one axis."""
    lo, hi = sorted((draw(limit), draw(limit)))
    s = draw(step)
    # A descending axis needs first > second; with equal limits the axis is
    # ascending by definition (the convention cannot express direction).
    descending = lo < hi and draw(st.booleans())
    return lo - 1, hi - 1, -s if descending else s


@given(x=axis_spec(), y=axis_spec())
def test_format_then_parse_returns_limits_and_strides(x, y):
    # source: analytic - IRAF section convention (module docstring)
    sec = (y[0], y[1], x[0], x[1])
    stride = (y[2], x[2])
    got_sec, got_stride = parse_imsec(_format(sec, stride))
    assert tuple(got_sec) == sec
    assert tuple(got_stride) == stride


@given(x1=limit, x2=limit, y1=limit, y2=limit,
       xs=st.one_of(st.none(), step), ys=st.one_of(st.none(), step))
def test_parse_then_format_returns_canonical_section(x1, x2, y1, y2, xs, ys):
    # source: analytic - IRAF section convention; the canonical form omits a
    # step of 1 (``[1:10:1,...]`` and ``[1:10,...]`` are the same section).
    def ax(a, b, s):
        return "%d:%d" % (a, b) + ("" if s is None else ":%d" % s)

    text = "[%s,%s]" % (ax(x1, x2, xs), ax(y1, y2, ys))
    canonical = "[%s,%s]" % (ax(x1, x2, None if xs in (None, 1) else xs),
                             ax(y1, y2, None if ys in (None, 1) else ys))
    assert _format(*parse_imsec(text)) == canonical
