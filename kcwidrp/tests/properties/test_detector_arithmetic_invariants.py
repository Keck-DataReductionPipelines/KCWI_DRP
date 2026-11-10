"""Properties of the detector arithmetic primitives (task T.5).

Run through the T.2 harness (``run_primitive`` fixture) on frames from the
``synthetic_ccddata`` fixture: the smallest callable units, since
SubtractOverscan, CorrectGain and CreateUncertaintyImage have no helper
functions; each is one ``_perform``.

1. Overscan: the overscan strip measures the electronic pedestal of the
   frame, so adding the same constant pedestal to every pixel (data and
   overscan) must leave the overscan-subtracted image unchanged.  Raw KCWI
   pixels are integer DN, so the pedestal is an integer number of DN (the
   primitive rounds the overscan level to an integer DN, so a fractional
   pedestal is not a physical input).  The overscan-subtracted image is
   compared after TrimOverscan, i.e. the image the next primitive sees
   (untrimmed, the overscan columns themselves keep the pedestal).
2. Gain: converting DN to electrons multiplies each amplifier by its gain,
   so the correction is linear: G(a X + b Y) = a G(X) + b G(Y) for any
   frames X, Y and numbers a, b.
3. Uncertainty: a standard deviation is non-negative, and it is finite
   wherever the data are finite (Poisson term from the data, read noise
   from the header), whatever the sign of the data (bias-subtracted pixels
   can be negative).

The frames are small (4 amplifiers) and each example runs 2-6 primitives,
so these tests use a fraction of the profile's examples (``scaled``).
"""
import numpy as np
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra import numpy as hnp

from kcwidrp.primitives.CorrectGain import CorrectGain
from kcwidrp.primitives.CreateUncertaintyImage import CreateUncertaintyImage
from kcwidrp.primitives.SubtractOverscan import SubtractOverscan
from kcwidrp.primitives.TrimOverscan import TrimOverscan
from kcwidrp.tests.properties._profiles import scaled

# Default fixture frame: BLUE 2x2, 4 amplifiers with the real overscan widths
# (136 x 316 raw, 128 x 64 after TrimOverscan).
RAW = (136, 316)
TRIMMED = (128, 64)


@given(seed=st.integers(0, 2**32 - 1),
       level=st.integers(0, 3000), noise=st.integers(1, 20),
       pedestal=st.integers(-1000, 20000))
@scaled(0.6)
def test_constant_pedestal_does_not_change_overscan_subtracted_image(
        run_primitive, synthetic_ccddata, seed, level, noise, pedestal):
    # source: analytic - (D + c) - (O + c) = D - O for every pixel.
    rng = np.random.default_rng(seed)
    data = np.round(level + rng.normal(0.0, noise, RAW))
    out = []
    for c in (0, pedestal):
        frame = synthetic_ccddata()
        assert frame.shape == RAW
        frame.data = data + c          # integer-valued float64, exact
        sub = run_primitive(SubtractOverscan, frame)
        out.append(run_primitive(TrimOverscan, sub).ccddata.data)
    assert np.array_equal(out[0], out[1])


def _trimmed_args(run_primitive, synthetic_ccddata, data):
    """TrimOverscan output (ATSECn and map_ccd set) carrying ``data``."""
    args = run_primitive(TrimOverscan, synthetic_ccddata())
    args.ccddata.data = data.copy()
    return args


finite = st.floats(min_value=-1e5, max_value=1e5, allow_nan=False,
                   allow_infinity=False)
images = hnp.arrays(np.float64, TRIMMED, elements=finite)
coeff = st.floats(min_value=-8.0, max_value=8.0, allow_nan=False)


@given(x=images, y=images, a=coeff, b=coeff)
@scaled(0.6)
def test_gain_correction_is_linear(run_primitive, synthetic_ccddata,
                                   x, y, a, b):
    # source: analytic - G is a per-pixel multiplication by GAINn, so
    # G(aX + bY) = a G(X) + b G(Y).
    combo = a * x + b * y
    g_combo = run_primitive(CorrectGain, _trimmed_args(
        run_primitive, synthetic_ccddata, combo)).ccddata.data
    g_x = run_primitive(CorrectGain, _trimmed_args(
        run_primitive, synthetic_ccddata, x)).ccddata.data
    g_y = run_primitive(CorrectGain, _trimmed_args(
        run_primitive, synthetic_ccddata, y)).ccddata.data
    expected = a * g_x + b * g_y
    # Tolerance: both sides are the same real number rounded differently in
    # float64: G(aX + bY) has 4 roundings (a*x, b*y, the sum, *gain),
    # a G(X) + b G(Y) has 5, each <= eps/2 of a term bounded by
    # (|aX| + |bY|) * gain.  So |difference| <= 4.5 eps (|aX| + |bY|) gain;
    # the bound uses 8 eps and gain <= 0.2 e-/DN (fixture: GAINMUL 10).  A
    # gain applied twice, to the wrong amplifier, or an added offset breaks
    # it by many orders of magnitude.  Below 2.2e-308 (subnormal results)
    # rounding is absolute, <= half the smallest subnormal (4.9e-324) per
    # operation, hence the absolute term for 9 roundings (found by the
    # thorough profile: y ~ 1e-308, b = 4).
    scale = (np.abs(a * x) + np.abs(b * y)) * 0.2
    tiny = np.finfo(np.float64).smallest_subnormal
    assert np.all(np.abs(g_combo - expected)
                  <= 8 * 2.23e-16 * scale + 9 * tiny)


@given(data=hnp.arrays(np.float64, TRIMMED,
                       elements=st.one_of(
                           st.floats(min_value=-1e6, max_value=1e7),
                           st.just(np.nan))))
@scaled(0.6)
def test_uncertainty_is_nonnegative_and_finite_where_data_are(
        run_primitive, synthetic_ccddata, data):
    # source: analytic - sqrt(|D| + RN^2) is real, >= 0 and finite for
    # finite D and finite read noise.
    args = _trimmed_args(run_primitive, synthetic_ccddata, data)
    out = run_primitive(CreateUncertaintyImage, args).ccddata
    unc = out.uncertainty.array
    good = np.isfinite(data)
    assert unc.shape == data.shape
    assert np.all(np.isfinite(unc[good]))
    assert np.all(unc[good] >= 0)
