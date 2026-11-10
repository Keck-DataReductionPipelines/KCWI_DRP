"""Characterization tests for GetAtlasLines.gaus, get_line_window and
findpeaks (task T.3).

Synthetic spectra: Gaussian emission lines of known amplitude, centre and
sigma on a zero background, sampled at integer pixels.
"""
import numpy as np
import pytest

from kcwidrp.primitives.GetAtlasLines import findpeaks, gaus, get_line_window

X = np.arange(400.0)
SIGMA = 2.0
# (centre, amplitude): centres are deliberately not on pixel centres, because
# a peak exactly on a pixel centre of a noiseless symmetric line is reported
# twice by findpeaks (see the report for T.3).
LINES = [(50.3, 1000.0), (120.4, 500.0), (200.6, 2000.0), (300.25, 800.0)]


def _spectrum(lines=LINES, sigma=SIGMA):
    y = np.zeros_like(X)
    for c, a in lines:
        y += a * np.exp(-(X - c) ** 2 / (2.0 * sigma ** 2))
    return y


# ---------------------------------------------------------------- gaus
def test_gaus_peak_value_is_amplitude():
    # source: analytic - at x = mu the exponent is 0
    assert gaus(3.0, 7.5, 3.0, 2.0) == 7.5


def test_gaus_one_sigma_and_two_sigma():
    # source: analytic - exp(-1/2) = 0.6065306597126334 at mu +/- sigma;
    # exp(-2) = 0.1353352832366127 at mu +/- 2 sigma
    assert gaus(5.0, 10.0, 3.0, 2.0) == pytest.approx(
        10.0 * 0.6065306597126334, rel=1e-15, abs=0)  # 1 ulp rounding of exp
    assert gaus(-1.0, 10.0, 3.0, 2.0) == pytest.approx(
        10.0 * 0.1353352832366127, rel=1e-15, abs=0)  # 1 ulp rounding of exp


def test_gaus_is_symmetric_about_mu():
    # source: analytic - G(mu + d) == G(mu - d)
    d = np.array([0.5, 1.0, 3.0, 7.0])
    assert np.array_equal(gaus(10.0 + d, 2.0, 10.0, 1.5),
                          gaus(10.0 - d, 2.0, 10.0, 1.5))


def test_gaus_accepts_arrays_and_sign_of_sigma_irrelevant():
    # source: analytic - sigma enters squared
    x = np.linspace(-3, 3, 7)
    assert np.array_equal(gaus(x, 1.0, 0.0, 2.0), gaus(x, 1.0, 0.0, -2.0))


# ------------------------------------------------------ get_line_window
def test_window_half_maximum_hand_computed_line_1():
    # source: analytic - line at 50.3, sigma 2, nominal centre 50.  Half
    # maximum is reached at 50.3 -/+ 2.355: y[48]=0.52 > 0.5, y[47]=0.26 < 0.5
    # -> x0 = 47;  y[52]=0.70 > 0.5, y[53]=0.40 < 0.5 -> x1 = 53.  Count: the
    # starting 5 pixels plus one added per side.
    assert get_line_window(_spectrum(), 50) == (47, 53, 7)


def test_window_half_maximum_hand_computed_line_2():
    # source: analytic - line at 120.4: y[118]=0.487 < 0.5 -> x0 = 118;
    # y[122]=0.73 > 0.5, y[123]=0.43 < 0.5 -> x1 = 123; count 5 + 1 = 6
    assert get_line_window(_spectrum(), 120) == (118, 123, 6)


def test_window_contains_the_line_centre_for_all_lines():
    # source: analytic - the half-maximum window contains the centre
    y = _spectrum()
    for c, _ in LINES:
        x0, x1, n = get_line_window(y, int(round(c)))
        assert x0 < c < x1
        assert n > 0


def test_window_frac_max_widens_window():
    # source: analytic - a lower fraction of the maximum (0.1) gives a wider
    # window than 0.5 for a Gaussian
    y = _spectrum()
    w5 = get_line_window(y, 200, frac_max=0.5)
    w1 = get_line_window(y, 200, frac_max=0.1)
    assert w1[0] < w5[0] and w1[1] > w5[1]


def test_window_strict_accepts_centred_line():
    # source: analytic - the nominal centre is inside the final window
    assert get_line_window(_spectrum(), 50, strict=True) == (47, 53, 7)


def test_window_rejects_line_near_edges():
    # source: analytic - the nominal centre must be at least 2 pixels from the
    # edges; (None, None, 0) is the documented rejection
    y = _spectrum()
    assert get_line_window(y, 1) == (None, None, 0)
    assert get_line_window(y, len(y) - 1) == (None, None, 0)


def test_window_rejects_line_below_threshold():
    # source: analytic - the line peak is 1000; threshold 1500 rejects it
    assert get_line_window(_spectrum(), 50, thresh=1500.0) == (None, None, 0)
    assert get_line_window(_spectrum(), 50, thresh=500.0)[2] > 0


def test_window_rejects_wiggly_profile():
    # source: analytic - a profile that rises again on the low side before
    # reaching half maximum is rejected
    y = np.zeros(60)
    y[30] = 10.0           # peak
    y[29] = 8.0
    y[28] = 6.0
    y[27] = 7.0            # rises again while still above half maximum
    y[31], y[32], y[33] = 8.0, 6.0, 1.0
    assert get_line_window(y, 30) == (None, None, 0)


# ------------------------------------------------------------- findpeaks
# Tolerance (positions, sigma, heights): curve_fit (Levenberg-Marquardt, xtol
# 1.5e-8) is fitting noiseless Gaussian data with a Gaussian model, so it
# converges to the generating parameters; measured error < 1e-11, tolerance
# 1e-6 pixel / 1e-6 relative leaves a large margin without hiding a misplaced
# peak.
def test_findpeaks_recovers_line_centres_widths_and_heights():
    # source: analytic - the lines were generated at LINES with sigma 2
    y = _spectrum()
    pks, sgmn, hts = findpeaks(X, y, 5, 0.001, 50.0)
    assert len(pks) == len(LINES)
    np.testing.assert_allclose(pks, [c for c, _ in LINES], rtol=0, atol=1e-6)
    np.testing.assert_allclose(hts, [a for _, a in LINES], rtol=1e-6, atol=0)
    assert sgmn == pytest.approx(SIGMA, rel=1e-6, abs=0)


def test_findpeaks_amplitude_threshold_drops_faint_lines():
    # source: analytic - amplitudes 1000, 500, 2000, 800; ath = 900 keeps
    # only the lines of 1000 and 2000
    y = _spectrum()
    pks, _, hts = findpeaks(X, y, 5, 0.001, 900.0)
    np.testing.assert_allclose(pks, [50.3, 200.6], rtol=0, atol=1e-6)
    np.testing.assert_allclose(hts, [1000.0, 2000.0], rtol=1e-6, atol=0)


def test_findpeaks_x_axis_in_wavelength_units():
    # source: analytic - the same lines on x = 4000 + 0.5*pixel (Angstrom)
    # are reported in x units: centre = 4000 + 0.5*c, sigma = 0.5*2 = 1.0
    xw = 4000.0 + 0.5 * X
    y = _spectrum()
    pks, sgmn, _ = findpeaks(xw, y, 5, 0.001, 50.0)
    np.testing.assert_allclose(pks, [4000.0 + 0.5 * c for c, _ in LINES],
                               rtol=0, atol=1e-6)
    assert sgmn == pytest.approx(1.0, rel=1e-6, abs=0)


def test_findpeaks_flat_spectrum_returns_nothing(capsys):
    # source: analytic - a constant spectrum has no peaks: empty lists and no
    # mean sigma (the function prints "No peaks found!")
    pks, sgmn, hts = findpeaks(X, np.full_like(X, 100.0), 5, 0.001, 50.0)
    assert pks == [] and hts == [] and sgmn is None
    assert "No peaks found" in capsys.readouterr().out


def test_findpeaks_ignores_lines_within_window_of_edge():
    # source: characterization @ b0e2550 - the scan starts `pkg` (= wid)
    # pixels from either end, so a line centred at 3.3 is not reported while
    # a line at 20.3 is.
    y = _spectrum([(3.3, 1000.0), (20.3, 1000.0), (200.6, 2000.0)])
    pks, _, _ = findpeaks(X, y, 5, 0.001, 50.0)
    np.testing.assert_allclose(pks, [20.3, 200.6], rtol=0, atol=1e-6)
