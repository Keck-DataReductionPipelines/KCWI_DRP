"""SolveArcs prompt sites run through ``run_primitive`` (task T.4).

Prompt sites (docs/plot_inventory.md, section 3.1), at plot_level 2:

* P43 ``Next? <cr>: `` after the residuals figure of bar
  ``next_bar_to_plot``, then
* P44 ``Next? <int> or <cr>, q - quit: `` after its atlas/arc figure --
  <cr> next bar, integer n jumps to bar n, q stops the review;
* P45 ``Next? <cr>: `` after each coefficient figure (5 for the default
  order-4 fit), P46 and P47 ``Next? <cr>: `` after the two statistics
  figures.

Input: the arc-chain state FitCenter, GetAtlasLines and ReadAtlas leave
(``context.arcs``; args ``twkcoeff``, ``at_wave``, ``at_flux``,
``refwave``, ``reflux``, ``atminrow``, ``atmaxrow``, ``atminwave``,
``atmaxwave``, ``dichroic_fraction``), for 6 synthetic bars that share
one line list at 1 A/px; NBARS and REFBAR are set to match.  The bar of each
P44 prompt is read from the figure title (``ATLAS/ARC FIT for Bar <n>``),
recorded by replacing ``bokeh_plot`` (a no-op with bokeh_session None);
``save_plot`` (PNG export through a browser) is replaced by a no-op.

Known bug pinned as a strict xfail (D8): B4 (a bar number lower than or
equal to the current one, negative, or >= the number of bars silently ends
the review).
"""
import re

import numpy as np
import pytest

from kcwidrp.primitives import SolveArcs as solve_arcs_module
from kcwidrp.primitives.SolveArcs import SolveArcs

HEAD = "c0fed5c"
NBARS, NPIX, W0, DISP = 6, 300, 4450.0, 1.0
LINES = [4480.0, 4512.0, 4547.0, 4590.0, 4633.0, 4671.0, 4702.0, 4730.0]
P44 = "Next? <int> or <cr>, q - quit: "
NEXT = "Next? <cr>: "
TAIL = [NEXT] * 7               # 5 coefficient figures + 2 statistics


def chain_state():
    x = np.arange(NPIX, dtype=float)
    arc = sum(1000.0 * np.exp(-0.5 * ((x - (w - W0) / DISP) / 1.5) ** 2)
              for w in LINES) + 5.0
    refwave = np.arange(4400.0, 4800.0, 0.25)
    reflux = sum(np.exp(-0.5 * ((refwave - w) / 1.5) ** 2) for w in LINES)
    args = dict(twkcoeff={b: [DISP, W0] for b in range(NBARS)},
                at_wave=list(LINES), at_flux=[1.0] * len(LINES),
                refwave=refwave, reflux=reflux, atminrow=0,
                atmaxrow=len(refwave), atminwave=W0,
                atmaxwave=W0 + DISP * NPIX, dichroic_fraction=1.0)
    ctx = dict(arcs=[arc.copy() for _ in range(NBARS)])
    return args, ctx


@pytest.fixture
def shown(monkeypatch):
    bars = []

    def record(p, session):
        m = re.search(r"ATLAS/ARC FIT for Bar (\d+)", p.title.text)
        if m:
            bars.append(int(m.group(1)))
    monkeypatch.setattr(solve_arcs_module, "bokeh_plot", record)
    # PNG export needs a browser; not part of what is tested
    monkeypatch.setattr(solve_arcs_module, "save_plot",
                        lambda *a, **k: None)
    return bars


def run(run_primitive, synthetic_ccddata, plot_level, answers=None):
    args, ctx = chain_state()
    return run_primitive(
        SolveArcs, synthetic_ccddata(IMTYPE="MARC"),
        config_overrides={"plot_level": plot_level, "NBARS": NBARS,
                          "REFBAR": 0},
        context_extra=ctx, action_extra=args, answers=answers)


def test_solution_without_prompts(run_primitive, synthetic_ccddata, shown):
    out = run(run_primitive, synthetic_ccddata, 0)
    assert run_primitive.last.prompts == []
    # source: analytic @ c0fed5c - lines placed at pixel (lambda - W0) / DISP
    # are found, so every bar's solution maps pixel x to W0 + DISP * x.
    # Tolerance: Gaussian line centres on 1.5 px lines, ~1e-3 A.
    x = np.arange(NPIX)
    for coeffs in out.fincoeff:
        np.testing.assert_allclose(np.polyval(coeffs, x), W0 + DISP * x,
                                   rtol=0, atol=0.01)
    assert out.av_bar_nls == len(LINES)


def test_enter_shows_every_bar(run_primitive, synthetic_ccddata, shown):
    # source: characterization @ c0fed5c (inventory P43-P47)
    run(run_primitive, synthetic_ccddata, 2, ["", ""] * NBARS + [""] * 7)
    assert shown == list(range(NBARS))
    assert run_primitive.last.prompts == [NEXT, P44] * NBARS + TAIL


def test_number_jumps_and_q_quits(run_primitive, synthetic_ccddata, shown):
    # source: characterization @ c0fed5c (inventory P44)
    run(run_primitive, synthetic_ccddata, 2,
        ["", "3", "", "", "", "q"] + [""] * 7)
    assert shown == [0, 3, 4]


@pytest.mark.xfail(strict=True, raises=pytest.fail.Exception, reason=(
    "B4: a bar number <= the current bar, negative, or >= the number of "
    "bars never matches the loop index (SolveArcs.py:438-445), so the "
    "review silently ends and the answers meant for later bars are left "
    "over"))
@pytest.mark.parametrize("bad", ["1", "2", "-1", str(NBARS)])
def test_out_of_range_bar_number_is_asked_again(run_primitive,
                                                synthetic_ccddata, shown,
                                                bad):
    # source: analytic @ c0fed5c - plan section 8.2 (ask_next) - a number
    # outside current+1 ... n_items-1 re-prompts; "4" then jumps to bar 4
    run(run_primitive, synthetic_ccddata, 2,
        ["", "2", "", bad, "4", "", "q"] + [""] * 7)
    assert shown[0] == 0 and shown[-1] == 4
