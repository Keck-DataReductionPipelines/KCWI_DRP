"""FitCenter prompt sites run through ``run_primitive`` (task T.4).

(The T.3 file ``test_FitCenter.py`` covers the module function
pascal_shift.)

Prompt sites (docs/plot_inventory.md, section 3.1):

* P08 (plot_level >= 2), after the central-dispersion figure of bar
  ``next_bar_to_plot``: ``Next? <int> or <cr>, q to quit: `` -- <cr> next
  bar, integer n jumps to bar n, q stops the review.
* P09, P10 (plot_level >= 2): ``Next? <cr>: `` after the two summary
  figures.

Input: the arc-chain state that ExtractArcs, ArcOffsets, CalcPrelimDisp and
ReadAtlas leave (``context.arcs``, ``bar_offsets``, ``prelim_disp``; args
``xvals``, ``x0``, ``minrow``, ``maxrow``, ``refwave``, ``reflux``,
``refdisp``, ``offset_wave``), here for 8 synthetic bars sharing one
emission-line spectrum at 1 A/px, and a matching synthetic atlas.

Harness adaptations (they do not change the prompt site): the per-bar fits
run in a ``multiprocessing`` "spawn" pool, which takes seconds to start;
``get_context`` is replaced so the same ``bar_fit_helper`` runs in-process
through the builtin ``map``.  ``bokeh_plot`` is replaced by a recorder of
figure titles (it is a no-op with bokeh_session None), and ``save_plot``
(PNG export through a browser) by a no-op.

Known bug pinned as a strict xfail (D8): B4 (a bar number lower than or
equal to the current one, negative, or >= the number of bars silently ends
the review).
"""
import re
from types import SimpleNamespace

import numpy as np
import pytest

from kcwidrp.primitives import FitCenter as fit_center_module
from kcwidrp.primitives.FitCenter import FitCenter

HEAD = "c0fed5c"
NBARS, NPIX, CWAVE, DISP = 8, 200, 4600.0, 1.0
P08 = "Next? <int> or <cr>, q to quit: "
NEXT = "Next? <cr>: "
LINES = [4520.0, 4547.0, 4580.0, 4601.0, 4633.0, 4655.0, 4690.0]


def spectrum(wave):
    return sum(1000.0 * np.exp(-0.5 * ((wave - c) / 1.2) ** 2)
               for c in LINES) + 5.0


def chain_state():
    xvals = np.arange(NPIX) - NPIX // 2
    refdisp = 0.25
    refwave = np.arange(4400.0, 4800.0, refdisp)
    arcs = [spectrum(CWAVE + DISP * xvals) for _ in range(NBARS)]
    args = dict(xvals=xvals, x0=NPIX // 2, minrow=50, maxrow=150,
                refwave=refwave, reflux=spectrum(refwave), refdisp=refdisp,
                offset_wave=0.0)
    ctx = dict(arcs=arcs, bar_offsets=[0] * NBARS, prelim_disp=DISP)
    return args, ctx


class _InProcessContext:
    """Stands in for get_context('spawn'): Pool().map is the builtin map."""

    def Pool(self):
        return SimpleNamespace(map=lambda f, it: list(map(f, it)),
                               close=lambda: None)


@pytest.fixture
def shown(monkeypatch):
    bars = []

    def record(p, session):
        m = re.search(r"CENTRAL DISPERSION FIT for Bar: (\d+)", p.title.text)
        if m:
            bars.append(int(m.group(1)))
    monkeypatch.setattr(fit_center_module, "get_context",
                        lambda method: _InProcessContext())
    monkeypatch.setattr(fit_center_module, "bokeh_plot", record)
    # PNG export needs a browser; not part of what is tested
    monkeypatch.setattr(fit_center_module, "save_plot",
                        lambda *a, **k: None)
    return bars


def run(run_primitive, synthetic_ccddata, plot_level, answers=None):
    args, ctx = chain_state()
    return run_primitive(
        FitCenter, synthetic_ccddata(IMTYPE="MARC", BCWAVE=CWAVE),
        config_overrides={"plot_level": plot_level, "NBARS": NBARS},
        context_extra=ctx, action_extra=args, answers=answers)


def test_central_fit_without_prompts(run_primitive, synthetic_ccddata,
                                     shown):
    out = run(run_primitive, synthetic_ccddata, 0)
    assert run_primitive.last.prompts == []
    # source: characterization @ c0fed5c - one set of coefficients per bar,
    # and the central wavelength (value at pixel x0) close to CWAVE.
    # Tolerance: the shift is found to one atlas sample (0.25 A) plus the
    # interpolation of the cross-correlation peak; 0.5 A.  (The fit quality
    # itself is not the subject of this file.)
    assert sorted(out.twkcoeff) == list(range(NBARS))
    for b in range(NBARS):
        assert np.polyval(out.twkcoeff[b], NPIX // 2) == pytest.approx(
            CWAVE, abs=0.5)


def test_enter_shows_every_bar(run_primitive, synthetic_ccddata, shown):
    # source: characterization @ c0fed5c (inventory P08-P10)
    run(run_primitive, synthetic_ccddata, 2, [""] * NBARS + ["", ""])
    assert shown == list(range(NBARS))
    assert run_primitive.last.prompts == [P08] * NBARS + [NEXT, NEXT]


def test_number_jumps_and_q_quits(run_primitive, synthetic_ccddata, shown):
    # source: characterization @ c0fed5c (inventory P08)
    run(run_primitive, synthetic_ccddata, 2, ["3", "", "q", "", ""])
    assert shown == [0, 3, 4]


@pytest.mark.xfail(strict=True, raises=pytest.fail.Exception, reason=(
    "B4: a bar number <= the current bar, negative, or >= the number of "
    "bars never matches the loop index (FitCenter.py:339-346), so the "
    "review silently ends and the answers meant for later bars are left "
    "over"))
@pytest.mark.parametrize("bad", ["1", "3", "-2", str(NBARS)])
def test_out_of_range_bar_number_is_asked_again(run_primitive,
                                                synthetic_ccddata, shown,
                                                bad):
    # source: analytic @ c0fed5c - plan section 8.2 (ask_next) - a number
    # outside current+1 ... n_items-1 re-prompts; "6" then jumps to bar 6
    run(run_primitive, synthetic_ccddata, 2, ["3", bad, "6", "q", "", ""])
    assert shown[0] == 0 and shown[-1] == 6
