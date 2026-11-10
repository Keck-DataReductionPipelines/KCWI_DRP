"""ArcOffsets prompt sites run through ``run_primitive`` (task T.4).

Prompt sites (docs/dev/plot_inventory.md, section 3.1):

* P01 (plot_level >= 2), after the figure of bar ``next_bar_to_plot``:
  ``Next? <cr>, <int> to jump, q to quit (? for help): `` (P.4: was
  ``Next? <int> or <cr>, q to quit: ``) -- <cr> shows the next bar, an
  integer n jumps to bar n, q stops the review.
* P02 (plot_level >= 2): ``Next? <cr> (? for help): `` (P.4: was
  ``Next? <cr>: ``) after the bar-offsets figure.

Input: ``context.arcs`` as ExtractArcs leaves it, here 12 synthetic bar
spectra (the same line list shifted by a known number of pixels per bar);
NBARS and REFBAR are set to match.  Which bar each P01 prompt belongs to is
read from the figure title (``... BAR OFFSET for Arc: <n> ...``), recorded by
replacing ``bokeh_plot`` (a no-op here anyway: bokeh_session is None).
save_plot (PNG export through a browser) is replaced by a no-op.

B4 (a bar number lower than or equal to the current one, negative, or >=
the number of bars silently ended the review) and B5 (an answer containing
"q" quit) were pinned as strict xfails (D8); P.4 fixed both by converting
the site to ``kcwi_interactive.ask_next`` and removed the markers.
"""
import re

import numpy as np
import pytest

from kcwidrp.primitives import ArcOffsets as arc_offsets_module
from kcwidrp.primitives.ArcOffsets import ArcOffsets

HEAD = "c0fed5c"
NBARS, REFBAR, NPIX = 12, 5, 300
P01 = "Next? <cr>, <int> to jump, q to quit (? for help): "
P02 = "Next? <cr> (? for help): "
SHIFTS = [(b - REFBAR) % 4 - 1 for b in range(NBARS)]   # -1 ... 2 px


def arcs():
    x = np.arange(NPIX, dtype=float)
    lines = [40.0, 75.0, 120.0, 150.0, 190.0, 230.0, 260.0]
    out = []
    for s in SHIFTS:
        out.append(sum(1000.0 * np.exp(-0.5 * ((x - c - s) / 1.5) ** 2)
                       for c in lines) + 10.0)
    return out


@pytest.fixture
def shown(monkeypatch):
    """Bar numbers of the P01 figures, in display order."""
    bars = []

    def record(p, session):
        m = re.search(r"BAR OFFSET for Arc: (\d+)", p.title.text)
        if m:
            bars.append(int(m.group(1)))
    monkeypatch.setattr(arc_offsets_module, "bokeh_plot", record)
    # PNG export needs a browser; not part of what is tested
    monkeypatch.setattr(arc_offsets_module, "save_plot",
                        lambda *a, **k: None)
    return bars


def run(run_primitive, synthetic_ccddata, plot_level, answers=None):
    return run_primitive(
        ArcOffsets, synthetic_ccddata(IMTYPE="MARC"),
        config_overrides={"plot_level": plot_level, "NBARS": NBARS,
                          "REFBAR": REFBAR},
        context_extra={"arcs": arcs()}, answers=answers)


def test_offsets_are_found_without_prompts(run_primitive, synthetic_ccddata,
                                           shown):
    run(run_primitive, synthetic_ccddata, 0)
    # source: analytic @ c0fed5c - bar b is the reference bar shifted by
    # SHIFTS[b] px; the cross-correlation offset brings it back: -(SHIFTS[b] -
    # SHIFTS[REFBAR])
    expected = [-(s - SHIFTS[REFBAR]) for s in SHIFTS]
    assert list(run_primitive.last.context.bar_offsets) == expected
    assert run_primitive.last.prompts == [] and shown == []


def test_enter_shows_every_bar(run_primitive, synthetic_ccddata, shown):
    # source: characterization @ c0fed5c (inventory P01, P02)
    run(run_primitive, synthetic_ccddata, 2, [""] * NBARS + [""])
    assert shown == list(range(NBARS))
    assert run_primitive.last.prompts == [P01] * NBARS + [P02]


def test_number_jumps_and_q_quits(run_primitive, synthetic_ccddata, shown):
    # source: characterization @ c0fed5c (inventory P01): "" -> next bar,
    # "4" -> bar 4, "" -> bar 5, "q" -> stop
    run(run_primitive, synthetic_ccddata, 2, ["", "4", "", "q", ""])
    assert shown == [0, 1, 4, 5]
    assert run_primitive.last.prompts == [P01] * 4 + [P02]


def test_text_without_q_acts_as_enter(run_primitive, synthetic_ccddata,
                                      shown):
    # source: characterization @ c0fed5c (inventory P01: non-integer text
    # such as 2.5 acts as <cr>)
    run(run_primitive, synthetic_ccddata, 2, ["2.5", "q", ""])
    assert shown == [0, 1]


@pytest.mark.parametrize("bad", ["3", "5", "-1", str(NBARS)])
def test_out_of_range_bar_number_is_asked_again(run_primitive,
                                                synthetic_ccddata, shown,
                                                bad):
    # source: analytic @ c0fed5c - plan section 8.2 (ask_next) - a number
    # outside current+1 ... n_items-1 re-prompts; "7" then jumps to bar 7
    run(run_primitive, synthetic_ccddata, 2, ["5", bad, "7", "q", ""])
    assert shown[0] == 0 and shown[-1] == 7


def test_answer_containing_q_does_not_quit(run_primitive, synthetic_ccddata,
                                           shown):
    # source: analytic @ c0fed5c - plan section 8.2 - only q/quit quit; other
    # text containing a q is asked again
    run(run_primitive, synthetic_ccddata, 2, ["eq", "", "q", ""])
    assert shown[0] == 0 and shown[-1] == 1
