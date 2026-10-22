# KCWI DRP upgrade: master implementation plan

Date: 2026-10-22. Status: **plan only, nothing implemented.**
Branch at time of writing: `replace_bokeh` (`dfbcf9d`, based on `master` at `12fd4aa`).

This document coordinates three upgrades into one sequence of tasks, pull requests and human checkpoints:

| Area | Detailed design lives in | This document adds |
|---|---|---|
| **A.** Bokeh → Matplotlib plotting | [bokeh_port_plan.md](bokeh_port_plan.md) (tasks 0.1–5.4), rationale in [bokeh_replacement_proposal.md](bokeh_replacement_proposal.md) | Amendments, model remapping, shared baseline, PR gates |
| **P.** Descriptive, robust interactive prompts | This document (§7) | Full plan |
| **M.** Sparse-matrix cube construction | [matrix_cube_plan.md](matrix_cube_plan.md) (Stages 0–5), copied verbatim from `~/drp/kcwi/KCWI_matrix_cube_plan.md` | Task breakdown, concrete acceptance thresholds, default flip, PR gates |

Where this document and a detailed design disagree, **this document wins** and the difference is listed in that track's "Amendments" section.

---

## 1. How to read this plan

Every task has the same parts:

- **Model**: Opus 5.5 (`claude-opus-5-5`) or Sonnet 5.5 (`claude-sonnet-5-5`), with an effort level. No other models are used (see Q9 in §2).
- **Depends on**: tasks that must be finished first.
- **Do**: what the task produces.
- **Validation**: automated or mechanical checks. The task is not done until each one passes and the result is pasted into the PR description.
- **👤 Human check**: what *you* must do or look at. These cannot be delegated to a model: they need real data, an interactive session, a browser, scientific judgement, or authority. A task with a 👤 block is not done until you have done it and written the result down (in the PR, or in this file for decisions).

Tasks marked **Model: Human** are entirely yours.

All human steps are also collected into one checklist in §10 so nothing is missed.

---

## 2. Decisions record (questions and answers)

Questions asked while writing this plan, and the answers given. Later changes go in as new rows with a date, not as edits to old rows.

| # | Question | Answer (2026-10-22) | Consequence in this plan |
|---|---|---|---|
| Q1 | How should the three areas be sequenced and delivered? | **Sequential PRs.** Bokeh port first, then prompts, then matrix stages 2+. Matrix stages 0–1 run in parallel from the start. One PR per area or stage. | §4 ordering. PR1 (A), PR2 (M0–M1, parallel), PR3a/3b (P), PR4–PR7 (M) |
| Q2 | What should "more descriptive messages" cover? | **Explanatory text + input robustness + fix logic bugs found.** On-plot guidance in the viewer was *not* selected. | §7. Viewer is unchanged by Track P. Bug fixes are a separate PR (3b) because they change outputs |
| Q3 | Which matrix stages are in scope? | **Stages 0–5.** Stage 6 (joint inversion of dithers) is out of scope. | §8 covers M0–M5. Stage 6 stays in matrix_cube_plan.md as future work |
| Q4 | What reference data is on hand? | **Blue + red with standard stars, no nod-and-shuffle.** | S.3 uses BL/BM + KCRM sets. N&S paths are validated on synthetic data only; this is a known gap (§11) |
| Q5 | After validation, what is the shipped default for `cube_method`? | **Flip to `matrix`** in this update. | M.6 makes the flip, gated by Keck maintainer sign-off (G6). Invsens/FluxCalibrate method guard (M.2d) becomes mandatory |
| Q6 | Which covariance products should Stage 5 produce? | **A** (neighbour correlations), **B** (store `W`, opt-in) and **C** (noise-scaling curve in header). | M.5 implements all three. A and C on by default, B behind a config key |
| Q7 | Who makes maintainer-level decisions? | **You plus Keck maintainers.** You decide day-to-day; upstream merges and user-facing defaults need Keck sign-off. | Gates G1–G7 in §4 name what needs Keck sign-off |
| Q8 | Should prompt answers be recorded in the pipeline log? | **Yes, at INFO.** Explanations go to the terminal only. | P.1 helper logs every answer and its effect |
| Q9 | (Instruction, not a question) Which models may be used? | **Opus or Sonnet only.** | Bokeh-plan tasks assigned to Haiku move to Sonnet; the Fable review moves to Opus in a fresh session (§3) |

Decisions taken by the plan author rather than asked, because they have a conventional answer. Overrule any of them by adding a row above.

| # | Decision | Why |
|---|---|---|
| D1 | The matrix plan is copied into `docs/matrix_cube_plan.md` so it is versioned with the code and visible to maintainers. The copy in `~/drp/kcwi/` is no longer the source of truth. | A plan linked from a PR must be readable by reviewers |
| D2 | The reference reductions for Track A (bokeh 0.1) and Track M (Stage 0) are one human session (S.3), not two. | Same data, same runs; halves the human effort |
| D3 | In Track P, the number and order of prompts **for valid answers** does not change. Re-prompting on *invalid* input is allowed. | Keeps the S.3 answer files usable as a regression harness for every later PR |
| D4 | End-of-input (Ctrl-D / closed stdin) at a prompt: by default, log a warning and take the default answer for that and all later prompts. With `KCWI_STRICT_INPUT=1`, raise `EOFError` as today. | Robustness for users; the strict mode keeps the bokeh-plan replay check ("ends with `EOFError` → prompts changed") working |
| D5 | The `warp` path is kept after the default flip, selectable with `cube_method = "warp"`, for at least one release. | Rollback without a code revert; side-by-side comparison for users |
| D6 | New pull requests target `develop` (upstream merges `develop` → `master`, as in PR #207). | Matches existing upstream practice |

---

## 3. Model assignments

| Model | ID | Used for |
|---|---|---|
| Opus 5.5 | `claude-opus-5-5` | Design that many later tasks inherit; numerics (`W`, variance, covariance, WCS); concurrency (plot viewer); domain-heavy text (prompt explanations); first-pass review of real-data results; the final independent review |
| Sonnet 5.5 | `claude-sonnet-5-5` | Well-specified implementation; mechanical call-site conversion; tests; scripts; config, docs and changelog; inventories and sweeps |

Rules carried over from the bokeh plan and adapted to two models:

1. **Tests are written by a different model, or a fresh session, from the code they test**, from the specification and not from the implementation's source. Where Opus implements, Sonnet writes tests. Where Sonnet implements, tests come from Opus or from a separate Sonnet session with no access to the implementation conversation.
2. **The final review (R.2) is Opus 5.5 at high effort in a fresh session** with no access to the implementation conversations. This replaces the Fable 5.1 review in the bokeh plan.
3. Switch models with `/model` in Claude Code, or run a task as a subagent with the named model.
4. These assignments are judgements, not benchmarks. If a Sonnet task produces work that fails validation twice, rerun it with Opus rather than iterating a third time.

Remapping of bokeh-plan tasks:

| Bokeh plan task | Original model | Model in this plan |
|---|---|---|
| 0.2 Plot inventory | Haiku 5.5 | Sonnet 5.5 (low), merged into S.5 |
| 5.2 Sweep | Haiku 5.5 | Sonnet 5.5 (low), merged into R.1 |
| 5.3 Independent review | Fable 5.1 (or Opus) | Opus 5.5 (high), fresh session, merged into R.2 |
| All others | Opus / Sonnet / Human | Unchanged |

---

## 4. Sequencing, pull requests and gates

```
Phase S (shared setup, mostly human)
  S.1 repo hygiene ─┐
  S.2 environment ──┼─► S.3 reference reductions ─► S.5 inventory
  S.4 maintainer Qs ┘
        │
        ├──────────────► Track A: Bokeh port ─────────► PR1 ──► (G1)
        │                                                  │
        └─► M.0 harness ─► M.1a/b maps & pixcoords ─► PR2 ─┤ (G2)  [parallel with A]
                                                           │
                         Track P: prompts (after PR1 merged)
                           P.0–P.4, P.6 ──────────────► PR3a ─► (G3)
                           P.5 bug fixes ─────────────► PR3b ─► (G4)
                                                           │
                         Track M: matrix (after PR3a merged)
                           M.2 W in MakeCube (opt-in) ─► PR4 ─► (G5)
                           M.3 + M.4 λ corr + DAR ─────► PR5
                           M.5 covariance ─────────────► PR6
                           M.6 default flip ───────────► PR7 ─► (G6)
                                                           │
                         Release: R.1 sweep, R.2 review, R.3 clean machine,
                                  R.4 docs branch, R.5 release ─► (G7)
```

Why this order:

- **A before P.** Track P rewrites the code around `bokeh_plot(...)` calls at most prompt sites, and the MakeInvsens limits prompt depends on the new viewer's cursor readout replacing Bokeh tooltips. Doing P first would conflict with every file in A.
- **A before M.2.** `MakeCube` is touched by A (image plot, possibly the `Pool` start method) and rewritten by M.2. Bokeh-plan task 4.2 also requires cubes to be numerically identical to the reference, which only holds while the warp path is the only path.
- **P before M.2** (per Q1). The file overlap is small (one prompt in `MakeCube`), so if Track P stalls, M.2 may start after PR1 instead; record that as a new decision row if it happens.
- **M.0–M.1 in parallel.** They touch `GenerateMaps`, a new `core/cube_matrix.py` and new tests, none of which Track A changes.

Gates. Each needs the named sign-off before the next step starts.

| Gate | When | Who signs | What they agree to |
|---|---|---|---|
| G0 | End of S.4 | You + Keck maintainers | Answers to the maintainer questions in S.4 |
| G1 | PR1 review | Keck maintainers | Bokeh removal, new viewer, config key changes, one viewer tab per process |
| G2 | PR2 review | You (Keck review as normal) | Vectorized GenerateMaps; no behaviour change |
| G3 | PR3a review | Keck maintainers | New prompt wording and behaviour, answers in log |
| G4 | PR3b review | Keck maintainers | Each bug fix and its measured effect on outputs |
| G5 | PR4 review | Keck maintainers | New `cube_method` key, new header keywords and extensions (opt-in) |
| G6 | PR7 review | Keck maintainers | Flipping the default to `matrix`; users must re-derive invsens |
| G7 | Release | You + Keck maintainers | Version number, release notes |

### Task overview

| ID | Task | Model (effort) | Depends on | 👤 |
|---|---|---|---|---|
| S.1 | Repo hygiene and starting tag | Human | – | ✔ |
| S.2 | Environments | Human + Sonnet (low) | – | ✔ |
| S.3 | Reference reductions (shared baseline) | Human | S.1, S.2 | ✔ |
| S.4 | Maintainer questions | Human | – | ✔ |
| S.5 | Plot and prompt inventory | Sonnet (low) | – | spot-check |
| A.* | Bokeh port tasks 0.3, 1.1–5.1b | per bokeh plan | S.3, S.5 | ✔ |
| M.0 | Matrix test harness and regression script | Sonnet (medium) | S.3 | ✔ |
| M.1a | Vectorize GenerateMaps | Sonnet (medium) | M.0 | – |
| M.1b | Pixel-edge coordinate helper | Opus (medium) | M.0 | – |
| M.1t | Tests for M.1a/M.1b | Sonnet (medium), fresh session | M.1a, M.1b | – |
| P.0 | Prompt catalogue and explanation text | Opus (medium) | S.5 | ✔ |
| P.1 | Interactive helper module | Opus (medium) | P.0 | – |
| P.2 | Helper tests | Sonnet (medium) | P.1 | – |
| P.3 | Convert simple pause/next sites | Sonnet (medium) | P.1, P.2, PR1 merged | – |
| P.4 | Convert value-entry and bar-jump sites | Opus (medium) | P.3 | – |
| P.5 | Bug fixes (one commit each) | Opus (medium) + Sonnet tests | P.4 | ✔ |
| P.6 | Interactive acceptance session | Human | P.4 (and P.5 for PR3b) | ✔ |
| P.7 | Changelog and user docs for prompts | Sonnet (low) | P.6 | – |
| M.2a | Build `W` | Opus (high) | M.1b, PR3a merged | – |
| M.2b | Tests for `W` | Sonnet (medium) | M.2a | – |
| M.2c | `cube_method` in MakeCube | Opus (medium) | M.2a, M.2b | – |
| M.2d | Cube-method guard on invsens / FluxCalibrate | Sonnet (medium) | M.2c | – |
| M.2e | Real-data validation of M.2 | Opus (high), then Human | M.2c, M.2d | ✔ |
| M.3 | Fold wavelength corrections into `W` | Opus (medium) | M.2e | ✔ |
| M.4 | Fold DAR into `W` | Opus (high) | M.3 | ✔ |
| M.5a | Covariance mathematics | Opus (high) | M.4 | – |
| M.5b | Covariance I/O, config and tests | Sonnet (medium) | M.5a | – |
| M.5c | Covariance validation | Opus (high), then Human | M.5b | ✔ |
| M.6 | Default flip to `matrix` | Sonnet (low) | M.5c, G0 answer | ✔ |
| M.7 | Performance and memory report | Sonnet (low) | M.6 | – |
| R.1 | Leftover sweep | Sonnet (low) | all | – |
| R.2 | Independent review | Opus (high), fresh session | R.1 | – |
| R.3 | Clean-machine acceptance | Human | R.2 | ✔ |
| R.4 | Documentation branch | Sonnet (medium) | R.3 | ✔ |
| R.5 | Release | Human + Keck | R.4 | ✔ |

---

## 5. Phase S: shared setup

### S.1 Repo hygiene — Human

Supersedes bokeh-plan task 0.5, which is mostly done: the working tree on `replace_bokeh` is clean, and the earlier mode flips (644 → 755) plus an uncommitted geckodriver fallback in `kcwi_plotting.py` were stashed as `stash@{0}` ("On master: pre-replace_bokeh: …").

👤 Human check:
1. Decide what happens to `stash@{0}`. The mode flips should be dropped (`git config core.filemode false` if they keep coming back). The geckodriver fallback becomes moot once Track A removes Firefox; drop it unless a release goes out before PR1.
2. Decide whether branch `geckodriver_prefer_non_snap` (`d0b4d7b`) is merged upstream before PR1. It changes only `kcwidrp/core/kcwi_plotting.py`, which Track A rewrites.
3. Tag the starting point: `git tag upgrade-start dfbcf9d` (or the commit you settle on).

Validation:
- `git status --short` is empty on `replace_bokeh`.
- `git diff --summary upgrade-start` shows no `mode change` lines; this is re-checked in R.2.

### S.2 Environments — Human, with Sonnet (low) for the scripted parts

The shell environment active at time of writing (`pypeit`, Python 3.13) does **not** have `bokeh` installed, so it cannot run the current pipeline. The baseline must run in the `kcwidrp` conda environment.

Do:
1. Baseline environment: `conda activate kcwidrp`, then `pip install -e .` from `upgrade-start`. Record `pip freeze > docs/baseline/env_baseline.txt`. Add a `.gitignore` rule so that only text files (`*.md`, `*.txt`, `*.json`) under `docs/baseline/` can be committed; the FITS outputs and plots from S.3 stay local.
2. Development environment for Tracks A, P and M (created after PR1 changes dependencies): a fresh venv or conda env from the branch's `pyproject.toml`.

Validation:
- `python -c "import kcwidrp, bokeh, matplotlib; print(bokeh.__version__, matplotlib.__version__)"` succeeds in the baseline env and the versions are recorded.
- `reduce_kcwi --help` exits 0.
- `pytest kcwidrp/tests` result on unmodified code is recorded (pass/fail per test), so later failures can be told apart from pre-existing ones.

👤 Human check: confirm the baseline env can open the Bokeh browser tab on your Mac (needed for S.3 screenshots). If it cannot, S.3 must run on another machine; record which.

### S.3 Reference reductions — Human

Replaces bokeh-plan task 0.1 **and** matrix Stage 0's data capture (decision D2). Follow bokeh-plan 0.1 steps 1–4 exactly, with these additions:

- **Datasets:** one blue night (BL or BM) and one red (KCRM) night, each with bias, flat, arc, continuum bars, a standard star and at least one science frame with a point source. No nod-and-shuffle (Q4).
- **Passes per channel:** `plot_level = 0`, the usual `plot_level` (1), and `plot_level = 3` with answers captured via `tee answers_<channel>_3.txt`. Give at least one non-default *valid* answer at every value prompt (ReadAtlas offset, MakeInvsens limits / one line mask / fit order, MakeMasterFlat peak if it appears). Do **not** type invalid input in this pass: some invalid inputs crash the current code (§12), which would make the replay unusable.
- **Archive per channel**, under `docs/baseline/<channel>/` (data excluded by the S.2 `.gitignore` rule): `plots_ref_<channel>_<level>/`; `answers_<channel>_3.txt`; the full `redux/` directory including `*_geom.pkl`, `*_wavemap`, `*_slicemap`, `*_posmap`, `icube`, `icubew`, `icubed`, `icubes`, `invsens`, the `.proc` table; and the pipeline log for each pass.
- **Timings:** per-primitive wall time from the log for `GenerateMaps`, `MakeCube`, `WavelengthCorrections`, `CorrectDar` (matrix plan §3 baseline). Machine, CPU and core count recorded.
- **Screenshots** of the live-only plots during the level-3 pass, as bokeh 0.1 requires.

Validation (in addition to bokeh-plan 0.1's list):
- The level-3 replay (`reduce_kcwi … < answers_<channel>_3.txt` from a clean output dir) completes unattended and produces identical cubes (`numpy.array_equal` on every extension) and the same PNG names.
- Every output file listed above exists for both channels.
- Search both level-1 logs for `Picking the highest peak`. Record whether it occurs: it decides how bug B1 (§12) can be validated.

👤 Human check: this whole task is yours. Write a short `docs/baseline/README.md` with dataset IDs (KOA file names), dates, config used, command lines, machine, and anything unusual during the runs.

### S.4 Maintainer questions — Human

Replaces bokeh-plan task 0.4 and adds questions from Tracks P and M. Ask the Keck maintainers and record the answers as new rows in §2 (G0).

From the bokeh plan (unchanged):
1. Does anything downstream (RTI, KOA ingestion) read `plots/*.png`, and does it depend on more than name and format?
2. How do observers view plots when the pipeline runs remotely (VNC, SSH forward to 5006)?
3. Are several pipeline processes watched at once? Is one tab per process acceptable?
4. Should the browser tab auto-open on summit and RTI hosts?
5. Who can merge to the `documentation` branch, and how is ReadTheDocs rebuilt?

New for Track P:
6. Are the prompts driven by any script or wrapper at Keck (answers piped in, `expect`, etc.)? If so, the exact prompt strings are an interface and the new wording must be agreed first.
7. Is logging user answers at INFO acceptable for KOA/RTI log ingestion?

New for Track M:
8. Is flipping the default `cube_method` to `matrix` in this release acceptable (Q5)? Under what evidence?
9. Can new FITS extensions be added to `icube*` outputs (coverage cube, covariance planes) without breaking KOA ingestion or known downstream tools (CWITools, QFitsView workflows)? Or must they go in separate files?
10. Are new header keywords acceptable: `CUBEMETH` (warp|matrix), `WAVCORR`, and the noise-scaling-curve keywords from M.5 (names to be proposed in M.5a)?
11. Is a ~160 MB `W` file per frame acceptable as an opt-in product, and where should it be written?
12. Are invsens files from previous DRP versions reused across reductions in practice (e.g. a library of standards)? This decides how strict the M.2d guard should be.

Validation: every question has an answer row in §2 before the task that depends on it starts (questions 1–5 before A.3.1; 6–7 before P.0; 8–12 before M.2c).

👤 Human check: this whole task is yours.

### S.5 Plot and prompt inventory — Sonnet 5.5 (low)

Bokeh-plan task 0.2 as written (file `docs/plot_inventory.md`, 48 figure rows, feature summary), with the prompt table extended for Track P. For every `input(` call, one row with:

- file and line; the `plot_level` guard; the exact prompt string;
- which figure is on screen when it is asked (cross-reference to the figure row);
- accepted answers and what each does (next, quit, jump to bar *n*, new value);
- the default (what `<cr>` does);
- **whether the answer changes the reduction**, and which output it changes;
- the matching `time.sleep(plot_pause)` branch, if any;
- any input that raises or misbehaves (to be checked against §12).

Validation:
- 52 `input(` lines across 17 primitives, plus 8 in `scripts/reduce_kcwi.py` and `scripts/kcwi_rti.py` (counts from `grep` on `dfbcf9d`). `MakeMasterFlat.py:471` and `:473` are one logical prompt and may share a row.
- 36 `time.sleep(...plot_pause)` sites listed.
- The bokeh-plan 0.2 checks (48 rows, 49 `bokeh_plot(`, 20 `save_plot(`, the easy-to-miss features).

👤 Human check: spot-check five figure rows and five prompt rows against the source, including the "changes the reduction" column for ReadAtlas, MakeInvsens and MakeMasterFlat.

---

## 6. Track A: Bokeh → Matplotlib

Execute [bokeh_port_plan.md](bokeh_port_plan.md) tasks 0.3, 1.1, 1.2, 2.1, 2.2, 3.1, 3.2, 4.1, 4.2, 4.3, 5.1, 5.1b as written, with these amendments.

### Amendments

1. **Shared tasks.** Bokeh 0.1 → S.3; 0.2 → S.5; 0.4 → S.4; 0.5 → S.1. Bokeh 5.2 → R.1; 5.3 → R.2; 5.4 → R.3 (which also covers Track M).
2. **Models.** As in the remap table in §3.
3. **No prompt changes in Track A.** Prompt strings, order and guards stay byte-identical; Track P owns them. The one exception is already in the bokeh plan: the MakeInvsens "Hover the cursor…" instruction must be true of the new viewer's cursor readout (bokeh 3.1 verification, last bullet).
4. **Minimal MakeCube changes.** In bokeh 2.2, change only the imports, the `image` call support and (if 0.3 says so) the `Pool` start method. Do not refactor `make_cube_helper`; M.2c rewrites that code.
5. **The coverage scan (bokeh 1.2) becomes a permanent gate.** It stays in the test suite and runs on PR3, PR4–PR7. Any new plot added by Tracks P or M must use only wrapper-supported features or extend the wrapper with a test in the same PR.
6. **Matplotlib is the only plotting library in new code.** Validation scripts in Track M may use `matplotlib.pyplot` because they are not part of the package; package code follows wrapper rule 4 (no `pyplot`).
7. **PR1** is opened against `develop` after bokeh 4.3 sign-off; 5.1 and 5.1b may land in PR1 or as a follow-up, at the maintainers' preference (G1).

### Human checks in Track A (from the bokeh plan, collected here)

👤 Human check:
- **0.3 prototype:** the prototype's pass/fail list must be run on both macOS and Ubuntu 24.04. The browser-facing criteria (zoom, pan, cursor readout, tab close/reopen, two instances on two ports) need you at a browser.
- **3.1:** the interactive MakeInvsens session (cursor readout precise enough to type wavelength limits).
- **4.3:** review the side-by-side comparison page and the 4.2 table, step through one interactive `plot_level = 3` reduction, write the sign-off with accepted differences.
- **5.1b:** whoever S.4 identified merges the documentation branch.

---

## 7. Track P: descriptive and robust interactive prompts

### 7.1 Goals and constraints

Goals (Q2, Q8):
1. Before each prompt, the user is told **what the plot shows, what to look for, what each answer does, and what the default is.**
2. Input is validated consistently: bad input re-prompts with a reason, `?` re-prints the full explanation, quitting is unambiguous, end-of-input is handled (D4).
3. Every answer that changes the reduction is logged at INFO with its effect.
4. Logic bugs at prompt sites (§12) are fixed, each in its own commit, with its effect on outputs measured.

Constraints:
- **D3:** for valid answers, the number and order of prompts is unchanged, so `answers_<channel>_3.txt` from S.3 replays on every commit of PR3a and produces identical outputs.
- **No behaviour change in PR3a.** Everything that changes a reduction result is in PR3b.
- **No plotting changes.** On-plot guidance was not selected (Q2). The helper must not import the plotting modules.

### 7.2 Design (to be fixed in P.1)

A new module `kcwidrp/core/kcwi_interactive.py`. Proposed API; P.1 may refine names but not semantics:

| Function | Replaces | Behaviour |
|---|---|---|
| `pause(config, logger, what=None)` | The 36 `if plot_level >= 2: input("Next? <cr>: ") else: time.sleep(plot_pause)` blocks | Same branching. When interactive, prints a one-line description (`what`) then `Next? <cr> (? for help)` |
| `ask_next(config, logger, what, help, n_items=None, current=None)` | `Next? <cr>, q to quit` and `Next? <int> or <cr>, q to quit` loops (ArcOffsets, SolveArcs, FitCenter, ExtractArcs, FindBars, MakeCube, MakeInvsens slice loop) | Returns `NEXT`, `QUIT` or a validated item number. A number outside `current+1 … n_items-1` re-prompts with the valid range (today it silently ends the review: bug B4) |
| `ask_value(config, logger, prompt, help, parse, validate, default, effect)` | Value prompts: ReadAtlas offset, MakeInvsens limits / line mask / fit order, MakeMasterFlat peak, NandshuffSubtractSky scale | Shows current value and default; re-prompts on parse or validation error with the reason; logs the accepted value with `effect` |
| `ask_yes_quit(...)` | `Proceed? <cr>=yes or Q=quit` in `reduce_kcwi.py` and `kcwi_rti.py` | Same choices; exact match on `q`/`quit` |

Rules:
- **Quit** is an exact, case-insensitive match on `q` or `quit`. Today `'Q' in q.upper()` quits on any input containing a `q`.
- **Help:** `?` or `help` prints the long `help` text and asks again; it never counts as an answer.
- **End of input (D4):** default: log one WARNING naming the prompt, take the default, and treat every later prompt in the process as answered by its default (no repeated warnings). With environment variable `KCWI_STRICT_INPUT=1`: raise `EOFError` as today. All replay checks in this plan run with `KCWI_STRICT_INPUT=1`.
- **Ctrl-C:** unchanged. The framework's SIGINT handler is not touched (see bokeh plan 3.1).
- **Logging (Q8):** explanations are printed to stdout and not logged. Each *answer* is logged at INFO in the form `<Primitive>: <effect>`, e.g. `ReadAtlas: atlas offset set by user to 12 px (calculated 10 px)` or `MakeInvsens: user quit slice review at slice 7`. Pure "next" presses are logged at DEBUG.
- **Text location:** short `what` and long `help` strings live at the call site, as keyword arguments, so a reader of the primitive sees what the user is told. No central text registry.
- **Text template** for `help`: *What you see* (plot title, which colour is which) / *What to look for* (good vs. bad) / *Your options* (each answer and its effect on the reduction) / *Default*.

### P.0 Prompt catalogue and explanation text — Opus 5.5 (medium)

Depends on: S.5. Produces `docs/interactive_prompts.md`: one section per logical prompt, using the S.5 prompt table, with the drafted `what` and `help` text following the template. This is domain writing (what a good arc fit or atlas cross-correlation looks like), hence Opus.

Validation:
- One section per logical prompt in S.5; none missing (count check).
- Every "changes the reduction" prompt states which output is affected.
- Every colour named in a help text matches the colour in the figure code at that site.

👤 Human check: **read every explanation for scientific correctness and clarity.** This is the most important review in Track P: a model can describe the code correctly and the science wrongly. Mark each section approved or corrected. If S.4 Q6 said prompt strings are an interface, have the maintainers approve the wording too.

### P.1 Interactive helper module — Opus 5.5 (medium)

Depends on: P.0. Implement `kcwidrp/core/kcwi_interactive.py` per §7.2. Opus because every one of the 60 call sites inherits its semantics, and the EOF/strict/replay interaction is easy to get subtly wrong.

Validation:
- `python -c "import kcwidrp.core.kcwi_interactive, sys; assert not any(m.startswith(('matplotlib','bokeh')) for m in sys.modules)"` passes.
- P.2 tests pass.
- Docstrings state the EOF and strict-mode behaviour.

### P.2 Helper tests — Sonnet 5.5 (medium)

Depends on: P.1. Written from §7.2, not from P.1's source. Use `monkeypatch` on `builtins.input` with scripted answers and `caplog` for logging.

Cover: each function with valid input; invalid then valid (re-prompt, reason printed, one answer logged); `?` then valid (help printed, not counted); exact `q`/`quit` vs. a value containing `q`; out-of-range item numbers; EOF default mode (one warning, default used, later prompts silent defaults); EOF with `KCWI_STRICT_INPUT=1` raises `EOFError`; `plot_level < 2` path calls `time.sleep(plot_pause)` and never `input`; INFO log line present for value answers and absent for "next".

Validation:
- At least three deliberate mutations of the helper (substring quit restored; EOF swallowed in strict mode; out-of-range accepted) each make a test fail.
- Suite runs in under 5 s with no terminal attached.

### P.3 Convert simple pause/next sites — Sonnet 5.5 (medium)

Depends on: P.1, P.2, PR1 merged. Convert every `pause` and `ask_next` site **except** those listed in P.4. Add each site's `what`/`help` text from P.0.

Validation:
- New test `kcwidrp/tests/test_no_bare_input.py`: an AST scan that fails if `input(` or `time.sleep(...plot_pause...)` appears in `kcwidrp/primitives` outside the helper. (It starts with the P.4 sites on an allow-list, which P.4 empties.)
- `pytest kcwidrp/tests` passes, including the bokeh coverage scan.
- **Replay:** with `KCWI_STRICT_INPUT=1`, replay S.3's `answers_<channel>_3.txt` for both channels: completes, consumes the file exactly, and every cube and the `.proc` table are `array_equal` to the PR1 output. PNGs are pixel-identical to PR1's.
- `git diff` per file touches only prompt blocks and imports.

### P.4 Convert value-entry and bar-jump sites — Opus 5.5 (medium)

Depends on: P.3. Sites: `ReadAtlas` offset; `MakeInvsens` limits, line mask, fit order and slice loop; `MakeMasterFlat` peak; `NandshuffSubtractSky` scale (×2); the bar-jump prompts in `ArcOffsets`, `SolveArcs`, `FitCenter`; the `Proceed?` prompts in `reduce_kcwi.py` and `kcwi_rti.py`. Opus because each site is interleaved with reduction logic, and a validator that is stricter than today changes which answers are accepted.

Rules: **behaviour-preserving for every input accepted today**. Where today's behaviour is a bug (§12), keep it in P.4 and fix it in P.5. Re-prompting on input that today crashes the primitive (e.g. B2) is allowed in P.4, because no valid run can depend on a crash.

Validation:
- As P.3, and the AST allow-list is now empty.
- Synthetic tests for each value site: the value reaching the reduction code is identical to today's for a table of valid inputs (from P.0).
- NandshuffSubtractSky: test on a synthetic non-standard N&S configuration (no real N&S data, Q4).

### P.5 Bug fixes — Opus 5.5 (medium) for fixes, Sonnet 5.5 (medium) for the failing tests

Depends on: P.4. One commit per bug in §12. For each bug: Sonnet first writes a test that fails on the current code, from the bug description only; Opus then fixes it; the test passes.

Validation:
- Each commit message states the bug, the fix, and the measured effect on the S.3 reductions ("no change in outputs" or "invsens changes by X % at λ Y").
- After each output-changing fix, the S.3 replay is rerun and the diff explained.

👤 Human check: approve or reject each fix individually. For B1 (MakeMasterFlat peak), look at the ledge and derivative plots for the affected flat and confirm the new peak is the right one. If S.3 never hit the buggy branch, confirm you accept validation on the synthetic test alone.

### P.6 Interactive acceptance session — Human

Depends on: P.4 (and again after P.5 for PR3b). Run one `plot_level = 3` reduction per channel in the development env, with the viewer open.

👤 Human check, at every prompt:
- Could a new user answer this prompt from the printed text alone, without reading code or the manual? Note any that fail.
- Try `?` at least once per prompt type; try one invalid input per value prompt; try a bar number below the current one and above the last.
- Try Ctrl-D at one prompt (expect: warning, defaults for the rest, reduction finishes).
- Afterwards, read the log: every value answer is there with its effect; no explanations are in it.
- Capture this session's input as `answers_<channel>_3_v2.txt` (it includes invalid inputs and `?`), and confirm it replays with `KCWI_STRICT_INPUT=1`. It becomes the regression input for all later PRs.

### P.7 Changelog and docs — Sonnet 5.5 (low)

Changelog entries for PR3a (new prompt behaviour, logging, EOF handling, `KCWI_STRICT_INPUT`) and PR3b (each fix and its effect). A short "Interactive mode" section for the documentation branch, held for R.4.

Validation: every behaviour change listed in §7.2 rules and every §12 fix appears in the changelog.

---

## 8. Track M: sparse-matrix cube construction

Executes [matrix_cube_plan.md](matrix_cube_plan.md) Stages 0–5.

### Amendments

1. **Stage 0 data** come from S.3 (D2). No nod-and-shuffle real data (Q4): N&S cube building is validated on a synthetic frame only.
2. **Stage 5** implements options A, B (opt-in) and C (Q6). Stage 6 is out of scope (Q3).
3. **Default flips to `matrix`** at the end of the track (Q5), as a separate PR (M.6) behind gate G6.
4. **Method guard.** Because the default changes, an invsens built with one method must not silently calibrate a cube built with the other (matrix plan §6, "never mix methods"). M.2d adds this.
5. **Concrete thresholds.** The matrix plan's acceptance criteria are qualitative ("small fraction of a pixel"). Proposed numbers are given below; 👤 you confirm or change them in M.0 before any validation runs, so the thresholds are not tuned to the results.
6. **Plots.** Validation figures are produced by the regression script (outside the package), not by the pipeline. If a pipeline plot is added (e.g. a coverage map at `plot_level = 3`), it uses the wrapper and passes the coverage scan (Track A amendment 5).
7. **Line references** in the matrix plan were taken at `d0b4d7b`; the primitives are identical at `dfbcf9d` (only `kcwi_plotting.py` differs), so they remain valid until Track A merges.

### M.0 Harness and regression script — Sonnet 5.5 (medium)

Depends on: S.3. Matrix plan Stage 0 (tests and metrics only; data from S.3).

Do:
- `kcwidrp/tests/test_cube_matrix.py`: a fast synthetic fixture: a small detector, two or three slices, and **affine** per-slice transforms built as `skimage` `PolynomialTransform` objects so the exact answer is known analytically. Initially it tests only the fixture itself.
- `kcwidrp/tests/regression/cube_regression.py` (not collected by pytest): given two `redux/` directories, compute and write a report (Markdown + PNGs via matplotlib) of the metrics in matrix plan Stage 0: arc-line centroid and width per slice vs `geom['avwvsig']`; contbar positions along slice; standard-star total flux and invsens ratio; σ calibration in sky regions (single spaxel and binned apertures of 1, 2×2, 3×3, 5×5 spaxels); runtime per primitive.
- `docs/baseline/regions_<channel>.json`: sky-only regions and the standard-star aperture used by the σ and flux metrics.

Proposed thresholds (👤 confirm before M.2e):

| Metric | Proposed acceptance (matrix vs. warp) |
|---|---|
| Arc line centroid | mean |Δ| ≤ 0.05 px, max |Δ| ≤ 0.2 px, per slice |
| Arc line width | matrix width ≤ warp width × 1.05 |
| Contbar position along slice | mean |Δ| ≤ 0.05 px |
| Standard-star total flux in aperture | within 1 % over WAVGOOD |
| Invsens ratio (re-derived with each method) | within 2 % over WAVGOOD, no structure > 1 % on scales < 50 Å |
| σ calibration, single spaxel | reported/empirical ratio closer to 1 than warp's |
| σ calibration, 5×5 binned | with covariance (M.5), ratio within 10 % of 1 |
| MakeCube wall time | ≤ warp's on the same machine |
| Peak memory of MakeCube | ≤ 2 GB for a 2×2-binned frame |

Validation:
- Running the regression script with both arguments pointing at the same S.3 `redux/` gives zero difference on every metric.
- Running it on blue vs. red completes (red/blue are different data, so it only checks the script handles both geometries).

👤 Human check:
- **Choose the sky regions and the standard-star aperture** in the S.3 cubes (ds9 or QFitsView) and write them to `regions_<channel>.json`. Regions must be free of the target and of bright sky-line residual artefacts.
- **Confirm or change the thresholds** in the table and commit them before M.2e.

### M.1a Vectorize GenerateMaps — Sonnet 5.5 (medium)

Matrix plan Stage 1, first bullet. Replace the per-pixel loops in `GenerateMaps._perform` with one transform call per slice on all pixels, keeping the pixel-centre convention and the `0 <= x <= xsize` test.

Validation:
- `wavemap`, `slicemap`, `posmap` are `array_equal` (NaN positions included) to S.3 for both channels.
- `MakeMasterSky` output on S.3 data is `array_equal` to the baseline.
- Wall time of GenerateMaps before/after recorded in the PR.

### M.1b Pixel-edge coordinate helper — Opus 5.5 (medium)

Matrix plan Stage 1, second bullet. `kcwidrp/core/cube_matrix.py`: from `geom`, return per on-slice pixel `u_lo, u_hi, lam_lo, lam_hi`, slice id and flat pixel index; cache as `*_pixcoords.npz` next to `*_geom.pkl`. Opus because pixel-edge sign conventions, the `xl0` offset and the red-channel rotation are where off-by-half-pixel errors originate, and every later stage inherits them.

Validation:
- Midpoints of `(u_lo, u_hi)` and `(lam_lo, lam_hi)` reproduce M.1a's `posmap` and `wavemap` at pixel centres to 1e-6.
- `u_hi > u_lo` and `lam_hi > lam_lo` everywhere (or the reversed convention is documented and consistent).
- The cache is invalidated when `geom.pkl` changes (keyed on its hash), and a stale cache is never read.
- Synthetic affine fixture: edges exactly match the analytic answer.

### M.1t Tests for M.1a and M.1b — Sonnet 5.5 (medium), fresh session

From the matrix plan and the M.1b validation list, not from the source. Includes a red-channel geometry case (matrix plan open question 1).

PR2 = M.0 + M.1a + M.1b + M.1t, gate G2. No change to cube outputs.

### M.2a Build `W` — Opus 5.5 (high)

Depends on: M.1b, PR3a merged. `cube_matrix.build_W(geom, pixcoords, gpm, out_grid, mode='overlap'|'bilinear')` per matrix plan §4: vectorized, no Python loops over pixels or voxels; per-slice build then `sp.vstack`/`block_diag`; `ncand` per axis computed from data; mask then row-normalize; coverage kept; float32 values, int32 indices; flag bit-plane OR via the boolean structure of the unmasked overlap matrix.

Validation:
- `mode='bilinear'` vs `skimage.transform.warp(order=1)` on the synthetic fixture and on one real S.3 slice: max relative difference ≤ 1e-6 where coverage is full.
- Constant image → constant cube where coverage is full.
- Unit impulse at one pixel → nonzero voxels form the pixel's footprint; for an unmasked pixel, the unnormalized column sums to 1.
- White-noise Monte Carlo (≥ 500 realizations on the synthetic fixture): `sqrt((W∘W) @ σ²)` matches the empirical per-voxel standard deviation to within 3 %.
- A masked pixel contributes nothing to any voxel; a voxel whose only contributors are masked has coverage 0 and is masked.
- Build time for a 2×2-binned real frame recorded; peak memory (`tracemalloc` or `/usr/bin/time -l`) below the M.0 threshold.
- A test asserts `build_W` does not loop over pixels: build time on a frame 4× larger grows by no more than ~5×.

### M.2b Tests for `W` — Sonnet 5.5 (medium)

Independent tests of the M.2a validation list, plus: `ncand` large enough for the largest footprint/voxel ratio (a deliberately coarse output grid); flags of two different bits from two contributing pixels both appear in the voxel; float32/int32 dtypes.

### M.2c `cube_method` in MakeCube — Opus 5.5 (medium)

Depends on: M.2a, M.2b, S.4 answers 8–10. Matrix plan Stage 2: config key `cube_method = "warp" | "matrix"` (default `warp` in this PR); in matrix mode, every MakeCube output (data, uncertainty, mask, flags, noskysub, obj, sky, delta-wave) via `W`; all four MakeCube events (`arc_make_cube`, `flat_make_cube`, `object_make_cube`, `nandshuff_make_cube`) use the same method; red 180° rotation preserved; WCS/header code unchanged; new header keyword `CUBEMETH`. Coverage cube written as agreed in S.4 Q9.

Validation:
- `cube_method = "warp"`: every output `array_equal` to PR3a's on S.3 data (the switch alone changes nothing).
- `cube_method = "matrix"`: full reduction of both channels completes; every output file exists with the same HDU names, shapes and dtypes as the warp output (plus the agreed new extension).
- Synthetic N&S frame runs through `nandshuff_make_cube` in matrix mode.
- `CUBEMETH` is present in every cube and propagates to `icubew`, `icubed`, `icubes`, `invsens`.
- `pytest kcwidrp/tests` passes, coverage scan included.

### M.2d Cube-method guard — Sonnet 5.5 (medium)

Matrix plan §6: invsens and science cubes must use the same method. In `FluxCalibrate`, compare `CUBEMETH` of the invsens and the cube. A missing keyword means `warp` (files from older versions). On mismatch: behaviour as decided in S.4 Q12; the proposed default is to log an ERROR naming both files and the fix ("re-derive invsens with cube_method = …") and skip flux calibration for that frame, rather than calibrate wrongly.

Validation: tests for matching, mismatching and missing-keyword cases; a real S.3 run with a warp invsens and a matrix cube triggers the guard.

### M.2e Real-data validation of M.2 — Opus 5.5 (high), then Human

Depends on: M.2c, M.2d. Run both channels with `cube_method = "matrix"`, re-derive invsens with matrix, run the regression script against S.3. Opus writes a first-pass report: each threshold pass/fail with numbers; any failure investigated to a cause (geometry convention, `ncand`, masking) before proposing a fix.

Validation: every M.0 threshold has a recorded value; failures have a stated cause and either a fix (with a test) or an explicit request to relax the threshold.

👤 Human check:
- Blink warp and matrix `icube` in ds9/QFitsView for one science frame and one standard per channel. Look for: slice-edge artefacts or stripes, ringing around cosmic rays and bad pixels (should be gone), holes where coverage should be full, red channel orientation matching blue.
- Look at the coverage cube: smooth, ~1 inside slices, falling to 0 only at slice and wavelength edges.
- Read the regression report and sign off thresholds, or record why a relaxed threshold is accepted.

PR4 = M.2a–M.2e, gate G5. Default still `warp`.

### M.3 Fold wavelength corrections into `W` — Opus 5.5 (medium)

Matrix plan Stage 3. Compute radial-velocity and air→vacuum corrections from the header before building `W`, apply them to `lam_lo/hi`, set `WAVCORR = T`; `WavelengthCorrections` skips resampling when `WAVCORR` is set but still writes `VCORR`, `VCORRTYP`, `CTYPE3` etc. and still produces `icubew`. Reuse `a2v_conversion` and the velocity computation in `heliocentric`; do not duplicate the formulas.

Validation:
- The velocity correction computed in the helper equals the one `WavelengthCorrections.heliocentric` computes for the same header, to 1e-9 km/s, for both `heliocentric` and `barycentric` modes, and `none` gives zero.
- Sky-line centroids in matrix `icubew` vs. warp `icubew`: mean |Δ| ≤ 0.05 px; arc lines likewise.
- Line widths in matrix `icubew` ≤ warp × 1.05 (fewer resamplings).
- Header diff of `icubew` (warp vs. matrix) shows only expected keywords (`CUBEMETH`, `WAVCORR`).
- `air_to_vacuum = False` and `radial_velocity_correction = "none"` each run and give the expected uncorrected grid.

👤 Human check: plot one bright sky line's profile (warp vs. matrix `icubew`, same spaxel) and confirm the matrix profile is at the same wavelength and not broader.

### M.4 Fold DAR into `W` — Opus 5.5 (high)

Matrix plan Stage 4. Per-λ along-slice and cross-slice offsets from `CorrectDar.atm_disper` and the same projection; padded output grid per grating as `CorrectDar` does; `DARCOR = T`; `CorrectDar` skips shifts when set but keeps header updates, CRPIX padding (matrix plan open question 3), padding flags (value 128) and the `icubed` output. Opus at high because the WCS and padding must match exactly and cross-slice weights are new behaviour.

Validation:
- Standard-star centroid vs. λ (along and across slice): RMS about the mean ≤ the warp path's, both channels.
- RA/Dec of the standard-star centroid at three wavelengths (blue end, centre, red end) from matrix `icubed` WCS matches warp `icubed` within 0.1″.
- Padding region carries flag 128 and the same extent as warp.
- Standard-star flux within 1 %; invsens re-derived with matrix vs. warp within the M.0 threshold.
- CRPIX1/2 in matrix `icubed` equal warp's.

👤 Human check: in ds9, step through wavelength planes of the matrix `icubed` standard star and confirm it does not drift. Compare the two invsens curves on one plot and accept or reject the difference.

PR5 = M.3 + M.4 (may be split if review prefers).

### M.5a Covariance mathematics — Opus 5.5 (high)

Matrix plan Stage 5, options A, B, C (Q6). Per slice, `W diag(var) Wᵀ`; option A extracts correlation coefficients with ±1 neighbours along each axis (6 planes); option B serializes `W` (and the per-pixel variance needed to use it) as `.npz`; option C fits a binned-aperture noise-scaling curve and proposes header keyword names (for S.4 Q10). Includes a documented formula for users of option B: `var(aᵀcube) = (Wᵀa)ᵀ diag(var) (Wᵀa)`.

Validation:
- Monte Carlo on the synthetic fixture: for apertures of 1, 2×2, 3×3, 5×5 spaxels, predicted variance from B equals empirical to within MC error (≤ 3 % at 1000 realizations); from A's neighbour approximation, within 10 %; C's curve within 10 % of empirical at each aperture size.
- Covariance build time per frame recorded; does not exceed 10 s on the S.3 machine (matrix plan measured ~6 s).

### M.5b Covariance I/O, config and tests — Sonnet 5.5 (medium)

Config keys (names proposed, 👤 confirmed in G5/G6 review): `cube_covariance = True` (A and C), `cube_save_matrix = False` (B). Write A as agreed in S.4 Q9; C as header keywords; B as `.npz` next to the cube. Tests independent of M.5a's source.

Validation:
- Files open with `astropy.io.fits` and `scipy.sparse.load_npz`; a round-trip test reconstructs B's matrix exactly.
- With both keys off, outputs are identical to M.4's.
- File-size increase per cube recorded for A and B.

### M.5c Covariance validation — Opus 5.5 (high), then Human

Using the S.3 sky regions: binned noise (1 to 5×5) vs. predicted from A, B and C. Report.

👤 Human check: confirm the reported vs. empirical noise curves look right on the sky regions; if CWITools (or another tool you use) reads covariance, run it on one cube and confirm it does not break on the new extensions.

PR6 = M.5a–M.5c.

### M.6 Default flip — Sonnet 5.5 (low)

Depends on: M.5c and the S.4 Q8 answer. `cube_method = "matrix"` default in `kcwi.cfg` and `kcwi_koarti.cfg`; changelog entry stating that invsens must be re-derived with the new method and that `cube_method = "warp"` remains available (D5).

Validation:
- A reduction with the default config of both S.3 channels produces `CUBEMETH = 'matrix'` on every cube and invsens.
- Re-using an old (warp) invsens triggers the M.2d guard.
- `test_parse_configs.py` covers the new keys.

👤 Human check: present M.2e/M.4/M.5c reports to the Keck maintainers for G6; record the decision in §2.

PR7 = M.6.

### M.7 Performance and memory report — Sonnet 5.5 (low)

A table per channel: wall time and peak memory of GenerateMaps, MakeCube, WavelengthCorrections, CorrectDar, warp vs. matrix (with and without covariance). Goes in the PR7 description and the changelog summary.

Validation: numbers come from logged timings of real runs on the S.3 machine, not estimates.

---

## 9. Release

### R.1 Leftover sweep — Sonnet 5.5 (low)

Bokeh-plan 5.2 search list, plus: bare `input(` or `plot_pause` sleeps in primitives (should be none outside the helper); `TODO`/`FIXME` added by this work; config keys mentioned in docs but not in `kcwi.cfg` or vice versa. Classify every hit as intentional or leftover.

Validation: the intentional list matches bokeh-plan 5.2's expected list plus the `warp` code path (D5).

### R.2 Independent review — Opus 5.5 (high), fresh session

Replaces bokeh-plan 5.3. Over the diff from `upgrade-start` to the release candidate, with no access to implementation conversations. Bokeh-plan 5.3 focus list, plus:
- Track P: any path where a prompt can change a reduction differently from the logged answer; EOF and strict-mode behaviour; any prompt whose order changed.
- Track M: half-pixel conventions; `ncand` adequacy; masking before normalization; flag semantics; WCS/CRPIX with padding; variance and covariance formulas; mixing methods; memory at full-frame unbinned size.
- No mode changes in the diff (S.1).

Validation: every finding is fixed with a test or closed with a written reason.

### R.3 Clean-machine acceptance — Human

Bokeh-plan 5.4 in full (Ubuntu 24.04, no browser installed), plus, in the same session:

👤 Human check:
- Default-config reductions of both S.3 channels produce matrix cubes and pass the M.0 thresholds against S.3.
- One `plot_level = 3` interactive session over the remote-viewing route from S.4, replaying `answers_<channel>_3_v2.txt` with `KCWI_STRICT_INPUT=1`, then one live session answering by hand.
- `cube_method = "warp"` still runs end to end.

### R.4 Documentation branch — Sonnet 5.5 (medium)

Bokeh-plan 5.1b, plus: interactive-mode page (P.7), `cube_method` and covariance keys in `configuration.rst`, a short "cube construction" explanation and the invsens re-derivation note, `api.rst` entries for `cube_matrix` and `kcwi_interactive`.

Validation: bokeh-plan 5.1b checks; new pages build without warnings.

👤 Human check: the person identified in S.4 merges it; check the published pages.

### R.5 Release — Human, with Keck maintainers

Version bump, consolidated changelog, release notes listing: Bokeh removal and new viewer; new prompt behaviour; each bug fix with its effect; matrix default and the invsens re-derivation requirement; new keys and products. Gate G7.

---

## 10. Consolidated human checklist

| When | Task | What you do | Done when |
|---|---|---|---|
| Start | S.1 | Decide fate of `stash@{0}` and `geckodriver_prefer_non_snap`; tag `upgrade-start` | Tag exists, tree clean |
| Start | S.2 | Confirm baseline env runs Bokeh and opens the browser | Versions and test baseline recorded |
| Start | S.3 | Run blue and red reference reductions at levels 0, 1, 3; capture answers; archive outputs and screenshots; check replay | `docs/baseline/README.md` written, replay identical |
| Start | S.4 | Ask maintainers questions 1–12 | Answers in §2 |
| After S.5 | S.5 | Spot-check 5 figure rows and 5 prompt rows | Noted in inventory |
| Track A | A.0.3 | Run prototype criteria in a browser on macOS and Ubuntu | Decision written in bokeh plan |
| Track A | A.3.1 | Interactive MakeInvsens session with new cursor readout | Noted in PR1 |
| Track A | A.4.3 | Review comparison page; interactive level-3 session; sign-off | Sign-off with accepted differences |
| Before P.1 | P.0 | Review every prompt explanation for scientific accuracy | Each section marked approved |
| Track P | P.5 | Approve or reject each bug fix; inspect B1 plots | Decision per fix in PR3b |
| Track P | P.6 | Interactive session with invalid input, `?`, Ctrl-D; read the log; record `answers_*_v2.txt` | Notes in PR3a; v2 file replays |
| Before M.2e | M.0 | Choose sky regions and std aperture; confirm thresholds | `regions_*.json` and thresholds committed |
| Track M | M.2e | Blink warp vs. matrix cubes; check coverage; sign off report | Sign-off in PR4 |
| Track M | M.3 | Compare one sky-line profile | Noted in PR5 |
| Track M | M.4 | Check std star does not drift; compare invsens | Noted in PR5 |
| Track M | M.5c | Check noise curves; downstream tool smoke test | Noted in PR6 |
| Track M | M.6 | Present evidence to Keck for G6 | Decision in §2 |
| Release | R.3 | Clean-machine run, remote interactive session, warp fallback | Checklist complete |
| Release | R.4–R.5 | Docs merge, release notes, version | Released |

---

## 11. Risks, known gaps and rollback

| Risk / gap | Mitigation |
|---|---|
| No nod-and-shuffle real data (Q4) | N&S validated on synthetic frames in P.4 and M.2c. 👤 Optional: fetch one N&S dataset from KOA before M.6 and run it; if not done, the release notes say N&S under `matrix` is validated on synthetic data only |
| Matrix default changes science results for all users (Q5) | Gate G6, method guard M.2d, `warp` kept (D5), thresholds fixed before results are seen (M.0) |
| Prompt wording may be a de-facto interface at Keck | S.4 Q6 before P.0 |
| Bug fixes change outputs at default `plot_level` (B1) | Separate PR3b, measured effect per fix, maintainer gate G4 |
| Long-running sequence; upstream `develop` moves | Rebase each track branch on `develop` before its PR; rerun the S.3 replay after each rebase |
| Bokeh plan viewer decision branches (WebAgg or fallback) | Only Track A is affected; P and M are independent of the viewer |
| Memory for unbinned (1×1) frames in `build_W` | M.2a builds per slice; R.2 checks full-frame unbinned memory explicitly |

Rollback: each PR is independently revertable. Until M.6, the matrix code is opt-in. After M.6, `cube_method = "warp"` restores old cubes without a code change. Track A rollback is as in the bokeh plan.

---

## 12. Prompt-site bugs found while planning (input to P.5)

Found by reading the code at `dfbcf9d`. Each must be confirmed by a failing test in P.5 before it is fixed. Line numbers will move after Track A.

| ID | Where | Bug | Effect | Changes outputs? |
|---|---|---|---|---|
| B1 | `MakeMasterFlat.py:475-481` | When more than one peak is found and `plot_level < 2`, `ipk = properties['peak_heights'].argmax()` is an index into `peaks`, but is then used as a pixel index (`apk = xvals[ipk]`). It should be `peaks[argmax]`. | The flat ledge position is taken near the start of the trimmed range instead of at the highest peak. **This runs at the default `plot_level = 1`** | Yes, when multiple peaks occur |
| B1b | Same block | When **no** peak is found and `plot_level < 2`, `argmax()` of an empty array raises `ValueError` and stops the primitive. | Crash | Yes (crash → result) |
| B1c | `MakeMasterFlat.py:469-474` | Interactive peak: `isnumeric()` loop accepts any non-negative integer with no bounds check; out of range raises `IndexError`. No way to accept the automatic choice. | Crash on typo | No for valid input |
| B2 | `MakeInvsens.py:469-480` | Line-mask entry with only one number: `qstr.split()[1]` raises `IndexError`, but only `ValueError` is caught. | Crashes MakeInvsens | No for valid input |
| B3 | `MakeInvsens.py:655-680` | Fit order: `ford = int(qstr)` is assigned before the fit; if `np.polyfit` then fails (e.g. negative order), the error is caught but `ford` keeps the bad value while the fit is the old one. Very large orders are not rejected. | Logged/recorded fit order may not match the fit used | Possibly (logged value) |
| B4 | `ArcOffsets.py:92-99`, `SolveArcs.py:438-445`, `FitCenter.py:339-346` | "Next? <int>": a number ≤ the current bar, negative, or ≥ NBARS never matches the loop index, so the interactive review silently ends. | Confusing; user thinks they will see bar *n* | No |
| B5 | All `'Q' in q.upper()` sites | Any input containing `q` quits. | Accidental quit | No |
| B6 | `MakeInvsens.py:365-380` | An invalid limits entry resets *both* limits to WAVGOOD, discarding limits the user accepted earlier in the loop, without saying so. | Surprising reset | Possibly |
| B7 | `NandshuffSubtractSky.py:185-195, 269-` | Sky scale 0 is accepted (only negatives re-prompt), giving no sky subtraction with no warning. | Silent no-op sky subtraction | Yes, if 0 entered |
| B8 | `reduce_kcwi.py:436-469`, `kcwi_rti.py:412-440` | The channel check uses `'kr' in frame.lower()` on the full path, so a directory name containing `kr` or `kb` triggers a false "wrong channel" prompt; the prompt repeats for every file. | Spurious prompts | No |

B1 is the only one known to affect unattended reductions at default settings. S.3's level-1 logs show whether it occurred on the reference data.
