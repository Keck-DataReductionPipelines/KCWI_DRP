# KCWI DRP upgrade: master implementation plan

Date: 2026-10-22, revised the same day after review (decisions D10–D19). Status: **plan only, nothing implemented.**
Plan documents live on branch `drp_upgrade_plan`. Code facts in this plan (counts, line numbers, bug locations) were verified at `dfbcf9d` (`master` `12fd4aa` plus planning documents only). Upstream has since released `v1.3.2` (`9be28f5`, PR #218), which differs from `12fd4aa` only in `kcwidrp/core/kcwi_plotting.py`, so primitive line numbers and counts still hold at `v1.3.2`. The work starts from `v1.3.2` (D10).

This document coordinates three upgrades, and the test overhaul that protects them, into one sequence of tasks, pull requests and human checkpoints:

| Area | Detailed design lives in | This document adds |
|---|---|---|
| **T.** Test overhaul: unit, invariant, config/CLI and manual regression tests | This document (§6) | Full plan. Lands first, as the safety net for A, P and M |
| **A.** Bokeh → Matplotlib plotting | [bokeh_port_plan.md](bokeh_port_plan.md) (tasks 0.1–5.4), rationale in [bokeh_replacement_proposal.md](bokeh_replacement_proposal.md) | Amendments, model remapping, shared baseline, PR gates |
| **P.** Descriptive, robust interactive prompts | This document (§8) | Full plan |
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

All human steps are also collected into one checklist in §11 so nothing is missed.

**Standing rule for every PR after PR0** (Track T): before opening it, run the fast suite (`pytest`) and, if the PR touches reduction code, the regression suite (`pytest -m regression`, quick profile; full profile where T.7 requires it, D16), and paste both summaries into the PR description. A regression failure is resolved either by a fix or by a *re-bless* (§6, T.7) with a written reason. There is no CI (Q11), so this rule is the only thing that enforces the tests.

---

## 2. Decisions record (questions and answers)

Questions asked while writing this plan, and the answers given. Later changes go in as new rows with a date, not as edits to old rows.

| # | Question | Answer (2026-10-22) | Consequence in this plan |
|---|---|---|---|
| Q1 | How should the three areas be sequenced and delivered? | **Sequential PRs.** Bokeh port first, then prompts, then matrix stages 2+. Matrix stages 0–1 run in parallel from the start. One PR per area or stage. | §4 ordering. PR1 (A), PR2 (M0–M1, parallel), PR3a/3b (P), PR4–PR7 (M) |
| Q2 | What should "more descriptive messages" cover? | **Explanatory text + input robustness + fix logic bugs found.** On-plot guidance in the viewer was *not* selected. | §8. Viewer is unchanged by Track P. Bug fixes are a separate PR (3b) because they change outputs |
| Q3 | Which matrix stages are in scope? | **Stages 0–5.** Stage 6 (joint inversion of dithers) is out of scope. | §9 covers M0–M5. Stage 6 stays in matrix_cube_plan.md as future work |
| Q4 | What reference data is on hand? | **Blue + red with standard stars, no nod-and-shuffle.** | S.3 uses BL/BM + KCRM sets. N&S paths are validated on synthetic data only; this is a known gap (§12) |
| Q5 | After validation, what is the shipped default for `cube_method`? | **Flip to `matrix`** in this update. | M.6 makes the flip, gated by Keck maintainer sign-off (G6). Invsens/FluxCalibrate method guard (M.2d) becomes mandatory |
| Q6 | Which covariance products should Stage 5 produce? | **A** (neighbour correlations), **B** (store `W`, opt-in) and **C** (noise-scaling curve in header). | M.5 implements all three. A and C on by default, B behind a config key |
| Q7 | Who makes maintainer-level decisions? | **You plus Keck maintainers.** You decide day-to-day; upstream merges and user-facing defaults need Keck sign-off. | Gates G1–G7 in §4 name what needs Keck sign-off |
| Q8 | Should prompt answers be recorded in the pipeline log? | **Yes, at INFO.** Explanations go to the terminal only. | P.1 helper logs every answer and its effect |
| Q9 | (Instruction, not a question) Which models may be used? | **Opus or Sonnet only.** | Bokeh-plan tasks assigned to Haiku move to Sonnet; the Fable review moves to Opus in a fresh session (§3) |
| Q10 | (Instruction, 2026-10-22) Overhaul the tests: unit tests plus manually launched regression tests on larger datasets, and other recommended test types. | Added as Track T. | §6 |
| Q11 | Where should the fast automated tests run? | **Locally only.** No CI. | The standing pre-PR rule in §1; T.9 writes the checklist. Travis is removed in T.10 |
| Q12 | How do regression tests find their large datasets? | **A local directory named by an environment variable; no download script.** | `KCWI_REGRESSION_DATA` (T.7). You populate it by hand from S.3; a manifest with checksums describes it |
| Q13 | When does the test overhaul happen? | **First, as a safety net.** Infrastructure, characterization tests and the regression runner land before Track A; track-specific tests are added inside each track. | PR0 before PR1 (§4). M.0's regression script becomes part of T.7 |
| Q14 | Which extra test types besides unit and regression? | **Config & CLI tests, property/invariant tests.** Synthetic end-to-end pipeline runs and performance benchmarks were *not* selected. | T.5 (invariants, with `hypothesis`) and T.6 (config, event tables, entry points). Timings are still recorded in regression reports, but nothing fails on them |
| Q15 | (Instruction, 2026-10-22) Remove all Travis CI configuration and references. GitHub Actions will be used eventually, but for now there are **no automated tests**. | Added as T.10, in PR0. | No CI workflow is added in this upgrade. §6.3 records what keeps the suite ready for GitHub Actions later |
| Q16 | Does anything consume Travis CI (people, branch protection, apps, badges)? | **No consumers, human or machine.** Removal is trivial. | T.10 needs no maintainer question and no GitHub settings check |

Decisions taken by the plan author rather than asked, because they have a conventional answer. Overrule any of them by adding a row above.

| # | Decision | Why |
|---|---|---|
| D1 | The matrix plan is copied into `docs/matrix_cube_plan.md` so it is versioned with the code and visible to maintainers. The copy in `~/drp/kcwi/` is no longer the source of truth. | A plan linked from a PR must be readable by reviewers |
| D2 | The reference reductions for Track A (bokeh 0.1) and Track M (Stage 0) are one human session (S.3), not two. | Same data, same runs; halves the human effort |
| D3 | In Track P, the number and order of prompts **for valid answers** does not change. Re-prompting on *invalid* input is allowed. | Keeps the S.3 answer files usable as a regression harness for every later PR |
| D4 | End-of-input (Ctrl-D / closed stdin) at a prompt: by default, log a warning and take the default answer for that and all later prompts. With `KCWI_STRICT_INPUT=1`, raise `EOFError` as today. | Robustness for users; the strict mode keeps the bokeh-plan replay check ("ends with `EOFError` → prompts changed") working |
| D5 | The `warp` path is kept after the default flip, selectable with `cube_method = "warp"`, for at least one release. | Rollback without a code revert; side-by-side comparison for users |
| D6 | New pull requests target `develop` (upstream merges `develop` → `master`, as in PR #207). | Matches existing upstream practice |
| D7 | Tests stay inside the package at `kcwidrp/tests/`, split into `unit/`, `properties/`, `config_cli/` and `regression/`. pytest configuration moves to `pyproject.toml`. | Keeps today's location; one config file |
| D8 | Existing bugs are pinned by tests marked `xfail(strict=True, reason="<bug ID>")`, not by asserting the buggy output. The fix PR removes the marker. | A characterization test must not make a bug look intended; `strict` makes the test fail loudly once the bug is fixed and the marker is forgotten |
| D9 | `pytest` moves from a runtime dependency (currently pinned `>=5.4.1,<5.5.0`, too old for Python 3.12) to a Poetry `test` dependency group, with `hypothesis` and `pytest-cov`. Travis removal is a separate decision (Q15, T.10). *Amended by D13.* | Users should not install pytest to reduce data. Needs Keck agreement (S.4 Q13) |

Decisions added after review (2026-10-22). Where one changes an earlier row, the earlier row is marked "amended by".

| # | Decision | Why |
|---|---|---|
| D10 | Work starts from upstream `v1.3.2`, not `12fd4aa`. `git fetch`, check that `develop` contains `v1.3.2`, branch from it, and tag `upgrade-start` there. The `geckodriver_prefer_non_snap` question in S.1 is closed: it was merged as PR #218. | `v1.3.2` was released 2026-10-15 and is not an ancestor of `dfbcf9d`; PRs against `develop` must not undo it |
| D11 | Two existing crashes when the viewer is off are fixed in PR0 (T.11, bugs B9 and B10 in §13): `context.bokeh_session` is initialised to `None` by the launch scripts, and `bokeh_clear(None)` returns. Behaviour with a working Bokeh server is unchanged. | Without it, no reduction at `plot_level ≥ 1` runs with the viewer disabled, so the regression suite (T.7) cannot run unattended and headless on `upgrade-start` code |
| D12 | Track P's "no behaviour change" rule is: **for input that is valid today, PR3a must give identical prompts, logged results and outputs; input that is invalid today (crashes, silent resets, silent early exits, accidental quits) may re-prompt instead.** Bugs whose fix only affects invalid input (B1c, B2, B4, B5, B6) are fixed by the helper in PR3a. Bugs whose fix changes outputs, or the prompts seen with valid input, (B1, B1b, B3, B7, B8, and later B11–B13) are fixed in PR3b. Supersedes the "keep the bug in P.4" wording. | §8.2's helper design already fixes B4 and B5 by construction; the earlier text put the same fixes in both PRs |
| D13 | Test dependencies are a pip-installable extra (`pip install -e ".[test]"`), declared as PEP 621 `[project.optional-dependencies]` (supported by `poetry-core>=2`) or `[tool.poetry.extras]`, not a Poetry-only dependency group. Amends D9. | Users and maintainers install with pip and conda; a Poetry group cannot be installed by pip |
| D14 | Two environments, both Python 3.12: the **baseline env** (repaired `kcwidrp` conda env plus the test extra) runs `upgrade-start` code and all of PR0; the **development env** for PR1 onward has the same numeric pins (`numpy`, `scipy`, `astropy`, `scikit-image`, `astroscrappy`) as the baseline, minus Bokeh. The environment check in T.7 is done in the baseline env, before any code change. | `upgrade-start` cannot import without Bokeh, so it cannot be run in a post-PR1 env; `numpy==1.26.4` has no Python 3.13 wheels |
| D15 | The regression runner does not trust the exit code or stdin position. Failures are found by scanning the log; prompts are recorded by a `sitecustomize.py` that the runner places on `PYTHONPATH` and that wraps `builtins.input` to append each prompt and answer to a sidecar file. | The framework catches primitive exceptions, carries on, and always leaves through `os._exit(0)`; Python reads a redirected stdin file in 8 KB blocks, so file position shows nothing about how many answers were used |
| D16 | Two regression profiles. **Quick** (default for `-m regression`): `plot_level = 1`, both channels. **Full** (`--regression-profile=full`): levels 0, 1 and 3. Full is required for PRs that touch plotting or prompt code (PR1, PR3a, PR3b), and for release. | Levels 0 and 1 give the same reduction outputs, since both are below `plot_level` 2; six full reductions for every PR is unnecessary |
| D17 | From M.2c, the regression manifest has explicit `*_warp` datasets pinned with `cube_method = "warp"`. They keep their references through the M.6 default flip. | Without them, blessing the default datasets at M.6 leaves the `warp` rollback path (D5) untested |
| D18 | Tests and the regression runner set `astropy.utils.iers.conf.auto_download = False`; the IERS table versions in use are recorded in the regression provenance. | `radial_velocity_correction` may otherwise download IERS-A, which breaks the "no network" rule and can change `icubew` with no code change |
| D19 | PR6 (covariance) gets its own gate, G5b. S.4 Q10 is split: Q10a (`CUBEMETH`, `WAVCORR`, `DARCOR`) before M.2c; Q10b (noise-curve keywords, named in M.5a) before M.5b. | PR6 adds user-facing keys, extensions and header keywords (Q7), and the M.5 names do not exist when M.2c starts |
| D20 | (Maintainer decision, 2026-10-22.) Upstream `develop` was fast-forwarded to `master` (`9be28f5`, `v1.3.2`), so `upgrade-start` is on `develop` and PRs target `develop` (D6) with no extra commits. `stash@{0}` (an earlier draft of PR #218 plus mode flips) was dropped. S.1 is complete. | `develop` had stayed at `6058c37` while `master` received PR #207's merge and PR #218 |

Maintainer answers to the S.4 questions (gate G0), given 2026-10-22 by the maintainer. Numbered after the S.4 list (S4-Q*n*) so they do not clash with Q1–Q16 above. S4-Q10b is still open (asked after M.5a).

| # | Question (S.4) | Answer | Consequence in this plan |
|---|---|---|---|
| S4-Q1 | Does anything downstream read `plots/*.png`? | **Only file names and PNG format matter.** | Track A may change plot appearance; names, format and pixel size are preserved as the bokeh plan already requires |
| S4-Q2 | How are plots viewed when the pipeline runs remotely? | **VNC / remote desktop.** The browser runs on the pipeline host. | The `localhost`-only viewer binding works as designed. R.3's remote session is done over VNC. SSH port-forwarding stays documented but is not tested |
| S4-Q3 | Several processes watched at once; one tab per process acceptable? | **Yes, one tab per process.** | The default topology from bokeh plan 0.3 is confirmed; the shared-port alternative is dropped |
| S4-Q4 | Auto-open the browser tab on summit and RTI hosts? | **Off everywhere.** | `plot_viewer_open_browser = False` in both `kcwi.cfg` and `kcwi_koarti.cfg` (§7 amendment 13). This changes today's behaviour (`session.show` opens a tab); the changelog says so and the printed URL is the way in |
| S4-Q5 | Who merges the `documentation` branch; how is ReadTheDocs rebuilt? | **User documentation lives in a separate repository, [KCWI_DRP-Docs](https://github.com/Keck-DataReductionPipelines/KCWI_DRP-Docs) (`main`), maintained by the maintainer, and is updated under a separate plan.** | All `documentation`-branch and `docs/source/` work leaves this plan: T.10's second half, bokeh 5.1b and R.4 become a hand-off list for the docs plan (§10 R.4) |
| S4-Q6 | Are prompts driven by scripts at Keck? | **No, humans only.** | Prompt wording can change freely in Track P; P.0 needs no maintainer wording approval |
| S4-Q7 | Log user answers at INFO? | **Yes.** | As planned (Q8) |
| S4-Q8 | Flip the default `cube_method` to `matrix` in this release? | **Yes, if every M.0 threshold passes on both channels and the M.2e, M.4 and M.5c human checks are signed off; otherwise the release ships with `matrix` opt-in.** | M.6's go/no-go is mechanical; G6 checks the evidence against this rule |
| S4-Q9 | New extensions in `icube*` files, or separate files? | **Separate files.** | The coverage cube (M.2c) and covariance products (M.5b) go in new files next to the cube, e.g. `<frame>_icube_coverage.fits`, `<frame>_icube_corr.fits`; the HDU layout of every existing output is unchanged. Exact names are proposed in M.2c/M.5b and listed in the changelog |
| S4-Q10a | Header keywords `CUBEMETH`, `WAVCORR`, `DARCOR`? | **Yes, as proposed.** `CUBEMETH` is written in both modes. | M.2c, M.3, M.4 as written |
| S4-Q11 | Opt-in `W` file: acceptable, and where? | **Yes. Write a routine that saves `W` in COO form as a FITS file with a single binary-table extension, next to the cube.** | M.5b writes `<frame>_icube_W.fits` (option B) with `cube_matrix.write_matrix_fits` / `read_matrix_fits`; no `.npz` (§9 amendment 10) |
| S4-Q12 | Are old invsens files reused across reductions? | **Rarely.** | M.2d's proposed default stands: on a method mismatch, log an ERROR naming both files and the fix, and skip flux calibration for that frame |
| S4-Q13 | Move `pytest` into a `test` extra? | **Yes.** | T.1 as written (D13) |
| S4-Q14 | Flip the RTI config to `matrix` too? | **Yes, same as `kcwi.cfg`.** | M.6 changes both configs, under the same S4-Q8 rule |

---

## 3. Model assignments

| Model | ID | Used for |
|---|---|---|
| Opus 5.5 | `claude-opus-5-5` | Design that many later tasks inherit; numerics (`W`, variance, covariance, WCS); concurrency (plot viewer); domain-heavy text (prompt explanations); first-pass review of real-data results; the final independent review |
| Sonnet 5.5 | `claude-sonnet-5-5` | Well-specified implementation; mechanical call-site conversion; tests; scripts; config, docs and changelog; inventories and sweeps |

Rules carried over from the bokeh plan and adapted to two models:

1. **Tests are written by a different model, or a fresh session, from the code they test**, from the specification and not from the implementation's source. Where Opus implements, Sonnet writes tests. Where Sonnet implements, tests come from Opus or from a separate Sonnet session with no access to the implementation conversation. This includes tasks that both implement and test (T.11, M.2d, M.5b, M.6): their tests are written in a separate fresh session.
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
Track T: tests (T.0–T.6, T.10, T.11 can start at once; T.7 needs S.3 and T.11)
  T.0 strategy ─► T.1 infra ─► T.2 primitive harness ─► T.3–T.6 unit/invariant/config
                  T.1 ─► T.11 viewer-off fixes ─┐
                                           S.3 ─┴─► T.7 regression suite ─► T.8 ─► PR0 ─► (G0T)
        │
        ├──────────────► Track A: Bokeh port ─────────► PR1 ──► (G1)
        │                                                  │
        └─► (after T.7) M.0 ─► M.1a/b maps & pixcoords ─► PR2 ─┤ (G2)  [parallel with A]
                                                           │
                         Track P: prompts (after PR1 merged)
                           P.0–P.4, P.6, P.7 ─────────► PR3a ─► (G3)
                           P.5 bug fixes, P.6, P.7 ───► PR3b ─► (G4)
                                                           │
                         Track M: matrix (after PR3a merged)
                           M.2 W in MakeCube (opt-in) ─► PR4 ─► (G5)
                           M.3 + M.4 λ corr + DAR ─────► PR5
                           M.5 covariance ─────────────► PR6 ─► (G5b)
                           M.6 default flip + M.7 ─────► PR7 ─► (G6)
                                                           │
                         Release: R.1 sweep, R.2 review, R.3 clean machine,
                                  R.4 docs hand-off, R.5 release ─► (G7)
```

Why this order:

- **T before everything (Q13).** Every later PR is checked by the characterization tests and the regression suite. Without them, "outputs unchanged" in Tracks A and P can only be checked by hand. Work on T.0–T.6 needs no data and can start on day one, alongside S.1–S.4. Only T.7 waits for S.3.
- **Track A may start before PR0 merges**, once T.1–T.3 are committed on the PR0 branch, because the bokeh plan has its own wrapper tests. PR1 must not merge before PR0. (The task table's dependency for A.* says the same.)
- **A before P.** Track P rewrites the code around `bokeh_plot(...)` calls at most prompt sites, and the MakeInvsens limits prompt depends on the new viewer's cursor readout replacing Bokeh tooltips. Doing P first would conflict with every file in A.
- **A before M.2.** `MakeCube` is touched by A (image plot, possibly the `Pool` start method) and rewritten by M.2. Bokeh-plan task 4.2 also requires cubes to be numerically identical to the reference, which only holds while the warp path is the only path.
- **P before M.2** (per Q1). The file overlap is small (one prompt in `MakeCube`), so if Track P stalls, M.2 may start after PR1 instead; record that as a new decision row if it happens.
- **M.0–M.1 in parallel.** They touch `GenerateMaps`, a new `core/cube_matrix.py` and new tests, none of which Track A changes. M.0 needs T.7's regression tooling, so it starts once T.7 is committed on the PR0 branch.

Branches. Each PR has its own branch, named `upgrade/pr<N>-<topic>` (for example `upgrade/pr0-tests`, `upgrade/pr1-bokeh`). PR0's branch starts at `upgrade-start` (D10). Each later branch starts from `develop` after the PR it depends on has merged upstream. If that PR has not merged yet, the branch starts from that PR's branch and is rebased onto `develop` once it merges. PR2 starts from the PR0 branch. Planning documents stay on `drp_upgrade_plan` and are merged with PR0. The task table below is the authority on dependencies; the diagram above is a summary.

Gates. Each needs the named sign-off before the next step starts.

| Gate | When | Who signs | What they agree to |
|---|---|---|---|
| G0 | End of S.4 | You + Keck maintainers | Answers to the maintainer questions in S.4 |
| G0T | PR0 review | Keck maintainers | New test layout, pytest moved to a `test` extra (D13), all Travis configuration and references removed (T.10), the viewer-off fixes B9 and B10 (T.11), the local pre-PR testing rule replacing "a pass from the automatic TravisCI tests" as a PR requirement |
| G1 | PR1 review | Keck maintainers | Bokeh removal, new viewer, config key changes, one viewer tab per process |
| G2 | PR2 review | You (Keck review as normal) | Vectorized GenerateMaps; no behaviour change |
| G3 | PR3a review | Keck maintainers | New prompt wording and behaviour, answers in log, EOF handling, the invalid-input fixes listed in D12 |
| G4 | PR3b review | Keck maintainers | Each bug fix and its measured effect on outputs |
| G5 | PR4 review | Keck maintainers | New `cube_method` key, new header keywords and extensions (opt-in) |
| G5b | PR6 review | Keck maintainers | Covariance config keys (`cube_covariance`, `cube_save_matrix`), the option A extension and option C header keywords, the opt-in `W` file (D19) |
| G6 | PR7 review | Keck maintainers | Flipping the default to `matrix`; users must re-derive invsens |
| G7 | Release | You + Keck maintainers | Version number, release notes |

### Task overview

| ID | Task | Model (effort) | Depends on | 👤 |
|---|---|---|---|---|
| S.1 | Repo hygiene and starting tag | Human | – | ✔ |
| S.2 | Environments (repair baseline env; D14) | Human + Sonnet (low) | S.1 | ✔ |
| S.3 | Reference reductions (shared baseline) | Human | S.1, S.2 | ✔ |
| S.4 | Maintainer questions | Human | – | ✔ |
| S.5 | Plot and prompt inventory | Sonnet (low) | – | spot-check |
| T.0 | Test strategy and conventions | Opus (medium) | – | ✔ |
| T.1 | Test infrastructure and packaging | Sonnet (medium) | T.0 | – |
| T.2 | Primitive test harness | Opus (medium) | T.1 | – |
| T.3 | Characterization unit tests: core and helper functions | Sonnet (medium) | T.1 | – |
| T.4 | Primitive-level unit tests | Sonnet (medium); Opus (medium) for the numeric ones | T.2 | – |
| T.5 | Property and invariant tests | Opus (medium) | T.3 | – |
| T.6 | Config, event-table and CLI tests | Sonnet (medium) | T.1 | – |
| T.7 | Manual regression suite | Opus (medium) design + Sonnet (medium) build | T.1, T.11, S.3 | ✔ |
| T.8 | Suite quality check (mutations, coverage) | Sonnet (medium), fresh session | T.3–T.7 | ✔ |
| T.9 | Testing guide and pre-PR checklist | Sonnet (low) | T.8 | – |
| T.10 | Remove Travis CI | Sonnet (low) | – | – |
| T.11 | Fix viewer-off crashes B9, B10 | Sonnet (low); tests in a fresh session | T.1 | – |
| A.* | Bokeh port tasks 0.3, 1.1–5.1 (not 5.1b) | per bokeh plan | S.3, S.5, T.1–T.3 committed (PR1 merges after PR0) | ✔ |
| M.0 | Matrix fixture, metrics thresholds and regions | Sonnet (medium) | T.7 | ✔ |
| M.1a | Vectorize GenerateMaps | Sonnet (medium) | M.0 | – |
| M.1b | Pixel-edge coordinate helper | Opus (medium) | M.0 | – |
| M.1t | Tests for M.1a/M.1b | Sonnet (medium), fresh session | M.1a, M.1b | – |
| P.0 | Prompt catalogue and explanation text | Opus (medium) | S.5 | ✔ |
| P.1 | Interactive helper module | Opus (medium) | P.0 | – |
| P.2 | Helper tests | Sonnet (medium) | P.1 | – |
| P.3 | Convert simple pause/next sites | Sonnet (medium) | P.1, P.2, PR1 merged | – |
| P.4 | Convert value-entry and bar-jump sites | Opus (medium) | P.3 | – |
| P.5 | Output-changing bug fixes (one commit each) | Opus (medium) + Sonnet tests | P.4, PR3a merged | ✔ |
| P.6 | Interactive acceptance session | Human | P.4 (and P.5 for PR3b) | ✔ |
| P.7 | Changelog and user docs for prompts (PR3a part, then PR3b part) | Sonnet (low) | P.6 | – |
| M.2a | Build `W` | Opus (high) | M.1b, PR3a merged | – |
| M.2b | Tests for `W` | Sonnet (medium) | M.2a | – |
| M.2c | `cube_method` in MakeCube | Opus (medium) | M.2a, M.2b, S.4 Q8, Q9, Q10a | – |
| M.2d | Cube-method guard on invsens / FluxCalibrate | Sonnet (medium); tests in a fresh session | M.2c, S.4 Q12 | – |
| M.2e | Real-data validation of M.2 | Opus (high), then Human | M.2c, M.2d | ✔ |
| M.3 | Fold wavelength corrections into `W` | Opus (medium) | M.2e | ✔ |
| M.4 | Fold DAR into `W` | Opus (high) | M.3 | ✔ |
| M.5a | Covariance mathematics | Opus (high) | M.4 | – |
| M.5b | Covariance I/O and config; tests in a fresh session | Sonnet (medium) | M.5a, S.4 Q10b | – |
| M.5c | Covariance validation | Opus (high), then Human | M.5b | ✔ |
| M.6 | Default flip to `matrix` | Sonnet (low); tests in a fresh session | M.5c, S.4 Q8 and Q14 | ✔ |
| M.7 | Performance and memory report | Sonnet (low) | M.6 | – |
| R.1 | Leftover sweep (code and tests) | Sonnet (low) | all | – |
| R.2 | Independent review | Opus (high), fresh session | R.1 | – |
| R.3 | Clean-machine acceptance | Human | R.2 | ✔ |
| R.4 | Documentation hand-off list (docs repo is out of scope, S4-Q5) | Sonnet (medium) | R.3 | ✔ |
| R.5 | Release | Human + Keck | R.4 | ✔ |

---

## 5. Phase S: shared setup

### S.1 Repo hygiene — Human

Supersedes bokeh-plan task 0.5, which is mostly done: the working tree on `replace_bokeh` is clean, and the earlier mode flips (644 → 755) plus an uncommitted geckodriver fallback in `kcwi_plotting.py` were stashed as `stash@{0}` ("On master: pre-replace_bokeh: …").

👤 Human check:
1. Decide what happens to `stash@{0}`. The mode flips should be dropped (`git config core.filemode false` if they keep coming back). The geckodriver fallback was released upstream as PR #218 / `v1.3.2`, so the stashed copy can be dropped.
2. `git fetch origin --tags` (the local `origin/*` refs date from March 2026). Confirm that `origin/develop` contains `v1.3.2`; if it does not, ask the maintainers which branch new PRs should start from, and record the answer as a decision row (D6, D10).
3. Tag the starting point on that commit: `git tag upgrade-start v1.3.2` (or the `develop` commit containing it). Merge `drp_upgrade_plan` into the PR0 branch so the plans travel with PR0.
4. Note that `pyproject.toml` at `v1.3.2` still says `version = "1.3.1"`; R.5 sets the version.

Validation:
- `git merge-base --is-ancestor v1.3.2 upgrade-start` succeeds.
- `git status --short` is empty on the PR0 branch.
- `git diff --summary upgrade-start` shows no `mode change` lines; this is re-checked in R.2.

### S.2 Environments — Human, with Sonnet (low) for the scripted parts

The shell environment active at time of writing (`pypeit`, Python 3.13) does **not** have `bokeh` installed, so it cannot run the current pipeline. The baseline must run in the `kcwidrp` conda environment (Python 3.12).

**That environment is currently broken** (checked 2026-10-22): `import psutil` fails with `AttributeError: module 'psutil._psutil_osx' has no attribute 'getpagesize'`, because `site-packages/psutil/` holds both the 5.7.3 extension and an `abi3` extension left by a newer psutil. As a result `reduce_kcwi --help` fails. It must be repaired before anything else in S.2.

Do (environments as in D14):
1. Repair the baseline environment: `conda activate kcwidrp`, then `pip uninstall -y psutil` (repeat until it reports nothing installed), delete any leftover `site-packages/psutil/` directory, and `pip install "psutil>=5.7.0,<5.8.0"`. If the pinned psutil will not build or import on this Mac, rebuild the env from scratch with Python 3.12 and record the reason; the `psutil` pin is then a question for the maintainers.
2. Baseline environment: `pip install -e .` from `upgrade-start`. Record `pip freeze > docs/baseline/env_baseline.txt`. Add a `.gitignore` rule so that only text files (`*.md`, `*.txt`, `*.json`) under `docs/baseline/` can be committed; the FITS outputs and plots from S.3 stay local. After T.1 lands, `pip install -e ".[test]"` into the same env; PR0 is developed and tested here.
3. Development environment for PR1 onward: a fresh Python 3.12 env from the branch's `pyproject.toml`, with the numeric packages (`numpy`, `scipy`, `astropy`, `scikit-image`, `astroscrappy`) pinned to the exact versions in `env_baseline.txt`. Record its `pip freeze` as `docs/baseline/env_dev.txt`, and diff the two files in the PR1 description. Python 3.13 is not usable while `numpy==1.26.4` is pinned.

Validation:
- `python -c "import psutil, kcwidrp, bokeh, matplotlib; print(bokeh.__version__, matplotlib.__version__)"` succeeds in the baseline env and the versions are recorded.
- `reduce_kcwi --help` exits 0.
- `pytest kcwidrp/tests` result on unmodified code is recorded (pass/fail per test), so later failures can be told apart from pre-existing ones.

👤 Human check: confirm the baseline env can open the Bokeh browser tab on your Mac (needed for S.3 screenshots). If it cannot, S.3 must run on another machine; record which.

### S.3 Reference reductions — Human

Replaces bokeh-plan task 0.1 **and** matrix Stage 0's data capture (decision D2). Follow bokeh-plan 0.1 steps 1–4 exactly, with these additions:

- **Datasets:** one blue night (BL or BM) and one red (KCRM) night, each with bias, flat, arc, continuum bars, a standard star and at least one science frame with a point source. No nod-and-shuffle (Q4). Keep the raw data in directories whose full path contains neither `kr` nor `kb` other than in the file names, so the channel check (bug B8) does not insert extra prompts that would consume answer lines.
- **Viewer:** S.3 runs with `enable_bokeh = True` and a working Bokeh server, as users run today. Runs with the viewer off are not possible on this code (B9, B10).
- **Passes per channel:** `plot_level = 0`, the usual `plot_level` (1), and `plot_level = 3` with answers captured via `tee answers_<channel>_3.txt`. Give at least one non-default *valid* answer at every value prompt (ReadAtlas offset, MakeInvsens limits / one line mask / fit order, MakeMasterFlat peak if it appears). Do **not** type invalid input in this pass: some invalid inputs crash the current code (§13), which would make the replay unusable.
- **Archive per channel**, under `docs/baseline/<channel>/` (data excluded by the S.2 `.gitignore` rule): `plots_ref_<channel>_<level>/`; `answers_<channel>_3.txt`; the full `redux/` directory including `*_geom.pkl`, `*_wavemap`, `*_slicemap`, `*_posmap`, `icube`, `icubew`, `icubed`, `icubes`, `invsens`, the `.proc` table; and the pipeline log for each pass.
- **Timings:** per-primitive wall time from the log for `GenerateMaps`, `MakeCube`, `WavelengthCorrections`, `CorrectDar` (matrix plan §3 baseline). Machine, CPU and core count recorded.
- **Screenshots** of the live-only plots during the level-3 pass, as bokeh 0.1 requires.

Validation (in addition to bokeh-plan 0.1's list):
- The level-3 replay (`reduce_kcwi … < answers_<channel>_3.txt` from a clean output dir) completes unattended and produces identical cubes (`numpy.array_equal` on every extension) and the same PNG names.
- Every output file listed above exists for both channels.
- No log contains `Failed executing primitive` or `while invoking`. The pipeline exits with status 0 even when a primitive fails (the framework logs the exception, carries on with the next event and leaves through `os._exit(0)`), so the log is the only evidence of a failure.
- Search both level-1 logs for `Picking the highest peak`. Record whether it occurs: it decides how bug B1 (§13) can be validated.

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
10. (Split by D19.)
    - 10a. Are new header keywords `CUBEMETH` (warp|matrix), `WAVCORR` and `DARCOR` acceptable? `CUBEMETH` is also written on warp cubes.
    - 10b. Are the noise-scaling-curve keywords from M.5 acceptable? Asked once M.5a has proposed the names.
11. Is a ~160 MB `W` file per frame acceptable as an opt-in product, and where should it be written?
12. Are invsens files from previous DRP versions reused across reductions in practice (e.g. a library of standards)? This decides how strict the M.2d guard should be.

New for Track T:
13. Is it acceptable to move `pytest` out of the runtime dependencies into a `test` extra (`pip install -e ".[test]"`, D13), and add `hypothesis` and `pytest-cov` there?

Added after review (Track M):
14. Should the RTI configuration (`kcwi_koarti.cfg`) also default to `cube_method = "matrix"` at M.6, or stay on `warp` for quick-look? (Ask together with Q8.)

Status: answered 2026-10-22 (§2, S4-Q1–Q14), except Q10b, which is asked after M.5a.

Validation: every question has an answer row in §2 before the task that depends on it starts (question 13 before T.1; 1–5 before A.3.1; 6–7 before P.0; 8, 9, 10a, 11 and 12 before M.2c; 10b before M.5b; 14 before M.6).

👤 Human check: this whole task is yours.

### S.5 Plot and prompt inventory — Sonnet 5.5 (low)

Bokeh-plan task 0.2 as written (file `docs/plot_inventory.md`, 48 figure rows, feature summary), with the prompt table extended for Track P. For every `input(` call, one row with:

- file and line; the `plot_level` guard; the exact prompt string;
- which figure is on screen when it is asked (cross-reference to the figure row);
- accepted answers and what each does (next, quit, jump to bar *n*, new value);
- the default (what `<cr>` does);
- **whether the answer changes the reduction**, and which output it changes;
- the matching `time.sleep(plot_pause)` branch, if any;
- any input that raises or misbehaves (to be checked against §13).

Validation:
- 52 `input(` lines across 17 primitives, plus 8 in `scripts/reduce_kcwi.py` and `scripts/kcwi_rti.py` (counts from `grep` on `dfbcf9d`). `MakeMasterFlat.py:471` and `:473` are one logical prompt and may share a row.
- 36 `time.sleep(...plot_pause)` sites listed.
- The bokeh-plan 0.2 checks (48 rows, 49 `bokeh_plot(`, 20 `save_plot(`, the easy-to-miss features).

👤 Human check: spot-check five figure rows and five prompt rows against the source, including the "changes the reduction" column for ReadAtlas, MakeInvsens and MakeMasterFlat.

---

## 6. Track T: test overhaul

### 6.1 Where the tests are today

At `dfbcf9d` the suite has 249 lines in 5 files and tests very little:

| File | What it does |
|---|---|
| `test_import_primitives.py` (188 lines) | One `import` per primitive, generated by `make_import_test.csh` |
| `test_import_pipeline.py` | Imports the pipeline module |
| `test_parse_configs.py` | Checks that three config files parse; nothing about their contents |
| `test_bokeh_plotting.py` | Checks that a process called `python` is running |
| `make_import_test.csh` | Generator for the import test |

Other problems:
- `kcwidrp/pytest.ini` passes `--noconftest`, which disables `conftest.py`, so shared fixtures cannot be used.
- `pytest` is a **runtime** dependency in `pyproject.toml`, pinned to `>=5.4.1,<5.5.0`, which predates Python 3.12 (the package minimum).
- `.travis.yml` targets Python 3.7 and a CI service that is no longer running.
- No numerical behaviour, header handling, file I/O, config consistency or event-table wiring is tested, and nothing runs a primitive.

Facts about the framework that every test and replay check must allow for (checked against `keckdrpframework` 1.0.0):
- `BasePrimitive.apply` logs a primitive's exception and re-raises it; `Framework.execute` catches it, logs `Exception … while invoking <action>`, and the default `on_state` hook resets the state, so **the reduction carries on with the next event**. Only the failed frame's chain stops.
- The action loop ends with `on_exit(0)`, i.e. `os._exit(0)`: **the exit status is always 0** once the loop has started.
- `Framework.testing = True` makes `execute` re-raise instead. This is useful for in-process tests, but `reduce_kcwi` does not set it.
- With `enable_bokeh = False`, `context.bokeh_session` is never set, so every `bokeh_plot(p, self.context.bokeh_session)` reached at `plot_level ≥ 1` raises `AttributeError` (bug B9); `bokeh_clear(None)` in `SubtractOverscan` raises at `plot_level ≥ 3` (bug B10). T.11 fixes both.

### 6.2 Goals

1. **Fast suite** (`pytest`): unit, invariant, config and CLI tests. Under 60 s on a laptop, with no display, no network and no data outside the repo. Run before every PR (§1).
2. **Regression suite** (`pytest -m regression`): full reductions of real datasets from a local directory (`KCWI_REGRESSION_DATA`, Q12), compared with blessed reference outputs. Launched by hand; skipped cleanly when the variable is unset. Quick and full profiles (D16).
3. **Characterization first.** Tests written in PR0 pin the behaviour of the code at `upgrade-start`. Known bugs are pinned as strict `xfail`s (D8), so later fixes show up as test changes.
4. **Every later track adds to this structure** instead of inventing its own (§6.4).

Layout (D7):

```
kcwidrp/tests/
  conftest.py            shared fixtures (config, synthetic CCDData, tmp redux dir)
  harness.py             run_primitive(): run one primitive without the framework queue (T.2)
  unit/                  core/, primitives/, plotting/ (Track A), interactive/ (Track P), cube_matrix/ (Track M)
  properties/            hypothesis-based invariant tests (T.5)
  config_cli/            config, event tables, entry points (T.6)
  regression/            marker "regression", runner, comparison, report, bless (T.7)
```

### T.0 Test strategy and conventions — Opus 5.5 (medium)

Draft `docs/testing.md` covering: the layout above; markers (`regression`, plus `slow` for fast-suite tests over 2 s, which should be rare); what qualifies as a unit test (writes only to `tmp_path`, no network, no display, under 2 s); the xfail policy (D8); the independence rule (§3 rule 1); the tolerance policy (exact equality for integers, flags and masks; floating-point tolerances only with a written reason next to them); where an expected value comes from (an independent source such as astropy, an analytic result or a published value, or a "characterization" value captured from `upgrade-start` and labelled as such); and the bless policy (T.7). Opus because every later test follows these rules.

Validation: every rule in the document is concrete enough to be checked by reading a test (no "should be reasonable").

👤 Human check: read and approve `docs/testing.md`. In particular, approve the tolerance policy and the rule that floating-point regression comparisons start out exact (T.7).

### T.1 Test infrastructure and packaging — Sonnet 5.5 (medium)

Depends on: T.0, S.4 Q13.

Do:
- Move the pytest configuration into `pyproject.toml` (`[tool.pytest.ini_options]`): `testpaths`, the registered markers, `--strict-markers`, and `-m "not regression"` as the default. Delete `kcwidrp/pytest.ini` and its `--noconftest`. Narrow the blanket `ignore::DeprecationWarning` to third-party modules, so deprecations in our own code are visible.
- A `test` extra (D13): `pytest>=8`, `hypothesis`, `pytest-cov`, installed with `pip install -e ".[test]"`. Remove `pytest` from the runtime dependencies. Regenerate `requirements.txt` with the tool that produced it (`pip-compile`; its header says it last ran under Python 3.11, so run it under 3.12 and note that in the PR).
- Replace the 188-line generated import test and `make_import_test.csh` with one parametrized test that finds every module under `kcwidrp/primitives`, `core`, `pipelines` and `scripts` with `pkgutil.walk_packages` and imports it. `kcwidrp.scripts.smart_reduce_kcwi` fails to import today (bokeh plan 3.2), so its case is a strict `xfail` with reason "A.3.2", as in T.6.
- `conftest.py` fixtures: `kcwi_config(**overrides)` (parsed `kcwi.cfg`, `plot_level = 0`, `plot_pause = 0` by default); `synthetic_ccddata(camera, binning, ampmode, ...)`, which builds a `KCCDData` with the header keywords primitives read (the list compiled by AST scan of `header['…']` reads); `tmp_redux` (a temporary output directory with an empty proctab).
- `conftest.py` sets `astropy.utils.iers.conf.auto_download = False` for the whole session (D18).
- Record baseline coverage (`pytest --cov=kcwidrp --cov-report=term`) per module in the PR description.

Validation:
- `pytest` passes on `upgrade-start` code in under 60 s with `DISPLAY` unset and network disabled.
- Without `KCWI_REGRESSION_DATA`, `pytest -m regression -rs` reports every regression test as skipped, with the reason "KCWI_REGRESSION_DATA not set".
- Adding a dummy module to `kcwidrp/primitives` makes the import test pick it up without editing the test; a module with a syntax error makes it fail.
- `pip install .` in a fresh environment does not install `pytest`; `pip install ".[test]"` does, with `hypothesis` and `pytest-cov`.

### T.2 Primitive test harness — Opus 5.5 (medium)

Depends on: T.1. `kcwidrp/tests/harness.py`, exposed as the fixture `run_primitive(PrimitiveClass, ccddata, config_overrides=None, context_extra=None, answers=None)`. It builds a stand-in `action` (with `args` carrying `ccddata` and the other attributes primitives read) and `context` (`config.instrument`, `logger`, `pipeline_logger`, a `proctab` in `tmp_redux`, `bokeh_session = None`, plus any extras), calls `apply()` and returns the output. `keckdrpframework`'s `BasePrimitive` reads only `action.args`, `context.logger` and `context.config`, so no framework queue is needed.

Defaults: `plot_level = 0`, `plot_pause = 0`. `input()` is patched to fail the test unless `answers` is given, in which case it feeds them in order and fails if any are left unused. `time.sleep` is patched to a no-op. Opus because every primitive test sits on these stand-ins, and stand-ins that differ from the real framework objects give false confidence.

Validation:
- An AST scan lists every `self.context.<attr>` and `self.action.args.<attr>` read anywhere in `kcwidrp/primitives`; the harness documents which ones it provides. A primitive that reads an attribute the harness does not provide fails with a message naming the attribute, not with an `AttributeError` deep inside the primitive.
- `SubtractOverscan`, `TrimOverscan` and `CorrectGain` run under the harness on synthetic data.
- With `plot_level = 0`, a run calls no `input()`, no `time.sleep()` and no plotting function (checked by patching each).
- `answers` with one answer too many, or one too few, fails the test with a clear message.

### T.3 Characterization unit tests: core and helper functions — Sonnet 5.5 (medium)

Depends on: T.1. Unit tests for the stand-alone functions, each expected value labelled with its source (T.0):

| Module | Functions | Expected values from |
|---|---|---|
| `core/geometric.py` | `AsymmetricPolynomialTransform`, `estimate_transform` | Points generated by a known transform; fit recovers it; inverse round-trip |
| `primitives/kcwi_file_primitives.py` | `parse_imsec`, `strip_fname`, `get_master_name`, `get_unique_CCD_master_name`, `get_unique_STATEID_master_name`, `fix_header`, `plotlabel`, `kcwi_fits_writer` ↔ `kcwi_fits_reader` (data, uncertainty, mask, flags, header, table), `write_table` ↔ `read_table` | Hand-built cases; round-trips |
| `core/kcwi_proctab.py` | `Proctab` create, add, `search_proctab` (including `nearest`), write and read back | Hand-built tables |
| `core/kcwi_get_std.py` | `kcwi_get_std`, `is_file_kcwi_std` | Names present in `data/stds`; case variants; non-standards |
| `core/kcwi_correct_extin.py` | `kcwi_correct_extin` | Ratio at two airmasses computed independently from `data/extin` |
| `primitives/CorrectDar.py` | `atm_disper` | Published values (e.g. Filippenko 1982) for stated conditions |
| `primitives/WavelengthCorrections.py` | `a2v_conversion`; heliocentric and barycentric velocity | Standard air→vacuum values for known lines (to 0.01 Å); `astropy` `SkyCoord.radial_velocity_correction` for the same header (to 1 m/s) |
| `primitives/FitCenter.py` | `pascal_shift` | `numpy.polynomial` shift of the same coefficients |
| `primitives/GetAtlasLines.py` | `gaus`, `get_line_window`, `findpeaks` | Synthetic spectra with lines at known positions |
| `primitives/MakeMasterFlat.py` | `bm_ledge_position` | Characterization table |
| `primitives/StackFlats.py` | `inv_median` | Hand-computed |
| `core/bspline/` | fit and evaluate | Smooth analytic function |
| `core/kcwi_plotting.py` | `get_plot_lims` | Hand-computed (survives Track A) |
| `scripts/check_cals.py` | `check_cal_type` | Hand-built proctab |

Validation:
- All pass on `upgrade-start` code, apart from strict `xfail`s that each carry a bug ID.
- Every expected value has a source comment (independent or characterization).
- Each test file runs in under 3 s.

### T.4 Primitive-level unit tests — Sonnet 5.5 (medium); Opus 5.5 (medium) for the numeric and prompt cases

Depends on: T.2. Synthetic inputs with known answers, run through `run_primitive`. Two fixture gaps found in T.2 come first: `synthetic_ccddata()` defaults to 8 overscan columns, below `minoscanpix` (75 for BLUE), so overscan tests must pass a wider overscan or the default must change; and RED frames must number amplifiers 0–3 (`ASEC0`…, `GAIN0`…) as `ingest_file.map_ccd` and `fix_header` expect, which the fixture does not do yet.

Detector level (Sonnet):
- `SubtractOverscan`: constant and sloped overscan per amplifier removed; each amplifier mode in the config. Also at `plot_level = 3` with `bokeh_session = None` (B10, from T.11).
- `TrimOverscan`: output shape and position from `DATASEC`, for 1×1 and 2×2 binning.
- `CorrectGain`: per-amplifier gains from the header applied.
- `FlagSaturation`, `CorrectDefects`: flags set exactly at the expected pixels; defect file chosen by binning and amp mode.
- `CreateUncertaintyImage`: matches the formula from counts, gain and read noise.
- `SubtractBias`, `SubtractDark`: master found through the proctab in `tmp_redux` and subtracted; missing master handled as today.
- `RemoveCosmicRays`: an injected cosmic ray on a smooth image is flagged; smooth pixels are untouched.
- `MakeMasterBias`, `StackFlats`: combination of N synthetic frames; output uncertainty.

Numeric and prompt cases (Opus):
- `CorrectDar` on a synthetic cube with a point source: the shift applied matches `atm_disper` (also the warp-path reference for M.4).
- `WavelengthCorrections` on a synthetic cube with a line at a known air wavelength: the line ends up at the vacuum, velocity-corrected wavelength.
- `FluxCalibrate` with a synthetic invsens: output equals input × invsens / exposure time.
- `CalcPrelimDisp`: the formula, for each grating.
- **Known prompt-site bugs as strict `xfail`s** (D8), driven through `answers` or `plot_level`: B1 and B1b (`MakeMasterFlat` peak selection with several and with zero peaks at `plot_level = 1`), B1c (out-of-range peak number), B2 (one-number line mask), B3 (fit order after a failed fit), B4 (bar number lower than the current one), B5 (an answer containing `q` quits), B6 (invalid limits reset accepted limits), B7 (sky scale 0, `nan`, `inf`), B12 (non-standard N&S branch, synthetic frame), B13 (`.lmsk` round trip). The markers are removed by the PR that fixes each bug, as listed in §13's "Fixed in" column (D12): B1c, B2, B4, B5 and B6 in PR3a (P.4), B1, B1b, B3, B7, B12 and B13 in PR3b (P.5).

Validation:
- Each test builds its inputs in memory or in `tmp_path`, and runs in under 2 s.
- All pass on `upgrade-start` code except the bug `xfail`s, and each of those fails for the stated reason (checked by reading the `xfail` output with `-rx`).

### T.5 Property and invariant tests — Opus 5.5 (medium)

Depends on: T.3. `hypothesis` tests in `kcwidrp/tests/properties/`. Opus because choosing invariants that really hold (and do not merely restate the code) needs care.

Invariants:
- **Transforms:** `inverse(forward(x)) ≈ x` over the detector domain for random well-conditioned coefficients.
- **Air → vacuum:** monotonic increasing; vacuum > air over 3000–11000 Å; the difference is within the physical range.
- **DAR:** `atm_disper(w, w, …) = 0`; swapping the two wavelengths flips the sign; the magnitude grows with airmass and with wavelength separation.
- **Detector arithmetic:** adding a constant to both data and overscan leaves the overscan-subtracted image unchanged; gain correction is linear; the uncertainty image is non-negative and finite wherever the data are.
- **I/O:** `kcwi_fits_writer` → `kcwi_fits_reader` round-trips random arrays (with NaN and masked pixels) and random valid header cards.
- **Sections:** `parse_imsec` round-trips randomly generated section strings.
- **Proctab:** a frame that was added is found by `search_proctab`.

Settings: a default `hypothesis` profile that is derandomized and small enough to keep the fast suite under 60 s, and a `thorough` profile (`--hypothesis-profile=thorough`) for manual runs. Later tracks add their own invariants here (§6.4).

Validation:
- For each invariant, one deliberate breakage of the code makes it fail, and `hypothesis` reports a minimal counterexample.
- The fast suite stays under 60 s.

### T.6 Config, event-table and CLI tests — Sonnet 5.5 (medium)

Depends on: T.1.

- **Config schema:** every key in `kcwi.cfg` and `kcwi_koarti.cfg` parses, and has the type and allowed values given in a schema table in the test (e.g. `plot_level` ∈ {0, 1, 2, 3}, `radial_velocity_correction` ∈ {"heliocentric", "barycentric", "none"}).
- **Config references:** an AST scan finds every `config.instrument.<key>` used in the code (59 distinct keys at `dfbcf9d`), and every `getattr(<…>.config.instrument, '<key>', …)` / `hasattr(…)` with a string literal (used in `StartBokeh.py` and `kcwi_rti.py`). A `getattr` with a non-literal name fails the test with file and line. Each attribute-style key must exist in `kcwi.cfg`; `getattr` keys with a default may be absent but are listed. Keys read only for backward compatibility (from bokeh 5.1: `enable_bokeh`) are on a commented "deprecated" allow-list. A key in `kcwi.cfg` that no code reads fails the test, unless it is on a commented allow-list. Every key in `kcwi_koarti.cfg` exists in `kcwi.cfg`, and the value differences are listed in an allow-list.
- **Event tables** (`kcwi_pipeline.py`, `keck_rti_pipeline.py`): every event's primitive resolves to a class with `_perform` (in `kcwidrp.primitives` or the framework); every "next event" exists in the table; there are no duplicate keys in the dict literal (checked by AST, since Python silently keeps the last one); unreachable events are reported.
- **Entry points:** every target in `[tool.poetry.scripts]` imports, and `--help` exits 0. `smart_reduce_kcwi` currently fails to import (bokeh plan 3.2), so it is a strict `xfail` with reason "A.3.2".
- **`reduce_kcwi` argument handling:** `-f`, `-l`, `-b`, `-r` parsed as documented; the channel check bug B8 and the ignored RTI quit B11 as strict `xfail`s (removed in PR3b, D12).
- **Viewer off** (with T.11): the launch scripts leave `context.bokeh_session` defined (as `None`) when `enable_bokeh = False`.

Validation: each of these breakages makes a test fail: a misspelled key reference (`config.instrument.plot_levl`); a duplicated event key; an event pointing to a primitive that does not exist; a "next event" that does not exist; `plot_level = "1"` (string) in a config.

### T.7 Manual regression suite — Opus 5.5 (medium) for design, Sonnet 5.5 (medium) to build

Depends on: T.1, T.11, S.3. Design (Opus): the runner, the comparison rules and the bless workflow. Build (Sonnet): the report and command-line tooling. The comparison rules are what decide whether every later PR "changed nothing", so they need Opus.

**Data directory** (`$KCWI_REGRESSION_DATA`, populated by hand, no download script, Q12):

```
manifest.json      datasets, file checksums, settings, tolerances, provenance, bless history
raw/<dataset>/     raw KOA frames
answers/           answers_<channel>_3.txt (S.3), later answers_<channel>_3_v2.txt (P.6)
reference/<dataset>/<plot_level>/   blessed outputs: redux/ and plots/
regions_<channel>.json              sky regions and standard aperture (M.0)
```

Each dataset entry in `manifest.json` holds: id, channel, description; raw file names with SHA-256; config overrides; command line; `plot_level`; profile (`quick` or `full`, D16); answers file; reference path; how the reference is compared (`self`, or `against: <dataset>` with tiers 3–4 informational, used by the `*_matrix` variants until their first bless); provenance (DRP commit, date, environment file, IERS table versions per D18); per-product tolerances (exact by default); and a history of blesses with reasons. Initial datasets: S.3 blue and red at `plot_level` 1 (quick) and 0 and 3 (full). Later variants are added as their own entries: `*_matrix` with `cube_method = "matrix"` (M.2e), and `*_warp` pinned with `cube_method = "warp"` (M.2c, D17). Raw paths must not contain `kr` or `kb` outside the file names (bug B8 would insert extra prompts until PR3b).

**Running.** `pytest -m regression` runs the quick profile, one test per dataset and check, so `-k blue` or `-k structure` selects a subset; `--regression-profile=full` adds the full-profile datasets. Each dataset is reduced once per session (a session-scoped fixture), by running `reduce_kcwi` as a subprocess in a temporary directory with the plot viewer disabled (needs T.11 on pre-PR1 code) and `astropy` IERS auto-download off (D18). Level-3 runs take stdin from the answers file (with `KCWI_STRICT_INPUT=1` once Track P lands); level-0 and level-1 runs take stdin from `/dev/null`, so an unexpected prompt fails instead of hanging. Raw-file checksums are verified first, and a mismatch stops that dataset before it runs.

**Prompt record (D15).** The runner puts a directory containing a `sitecustomize.py` on the subprocess's `PYTHONPATH`. It wraps `builtins.input` so that every call appends `{prompt, answer or "EOF"}` as one JSON line to `$KCWI_PROMPT_RECORD`. This is the only reliable way to check that an answers file was used exactly: Python reads a redirected stdin file in 8 KB blocks, so the file position shows nothing. It also detects changed prompt text. Track P's helper calls `builtins.input`, so the record works before and after PR3a. The record is kept in the report.

**Checks**, in tiers:

| Tier | Check | Fails when |
|---|---|---|
| 1 Completion | No `Failed executing primitive` or `Exception … while invoking` lines; no other `ERROR` log lines beyond an allow-list; the prompt record has exactly as many entries as the answers file has lines, no `EOF` entry, and (until a planned prompt change is blessed) the same prompt strings as the reference record (level 3); the prompt record is empty at levels 0 and 1. The exit code is recorded but not checked: it is always 0 (§6.1) | Any of these |
| 2 Structure | Same set of output files; same HDU names, shapes, dtypes; same header keyword set (allow-list: `DATE*`, `HISTORY`, DRP version, file paths; the `.proc` table's `filename` column is compared by base name) | Any difference |
| 3 Pixels | Per product and extension: exact for integers, flags and masks; floats within the manifest tolerance (exact unless stated) | Any value outside tolerance. Report gives max abs/rel diff, fraction differing, location |
| 4 Header values | Numeric keywords within tolerance; strings equal | Any difference |
| 5 Science metrics | The matrix plan Stage 0 metrics (arc centroids and widths, contbar positions, standard flux, invsens ratio, σ calibration in `regions_*.json`), with the thresholds from M.0 | A threshold missed |
| 6 Plots | Same PNG names and pixel sizes; after PR1, a pixel difference against the post-port reference | A missing or extra file, or a size change |
| 7 Timings | Wall time per primitive | Never (informational, Q14) |

This absorbs the regression script from matrix plan Stage 0 (tier 5) and the replay checks from the bokeh plan and Track P (tier 1). Bokeh-plan 4.1's side-by-side page can be generated from the tier 6 data.

**Report.** A Markdown summary plus diagnostic PNGs (difference images for failing extensions, metric plots), written to `$KCWI_REGRESSION_OUT/<timestamp>/` (default: a temporary directory, printed at the end), with a summary table ready to paste into a PR.

**Bless.** `python -m kcwidrp.tests.regression bless --dataset <id> --reason "<text>"` copies the latest run's outputs (and prompt record) over the reference and appends the old and new commit, the date, the reason and the diff summary to the manifest history. It never happens automatically. A bless is allowed only for a reason the plan expects, quoted in the PR:

| Reason | When | Datasets |
|---|---|---|
| `environment` | T.7 first run (below), and the PR1 environment switch if it changes outputs | All |
| `T.11 viewer-off` | PR0, if references were captured with the viewer on and differ only in what the viewer affects | All |
| `bokeh port 4.3` | After 4.3 sign-off | Plots (tier 6) |
| `prompt text PR3a` | After G3 | Level-3 prompt records only; outputs must not change |
| `P.5 <bug ID>` | After G4 | Datasets changed by that fix |
| `M.2c header` | After G5: `CUBEMETH` (and any agreed coverage extension) added to warp cubes | All warp datasets; pixel data must not change |
| `M.2e initial`, `M.3`, `M.4`, `M.5` | After the stage's 👤 sign-off | `*_matrix` only |
| `matrix default` | After G6 | Default-config datasets; `*_warp` stay as they are (D17) |
| `release version` | R.5 | Only if a header value outside the allow-list carries the version |

**Environment drift.** Library differences (numpy, scipy, astropy, IERS tables) can change floating-point output with no code change. S.3's references come from the baseline environment, so before any code change the suite is run on `upgrade-start` code in the baseline env with the test extra installed (D14); any differences from S.3 are reviewed and, if explained by the environment, blessed with reason `environment`. When PR1 switches to the development environment, the suite's tiers 1–5 are run on the PR1 branch; a difference that is not explained by the port is treated as a port bug, not blessed.

Validation:
- Run against its own references: everything passes.
- Scale one extension of one output by (1 + 1e-6): tier 3 fails, naming the file and extension.
- Delete one output file: tier 2 fails. Change one header value: tier 4 fails.
- Change one byte of a raw file: the checksum check fails before the reduction starts.
- Make one primitive raise (a test build only): tier 1 fails and names the primitive, although the process exits 0.
- Delete the last line of a level-3 answers file: tier 1 fails on the prompt-record count or the `EOF` entry. Add one extra line: tier 1 fails on the count.
- With `KCWI_REGRESSION_DATA` unset: every test is skipped with the reason.
- A bless records its reason, and a bless without `--reason` is refused.
- The full suite's wall time is recorded in `docs/testing.md`.

👤 Human check:
- **Assemble the data directory** from S.3: copy the raw frames, the answers files and the reference outputs. Sonnet can draft `manifest.json` (with checksums) from the directory; you check every entry.
- **Run the suite on unmodified code** and confirm it passes. Review any environment drift and do the "environment" bless.
- **Do every bless yourself.** A model may prepare the command but you run it.
- **Keep a second copy** of `$KCWI_REGRESSION_DATA` on another disk (§12).

### T.8 Suite quality check — Sonnet 5.5 (medium), fresh session

Depends on: T.3–T.7. Checks that the tests find real faults. In a scratch branch, apply at least 10 deliberate breakages one at a time, and record which test catches each:

overscan subtracted with the wrong sign; gain applied inverted; read noise left out of the uncertainty; a header keyword renamed in the FITS writer; mask dropped by the FITS reader; proctab search ignoring binning; an air→vacuum coefficient changed; DAR sign flipped; an event table "next event" misspelled; a config key renamed in one place only.

Then report coverage per module, against the T.1 baseline.

Validation: every breakage is caught by the fast suite (or by the regression suite, for the ones only real data can show, stated explicitly). A breakage that nothing catches gets a new test before PR0.

👤 Human check: review the coverage table and set a minimum per module for `core/` and `kcwi_file_primitives.py` (proposed: 70 % line coverage), to be written in `docs/testing.md`. Review any breakage that only the regression suite catches and decide whether that is acceptable.

### T.9 Testing guide and pre-PR checklist — Sonnet 5.5 (low)

Finish `docs/testing.md`: installing the test dependencies (`pip install -e ".[test]"`, Python 3.12); running the fast suite; setting up the regression data directory; running, reading and blessing regressions; writing a new test with `run_primitive` (a worked example); the pre-PR checklist that enforces the §1 rule. Propose a GitHub pull-request template with sections for the two test summaries, subject to Keck agreement.

Validation: a fresh Sonnet session, given only a clean clone and `docs/testing.md`, installs the test dependencies and runs the fast suite and (with a data directory) the regression suite, without other help.

### T.10 Remove Travis CI — Sonnet 5.5 (low)

Depends on: nothing. Q15: all Travis configuration and references go. Nothing consumes them (Q16), so this is a plain deletion with no settings to change. No replacement workflow is added in this upgrade.

Do, on the PR0 branch:
- Delete `.travis.yml`.
- `README.rst`: remove the "Tests" section's Travis link and status badge (lines 30–38 at `dfbcf9d`). Replace them with one sentence pointing to `docs/testing.md` for running the tests locally. Keep the ReadTheDocs badge.
- `environment.yml`: `.travis.yml` is its only consumer in the repo. Do not delete or change it here, but note in the PR that its fate is decided in bokeh-plan 5.1 (§7 amendment 9).

Hand-off to the documentation plan (S4-Q5; not done in this upgrade): the user documentation in `KCWI_DRP-Docs` should replace the pull-request requirement "A pass from the automatic TravisCI tests" (if it still carries it) with the §1 rule: fast-suite and (for reduction-code changes) regression-suite summaries in the PR description, with a link to the testing guide. R.4 collects this item.

Validation:
- `git grep -i travis` on the PR0 branch returns hits only in `docs/` planning documents and `changelog.txt`.
- `README.rst` renders without warnings (`python -m docutils README.rst /dev/null`, or `twine check` on a built sdist).
- The changelog entry states that Travis CI is removed, that there is currently no CI, and how to run the tests locally.

### T.11 Fix viewer-off crashes (B9, B10) — Sonnet 5.5 (low), tests in a fresh session

Depends on: T.1. D11. The only reduction-code change in PR0, and it changes nothing when a Bokeh server is running.

Do:
- Add `init_plot_session(context)` to `kcwidrp/core/bokeh_plotting.py`: it sets `context.bokeh_session = None` unless the attribute already exists. Call it from `reduce_kcwi.py`, `kcwi_rti.py` and `smart_reduce_kcwi.py` right after the framework is created and before any event is queued, so primitives can always read the attribute (B9). One helper rather than three assignments, so the scripts cannot drift apart (as bokeh plan 3.2 does for the viewer switch) and so B9 can be unit-tested; `main()` in the scripts cannot be run in pieces.
- `bokeh_clear(session)` in `kcwidrp/core/bokeh_plotting.py` returns at once when `session is None`, as `bokeh_plot` already does (B10).
- Changelog entry: `enable_bokeh = False` now works at every `plot_level`.

Validation:
- Tests written first, from §13 B9/B10 and this section only, in a fresh session, and seen to fail before the fix:
  - B10: `bokeh_clear(None)` returns without raising. (The fuller test, `SubtractOverscan` under `run_primitive` at `plot_level = 3` with `bokeh_session = None`, is added in T.4 once the T.2 harness exists.)
  - B9: `init_plot_session` on an object without the attribute sets it to `None`, and on an object that has one leaves it unchanged; and a `config_cli/` AST test checks that each of the three scripts calls `init_plot_session(framework.context)` before its first `append_event`/`ingest_data` call.
- 👤 After S.3 (needs data): a full reduction of one S.3 channel at `plot_level = 1` with `enable_bokeh = False` has no `Failed executing primitive` lines and gives outputs identical to S.3's (`array_equal` on every extension). This check runs as part of T.7's first run, which uses the viewer off.
- `git diff` touches only the lines above, the tests and the changelog.

PR0 = T.0–T.11 plus the planning documents (S.1), gate G0T.

### 6.3 Later: GitHub Actions (not in this upgrade)

No workflow is added now (Q15). The suite is built so that adding one later is a small, separate PR:
- The fast suite needs no display, network or data outside the repo (T.1), and excludes regression tests by default (`-m "not regression"`), so a workflow can run plain `pytest`.
- Test dependencies are a pip extra (D13), so a workflow installs them with one command (`pip install ".[test]"`).
- The regression suite needs `$KCWI_REGRESSION_DATA` and stays manual, unless Keck later provides a runner with access to the data.
- When a workflow is added, the §1 rule and `development.rst` requirement are updated to point to it, and the PR template from T.9 is adjusted.

### 6.4 How the other tracks use Track T

| Track | Uses | Adds |
|---|---|---|
| A | Fast suite on every commit; regression tiers 1–4 confirm reductions unchanged; tier 6 for plots | Bokeh-plan 1.2 tests in `unit/plotting/` (replacing `test_bokeh_plotting.py`); the coverage scan; bless of the plots after bokeh 4.3 |
| P | `run_primitive(answers=…)` for every prompt site; tier 1 replay and prompt record | P.2 tests in `unit/interactive/`; the bare-`input(` AST test in `unit/interactive/`; PR3a and PR3b remove the B-bug `xfail`s as in §13; `answers_*_v2.txt` added to the manifest; prompt-text bless after G3, output bless after G4 |
| M | Tier 5 metrics; T.4's `CorrectDar` and `WavelengthCorrections` tests as the warp-path reference | M.1t and M.2b in `unit/cube_matrix/`; `W` invariants in `properties/`: each column of the unnormalized overlap matrix sums to 1 for a pixel fully inside the grid (flux conservation), each row of the row-normalized `W` sums to 1 where coverage > 0 (weighted mean), variance non-negative, bilinear equivalence; `*_matrix` and `*_warp` dataset variants; blesses as in the T.7 table |
| R | R.1 checks the suite's hygiene | – |

---

## 7. Track A: Bokeh → Matplotlib

Execute [bokeh_port_plan.md](bokeh_port_plan.md) tasks 0.3, 1.1, 1.2, 2.1, 2.2, 3.1, 3.2, 4.1, 4.2, 4.3 and 5.1 as written, with these amendments. Task 5.1b is out of scope (amendment 13).

### Amendments

1. **Shared tasks.** Bokeh 0.1 → S.3; 0.2 → S.5; 0.4 → S.4; 0.5 → S.1. Bokeh 5.2 → R.1; 5.3 → R.2; 5.4 → R.3 (which also covers Track M).
2. **Models.** As in the remap table in §3.
3. **No prompt changes in Track A.** Prompt strings, order and guards stay byte-identical; Track P owns them. The one exception is already in the bokeh plan: the MakeInvsens "Hover the cursor…" instruction must be true of the new viewer's cursor readout (bokeh 3.1 verification, last bullet).
4. **Minimal MakeCube changes.** In bokeh 2.2, change only the imports, the `image` call support and (if 0.3 says so) the `Pool` start method. Do not refactor `make_cube_helper`; M.2c rewrites that code.
5. **The coverage scan (bokeh 1.2) becomes a permanent gate.** It stays in the test suite and runs on PR3, PR4–PR7. Any new plot added by Tracks P or M must use only wrapper-supported features or extend the wrapper with a test in the same PR.
6. **Matplotlib is the only plotting library in new code.** Validation scripts in Track M may use `matplotlib.pyplot` because they are not part of the package; package code follows wrapper rule 4 (no `pyplot`).
7. **PR1** is opened against `develop` after bokeh 4.3 sign-off; 5.1 may land in PR1 or as a follow-up, at the maintainers' preference (G1).
8. **Tests go into the Track T layout** (§6.4). Bokeh 1.2's `test_kcwi_plotting.py` lives in `kcwidrp/tests/unit/plotting/`. Bokeh 4.2's re-runs are regression-suite runs (full profile): tiers 1–5 must pass unchanged (cubes and the `.proc` table must not change in Track A), and tier 6 produces the plot comparison. Once 4.3 is signed off, the plots are blessed with reason "bokeh port 4.3". Bokeh 4.1's comparison page is generated from tier 6 rather than written separately. The bokeh plan's replay check ("ends with `EOFError`, or has unread lines left") is done with the T.7 prompt record, since an `EOFError` does not end the process and stdin position cannot be read (§6.1, D15).
9. **Travis is already gone when 5.1 runs.** T.10 deletes `.travis.yml` in PR0, so the bokeh-plan 5.1 bullet about it does not apply. Removing Travis also removes the only in-repo consumer of `environment.yml`; take that into account when 5.1 decides whether to fix or delete it.
10. **Test files named in the bokeh plan have moved** (T.1, T.6). Read them as:

    | Bokeh plan says | Use instead |
    |---|---|
    | `pytest kcwidrp/tests/test_import_primitives.py` (2.1, 2.2) | The `walk_packages` import test from T.1 |
    | `kcwidrp/tests/test_parse_configs.py` (5.1) | The config schema and reference tests in `kcwidrp/tests/config_cli/` (T.6) |
    | `kcwidrp/tests/test_bokeh_plotting.py` → `test_kcwi_plotting.py` (1.2) | `kcwidrp/tests/unit/plotting/test_kcwi_plotting.py` |

11. **The viewer-off fix from PR0 must survive.** T.11 makes the launch scripts set `context.bokeh_session = None`. Bokeh 3.1 keeps that attribute name; with `enable_plot_viewer` false no `StartPlotViewer` event runs, so the default from T.11 is what primitives see. Bokeh 3.1's "With `enable_bokeh = False` … all PNGs still written" check covers it.
12. **The starting point includes PR #218.** `v1.3.2` changed `kcwi_plotting.py` (non-snap Firefox/geckodriver preference). Bokeh 1.1 rewrites that file and 5.1 removes Firefox support, so the change is dropped deliberately; say so in the PR1 description and changelog.
13. **Maintainer answers (S4-Q1–Q5).** PNG consumers depend only on file names and format, so plot appearance may change (Q1). Remote viewing is over VNC with the browser on the pipeline host, so the `localhost` binding is right and R.3 tests over VNC (Q2). One tab per process is accepted and the shared-port alternative is dropped (Q3). `plot_viewer_open_browser` defaults to `False` in both `kcwi.cfg` and `kcwi_koarti.cfg`, which changes today's behaviour: users open the printed URL. Bokeh 3.1's start-up message must make the URL easy to find, and the changelog must say the tab no longer opens by itself (Q4). Bokeh 5.1b is out of scope: the user documentation is in the separate `KCWI_DRP-Docs` repository and is updated under its own plan, fed by R.4's hand-off list (Q5).

### Human checks in Track A (from the bokeh plan, collected here)

👤 Human check:
- **0.3 prototype:** the prototype's pass/fail list must be run on both macOS and Ubuntu 24.04. The browser-facing criteria (zoom, pan, cursor readout, tab close/reopen, two instances on two ports) need you at a browser.
- **3.1:** the interactive MakeInvsens session (cursor readout precise enough to type wavelength limits).
- **4.3:** review the side-by-side comparison page and the 4.2 table, step through one interactive `plot_level = 3` reduction, write the sign-off with accepted differences.

---

## 8. Track P: descriptive and robust interactive prompts

### 8.1 Goals and constraints

Goals (Q2, Q8):
1. Before each prompt, the user is told **what the plot shows, what to look for, what each answer does, and what the default is.**
2. Input is validated consistently: bad input re-prompts with a reason, `?` re-prints the full explanation, quitting is unambiguous, end-of-input is handled (D4).
3. Every answer that changes the reduction is logged at INFO with its effect.
4. Logic bugs at prompt sites (§13) are fixed, each in its own commit, with its effect on outputs measured.

Constraints:
- **D3:** for valid answers, the number and order of prompts is unchanged, so `answers_<channel>_3.txt` from S.3 replays on every commit of PR3a and produces identical outputs.
- **D12: no change for input that is valid today in PR3a.** Prompts, logged results and outputs for input accepted today are unchanged. Input that is invalid today (it crashes, silently resets, silently ends a review or quits by accident) may re-prompt instead; this fixes B1c, B2, B4, B5 and B6 in PR3a. Everything that changes a reduction result for valid input or for an unattended run, or removes a prompt that valid runs see, is in PR3b (B1, B1b, B3, B7, B8, B11, B12, B13).
- **EOF is a deliberate behaviour change in PR3a** (D4): today an `EOFError` fails the primitive and the reduction carries on without its output; after PR3a the default is used. It is listed in the changelog and covered by G3.
- **No plotting changes.** On-plot guidance was not selected (Q2). The helper must not import the plotting modules.

### 8.2 Design (to be fixed in P.1)

A new module `kcwidrp/core/kcwi_interactive.py`. Proposed API; P.1 may refine names but not semantics:

| Function | Replaces | Behaviour |
|---|---|---|
| `pause(config, logger, what=None)` | The 36 `if plot_level >= 2: input("Next? <cr>: ") else: time.sleep(plot_pause)` blocks | Same branching. When interactive, prints a one-line description (`what`) then `Next? <cr> (? for help)` |
| `ask_next(config, logger, what, help, n_items=None, current=None)` | `Next? <cr>, q to quit` and `Next? <int> or <cr>, q to quit` loops (ArcOffsets, SolveArcs, FitCenter, ExtractArcs, FindBars, MakeCube, MakeInvsens slice loop) | Returns `NEXT`, `QUIT` or a validated item number. A number outside `current+1 … n_items-1` re-prompts with the valid range (today it silently ends the review: bug B4, fixed here in PR3a per D12) |
| `ask_value(config, logger, prompt, help, parse, validate, default, effect)` | Value prompts: ReadAtlas offset, MakeInvsens limits / line mask / fit order, MakeMasterFlat peak, NandshuffSubtractSky scale | Shows current value and default; re-prompts on parse or validation error with the reason, keeping the previously accepted value (fixes B6); logs the accepted value with `effect`. The validator accepts exactly what today's code accepts as valid (D12); stricter checks that reject input valid today (sky scale 0, B7; very large fit orders, B3) wait for PR3b |
| `ask_yes_quit(...)` | `Proceed? <cr>=yes or Q=quit` in `reduce_kcwi.py` and `kcwi_rti.py` | Same choices; exact match on `q`/`quit` |

Rules:
- **Quit** is an exact, case-insensitive match on `q` or `quit`, after stripping whitespace. Today `'Q' in q.upper()` quits on any input containing a `q` (bug B5, fixed here in PR3a per D12; other input containing a `q` re-prompts).
- **Defaults** are what the code does today at `plot_level < 2`. For the MakeMasterFlat peak prompt that is the automatic highest-peak choice, which carries bug B1 until PR3b; the help text must not describe it as reliable before then.
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

Depends on: P.0. Implement `kcwidrp/core/kcwi_interactive.py` per §8.2. Opus because every one of the 60 call sites inherits its semantics, and the EOF/strict/replay interaction is easy to get subtly wrong.

Validation:
- `python -c "import kcwidrp.core.kcwi_interactive, sys; assert not any(m.startswith(('matplotlib','bokeh')) for m in sys.modules)"` passes.
- P.2 tests pass.
- Docstrings state the EOF and strict-mode behaviour.

### P.2 Helper tests — Sonnet 5.5 (medium)

Depends on: P.1. Written from §8.2, not from P.1's source. Use `monkeypatch` on `builtins.input` with scripted answers and `caplog` for logging.

Cover: each function with valid input; invalid then valid (re-prompt, reason printed, one answer logged); `?` then valid (help printed, not counted); exact `q`/`quit` vs. a value containing `q`; out-of-range item numbers; EOF default mode (one warning, default used, later prompts silent defaults); EOF with `KCWI_STRICT_INPUT=1` raises `EOFError`; `plot_level < 2` path calls `time.sleep(plot_pause)` and never `input`; INFO log line present for value answers and absent for "next".

Validation:
- At least three deliberate mutations of the helper (substring quit restored; EOF swallowed in strict mode; out-of-range accepted) each make a test fail.
- Suite runs in under 5 s with no terminal attached.

### P.3 Convert simple pause/next sites — Sonnet 5.5 (medium)

Depends on: P.1, P.2, PR1 merged. Convert every `pause` and `ask_next` site **except** those listed in P.4. Add each site's `what`/`help` text from P.0.

Validation:
- New test `kcwidrp/tests/unit/interactive/test_no_bare_input.py`: an AST scan that fails if `input(` or `time.sleep(...plot_pause...)` appears in `kcwidrp/primitives` outside the helper. (It starts with the P.4 sites on an allow-list, which P.4 empties.)
- `pytest kcwidrp/tests` passes, including the bokeh coverage scan.
- **Replay:** with `KCWI_STRICT_INPUT=1`, the full regression profile on both channels: tier 1 passes on the prompt record (same count; prompt strings differ only by the new `(? for help)` suffix and wording from P.0, which are reviewed and blessed as `prompt text PR3a` after G3), and every cube and the `.proc` table are `array_equal` to the PR1 output. PNGs are pixel-identical to PR1's.
- `git diff` per file touches only prompt blocks and imports.

### P.4 Convert value-entry and bar-jump sites — Opus 5.5 (medium)

Depends on: P.3. Sites: `ReadAtlas` offset; `MakeInvsens` limits, line mask, fit order and slice loop; `MakeMasterFlat` peak; `NandshuffSubtractSky` scale (×2); the bar-jump prompts in `ArcOffsets`, `SolveArcs`, `FitCenter`; the `Proceed?` prompts in `reduce_kcwi.py` and `kcwi_rti.py`. Opus because each site is interleaved with reduction logic, and a validator that is stricter than today changes which answers are accepted.

Rules (D12): **identical results for every input that is valid today.** Input that is invalid today may re-prompt (B1c out-of-range peak, B2 one-number line mask, B4 bar number out of range, B5 accidental quit, B6 reset of accepted limits). Bugs that affect valid input or unattended runs (B1, B1b, B3, B7, B8, B11, B12, B13) keep today's behaviour in P.4 and are fixed in P.5.

Validation:
- As P.3, and the AST allow-list is now empty. The B1c, B2, B4, B5 and B6 `xfail` markers are removed (all their sites are converted by now).
- Synthetic tests for each value site: the value reaching the reduction code is identical to today's for a table of valid inputs (from P.0).
- NandshuffSubtractSky: test on a synthetic non-standard N&S configuration (no real N&S data, Q4).

### P.5 Bug fixes — Opus 5.5 (medium) for fixes, Sonnet 5.5 (medium) for the failing tests

Depends on: P.4, PR3a merged. One commit per bug marked "PR3b" in §13 (B1, B1b, B3, B7, B8, B11, B12, B13, B16). For each bug, the failing test usually already exists as a strict `xfail` from T.4 or T.6. Where it does not, Sonnet first writes one from the bug description only. Opus then fixes the bug and removes the `xfail` marker in the same commit; the test passes.

Validation:
- Each commit message states the bug, the fix, and the measured effect on the S.3 reductions ("no change in outputs" or "invsens changes by X % at λ Y").
- After each output-changing fix, the regression suite is rerun and the diff explained. Once G4 is signed off, the changed datasets are blessed with reason "P.5 <bug ID>".

👤 Human check: approve or reject each fix individually. For B1 (MakeMasterFlat peak), look at the ledge and derivative plots for the affected flat and confirm the new peak is the right one. If S.3 never hit the buggy branch, confirm you accept validation on the synthetic test alone.

### P.6 Interactive acceptance session — Human

Depends on: P.4 (and again after P.5 for PR3b). Run one `plot_level = 3` reduction per channel in the development env, with the viewer open.

👤 Human check, at every prompt:
- Could a new user answer this prompt from the printed text alone, without reading code or the manual? Note any that fail.
- Try `?` at least once per prompt type; try one invalid input per value prompt; try a bar number below the current one and above the last.
- Try Ctrl-D at one prompt (expect: warning, defaults for the rest, reduction finishes).
- Afterwards, read the log: every value answer is there with its effect; no explanations are in it.
- Capture this session's input as `answers_<channel>_3_v2.txt` (it includes invalid inputs and `?`), and confirm it replays with `KCWI_STRICT_INPUT=1`. It becomes the regression input for all later PRs.
- After P.5 (PR3b), replay `answers_<channel>_3_v2.txt` again. If a fix changed the prompts it meets (B3 or B7 rejecting an answer it contains, B8 removing a prompt), re-record the file in this session and bless with reason `P.5 <bug ID>`.

### P.7 Changelog and docs — Sonnet 5.5 (low)

Changelog entries for PR3a (new prompt behaviour, logging, EOF handling, `KCWI_STRICT_INPUT`, the invalid-input fixes B1c, B2, B4, B5 and B6), committed in PR3a, and for PR3b (each fix and its effect), committed in PR3b. A short "Interactive mode" section for the documentation branch, held for R.4.

Validation: every behaviour change listed in §8.2 rules and every §13 fix appears in the changelog.

---

## 9. Track M: sparse-matrix cube construction

Executes [matrix_cube_plan.md](matrix_cube_plan.md) Stages 0–5.

### Amendments

1. **Stage 0 data** come from S.3 (D2). No nod-and-shuffle real data (Q4): N&S cube building is validated on a synthetic frame only.
2. **Stage 5** implements options A, B (opt-in) and C (Q6). Stage 6 is out of scope (Q3).
3. **Default flips to `matrix`** at the end of the track (Q5), as a separate PR (M.6) behind gate G6.
4. **Method guard.** Because the default changes, an invsens built with one method must not silently calibrate a cube built with the other (matrix plan §6, "never mix methods"). M.2d adds this.
5. **Concrete thresholds.** The matrix plan's acceptance criteria are qualitative ("small fraction of a pixel"). Proposed numbers are given below; 👤 you confirm or change them in M.0 before any validation runs, so the thresholds are not tuned to the results.
6. **Plots.** Validation figures are produced by the regression suite's report (outside the installed package), not by the pipeline. If a pipeline plot is added (e.g. a coverage map at `plot_level = 3`), it uses the wrapper and passes the coverage scan (Track A amendment 5).
7. **Line references** in the matrix plan were taken at `d0b4d7b`; the primitives are identical at `dfbcf9d` and `v1.3.2` (only `kcwi_plotting.py` differs), so they remain valid until Track A merges. By M.2 they will have moved: task prompts and reviews refer to code by function and symbol name, not line number.
8. **Wavelength and DAR folding apply to object frames only.** `WavelengthCorrections` and `CorrectDar` run only after `object_make_cube` and `nandshuff_make_cube`. In matrix mode, `arc_make_cube` and `flat_make_cube` build `W` from geometry alone (M.2), with no `WAVCORR` or `DARCOR`. "All four MakeCube events use the same method" (M.2c) means the same `cube_method`, not the same corrections.
9. **DAR breaks the per-slice block structure.** Once M.4 adds cross-slice weights, a pixel in slice *k* contributes to voxels in slices *k*−1, *k* and *k*+1, so `W` is no longer block-diagonal by slice. The matrix plan's "per slice, then `sp.vstack`/`sp.block_diag`" (§4) would only hold for M.2–M.3. So M.2a builds `W` from the start from per-*source*-slice blocks of global (row, column) indices that are summed, and the covariance in M.5 is computed per output slice with its neighbouring slices as a halo.
10. **Output files (S4-Q9, S4-Q11).** New products never add HDUs to existing outputs; each goes in its own file next to the cube: the coverage cube (M.2c), the option A correlation planes (M.5b) and the option B matrix (M.5b). The matrix file holds `W` in COO form as a FITS file with a single binary-table extension, replacing the matrix plan's `.npz`.

### M.0 Matrix fixture, metric thresholds and regions — Sonnet 5.5 (medium)

Depends on: T.7. Matrix plan Stage 0 (tests and metrics only; data from S.3). The regression script that Stage 0 describes is now tier 5 of the regression suite (T.7); M.0 supplies what that tier needs and the matrix-specific fixture.

Do:
- `kcwidrp/tests/unit/cube_matrix/conftest.py`: a fast synthetic fixture: a small detector, two or three slices, and **affine** per-slice transforms built as `skimage` `PolynomialTransform` objects so the exact answer is known analytically. Initially it tests only the fixture itself.
- Confirm that tier 5 computes the matrix plan Stage 0 metrics as M.2e–M.5c need them: arc-line centroid and width per slice vs `geom['avwvsig']`; contbar positions along slice; standard-star total flux and invsens ratio; σ calibration in sky regions (single spaxel and binned apertures of 1, 2×2, 3×3, 5×5 spaxels). Add any that are missing.
- **Contbar cube.** The pipeline never builds a contbar cube: the contbar chain ends at `contbar_trace_bars`. Tier 5 therefore builds one itself, outside the pipeline: run `MakeCube` through the T.2 harness on the master contbars frame from the reduction's `redux/`, with that night's `*_geom.pkl`, once with each `cube_method`. Document this in the metric's docstring.
- `regions_<channel>.json` in `$KCWI_REGRESSION_DATA`: sky-only regions and the standard-star aperture used by the σ and flux metrics.
- Add the thresholds below to the manifest as tier 5 thresholds for the `*_matrix` dataset variants.

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
- Tier 5 run on the warp datasets against their own references gives zero difference on every metric, for both channels (both slice geometries handled).
- The synthetic fixture's own tests pass in under 2 s.

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

Depends on: M.1b, PR3a merged. `cube_matrix.build_W(geom, pixcoords, gpm, out_grid, mode='overlap'|'bilinear')` per matrix plan §4: vectorized, no Python loops over pixels or voxels; built per source slice as blocks of global (row, column) indices that are summed, so M.4's cross-slice weights need no change of structure (§9 amendment 9); `ncand` per axis computed from data; mask then row-normalize; coverage kept; float32 values, int32 indices; flag bit-plane OR via the boolean structure of the unmasked overlap matrix.

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

Depends on: M.2a, M.2b, S.4 answers 8, 9 and 10a. Matrix plan Stage 2: config key `cube_method = "warp" | "matrix"` (default `warp` in this PR); in matrix mode, every MakeCube output (data, uncertainty, mask, flags, noskysub, obj, sky, delta-wave) via `W`; all four MakeCube events (`arc_make_cube`, `flat_make_cube`, `object_make_cube`, `nandshuff_make_cube`) use the same method (geometry only for arcs and flats, §9 amendment 8); red 180° rotation preserved; WCS/header code unchanged; new header keyword `CUBEMETH`, written in both modes. Coverage cube written as a separate file next to the cube (S4-Q9), e.g. `<frame>_icube_coverage.fits`, in matrix mode only.

Also add the `blue_warp` and `red_warp` regression datasets (D17): copies of the default datasets with `cube_method = "warp"` in their config overrides, referenced to the same outputs.

Validation:
- `cube_method = "warp"`: every output's pixel data, mask and flags are `array_equal` to PR3a's on S.3 data, and headers differ only by `CUBEMETH = 'warp'` (the switch alone changes nothing else). Tier 2 fails on the new keyword as expected; it is blessed as `M.2c header` after G5.
- `cube_method = "matrix"`: full reduction of both channels completes; every output file exists with the same HDU names, shapes and dtypes as the warp output, plus the separate coverage file; no existing file gains an HDU.
- Synthetic N&S frame runs through `nandshuff_make_cube` in matrix mode.
- `CUBEMETH` is present in every cube and propagates to `icubew`, `icubed`, `icubes`, `invsens`.
- `pytest kcwidrp/tests` passes, coverage scan included.

### M.2d Cube-method guard — Sonnet 5.5 (medium)

Matrix plan §6: invsens and science cubes must use the same method. In `FluxCalibrate`, compare `CUBEMETH` of the invsens and the cube. A missing keyword means `warp` (files from older versions). On mismatch (as confirmed by S4-Q12): log an ERROR naming both files and the fix ("re-derive invsens with cube_method = …") and skip flux calibration for that frame, rather than calibrate wrongly.

Validation: tests for matching, mismatching and missing-keyword cases; a real S.3 run with a warp invsens and a matrix cube triggers the guard.

### M.2e Real-data validation of M.2 — Opus 5.5 (high), then Human

Depends on: M.2c, M.2d. Add the `*_matrix` datasets to the manifest, compared `against` the warp datasets (T.7). Run both channels with `cube_method = "matrix"`, re-derive invsens with matrix, run the regression suite on the `*_matrix` dataset variants against the warp references (tiers 1–2 must pass; tiers 3–4 are expected to differ and are reported, not failed; tier 5 thresholds decide). Opus writes a first-pass report: each threshold pass/fail with numbers; any failure investigated to a cause (geometry convention, `ncand`, masking) before proposing a fix.

Validation: every M.0 threshold has a recorded value; failures have a stated cause and either a fix (with a test) or an explicit request to relax the threshold.

👤 Human check:
- Blink warp and matrix `icube` in ds9/QFitsView for one science frame and one standard per channel. Look for: slice-edge artefacts or stripes, ringing around cosmic rays and bad pixels (should be gone), holes where coverage should be full, red channel orientation matching blue.
- Look at the coverage cube: smooth, ~1 inside slices, falling to 0 only at slice and wavelength edges.
- Read the regression report and sign off thresholds, or record why a relaxed threshold is accepted.
- Once signed off, bless the `*_matrix` datasets with reason `M.2e initial`; from then on they are compared with their own references.

PR4 = M.2a–M.2e, gate G5. Default still `warp`.

### M.3 Fold wavelength corrections into `W` — Opus 5.5 (medium)

Matrix plan Stage 3. For `object_make_cube` and `nandshuff_make_cube` only (§9 amendment 8), compute radial-velocity and air→vacuum corrections from the header before building `W`, apply them to `lam_lo/hi`, set `WAVCORR = T`; `WavelengthCorrections` skips resampling when `WAVCORR` is set but still writes `VCORR`, `VCORRTYP`, `CTYPE3` etc. and still produces `icubew`. Reuse `a2v_conversion` and the velocity computation in `heliocentric`; do not duplicate the formulas.

Validation:
- The velocity correction computed in the helper equals the one `WavelengthCorrections.heliocentric` computes for the same header, to 1e-9 km/s, for both `heliocentric` and `barycentric` modes, and `none` gives zero.
- Sky-line centroids in matrix `icubew` vs. warp `icubew`: mean |Δ| ≤ 0.05 px; arc lines likewise.
- Line widths in matrix `icubew` ≤ warp × 1.05 (fewer resamplings).
- Header diff of `icubew` (warp vs. matrix) shows only expected keywords (`CUBEMETH`, `WAVCORR`).
- `air_to_vacuum = False` and `radial_velocity_correction = "none"` each run and give the expected uncorrected grid.
- Arc and flat cubes in matrix mode are `array_equal` to PR4's (no `WAVCORR`).
- The `*_matrix` datasets are re-blessed with reason `M.3` after the 👤 check below; `*_warp` and default datasets pass unchanged.

👤 Human check: plot one bright sky line's profile (warp vs. matrix `icubew`, same spaxel) and confirm the matrix profile is at the same wavelength and not broader.

### M.4 Fold DAR into `W` — Opus 5.5 (high)

Matrix plan Stage 4. For object and N&S frames only (§9 amendment 8). Per-λ along-slice and cross-slice offsets from `CorrectDar.atm_disper` and the same projection; padded output grid per grating as `CorrectDar` does; cross-slice weights handled by the per-source-slice assembly from M.2a (§9 amendment 9); `DARCOR = T`; `CorrectDar` skips shifts when set but keeps header updates, CRPIX padding (matrix plan open question 3), padding flags (value 128) and the `icubed` output. Opus at high because the WCS and padding must match exactly and cross-slice weights are new behaviour.

Validation:
- Standard-star centroid vs. λ (along and across slice): RMS about the mean ≤ the warp path's, both channels.
- RA/Dec of the standard-star centroid at three wavelengths (blue end, centre, red end) from matrix `icubed` WCS matches warp `icubed` within 0.1″.
- Padding region carries flag 128 and the same extent as warp.
- Standard-star flux within 1 %; invsens re-derived with matrix vs. warp within the M.0 threshold.
- CRPIX1/2 in matrix `icubed` equal warp's.
- With DAR folded in, the M.2a checks still hold: constant image → constant cube where coverage is full; the unnormalized column of an unmasked pixel sums to 1, now across up to three slices.
- The `*_matrix` datasets are re-blessed with reason `M.4` after the 👤 check below.

👤 Human check: in ds9, step through wavelength planes of the matrix `icubed` standard star and confirm it does not drift. Compare the two invsens curves on one plot and accept or reject the difference.

PR5 = M.3 + M.4 (may be split if review prefers).

### M.5a Covariance mathematics — Opus 5.5 (high)

Matrix plan Stage 5, options A, B, C (Q6). `W diag(var) Wᵀ` per output slice, with the neighbouring slices as a halo so cross-slice covariance from DAR is included (§9 amendment 9); option A extracts correlation coefficients with the +1 neighbour along each axis (3 planes; the −1 neighbour of a voxel is the +1 neighbour of the voxel before it, so 6 planes would repeat information); option B serializes `W` as a FITS file (§9 amendment 10; the per-pixel variance needed to use it is the uncertainty of the 2-D frame named in the file's header, not a copy); option C fits a binned-aperture noise-scaling curve and proposes header keyword names (for S.4 Q10b). Includes a documented formula for users of option B: `var(aᵀcube) = (Wᵀa)ᵀ diag(var) (Wᵀa)`.

Validation:
- Monte Carlo on the synthetic fixture: for apertures of 1, 2×2, 3×3, 5×5 spaxels, predicted variance from B equals empirical to within MC error (≤ 3 % at 1000 realizations); from A's neighbour approximation, within 10 %; C's curve within 10 % of empirical at each aperture size.
- Covariance build time per frame recorded; does not exceed 10 s on the S.3 machine (matrix plan measured ~6 s).

### M.5b Covariance I/O, config and tests — Sonnet 5.5 (medium)

Depends on: M.5a, S.4 Q10b. Config keys (names proposed, 👤 confirmed at G5b): `cube_covariance = True` (A and C), `cube_save_matrix = False` (B). Write A as a separate file next to the cube (S4-Q9), e.g. `<frame>_icube_corr.fits`; C as header keywords with the names agreed in Q10b; B as `<frame>_icube_W.fits` (S4-Q11): a routine `cube_matrix.write_matrix_fits(W, path, header)` writes `W` in COO form as a FITS file whose single binary-table extension has columns `ROW` (int32, voxel index), `COL` (int32, detector pixel index) and `VAL` (float32), with header keywords giving the matrix shape, the cube's voxel ordering and pixel ordering, the input 2-D file whose uncertainty holds the per-pixel σ, and `CUBEMETH`; `cube_matrix.read_matrix_fits(path)` returns a `scipy.sparse.coo_matrix`. Tests are written in a separate fresh session, independent of M.5a's and M.5b's source (§3 rule 1).

Validation:
- Files open with `astropy.io.fits`; `write_matrix_fits` → `read_matrix_fits` reconstructs B's matrix exactly (same shape, indices and float32 values); the matrix file has exactly one extension, a binary table.
- With both keys off, outputs are identical to M.4's.
- File-size increase per cube recorded for A and B.

### M.5c Covariance validation — Opus 5.5 (high), then Human

Using the S.3 sky regions: binned noise (1 to 5×5) vs. predicted from A, B and C. Report.

👤 Human check: confirm the reported vs. empirical noise curves look right on the sky regions; if CWITools (or another tool you use) reads covariance, run it on one cube and confirm it does not break on the new extensions.

After the 👤 check, the `*_matrix` datasets are re-blessed with reason `M.5`.

PR6 = M.5a–M.5c, gate G5b.

### M.6 Default flip — Sonnet 5.5 (low)

Depends on: M.5c and the S.4 Q8 and Q14 answers. `cube_method = "matrix"` default in `kcwi.cfg`, and in `kcwi_koarti.cfg` as Q14 decided; changelog entry stating that invsens must be re-derived with the new method and that `cube_method = "warp"` remains available (D5).

Validation:
- A reduction with the default config of both S.3 channels produces `CUBEMETH = 'matrix'` on every cube and invsens.
- Re-using an old (warp) invsens triggers the M.2d guard.
- The T.6 config schema and reference tests in `config_cli/` cover `cube_method`, `cube_covariance` and `cube_save_matrix`, including the allowed values of `cube_method`.
- The `*_warp` datasets pass against their own, unchanged references (D17); the default datasets are blessed with reason `matrix default` after G6.

👤 Human check: present M.2e/M.4/M.5c reports to the Keck maintainers for G6; record the decision in §2.

PR7 = M.6 + M.7 (M.7's report goes in the PR7 description, so PR7 opens after M.7).

### M.7 Performance and memory report — Sonnet 5.5 (low)

A table per channel: wall time and peak memory of GenerateMaps, MakeCube, WavelengthCorrections, CorrectDar, warp vs. matrix (with and without covariance). Goes in the PR7 description and the changelog summary.

Validation: numbers come from logged timings of real runs on the S.3 machine, not estimates.

---

## 10. Release

### R.1 Leftover sweep — Sonnet 5.5 (low)

Bokeh-plan 5.2 search list, plus: bare `input(` or `plot_pause` sleeps in primitives (should be none outside the helper); `TODO`/`FIXME` added by this work; config keys mentioned in docs but not in `kcwi.cfg` or vice versa. Test hygiene: no `xfail` without a bug or task ID, and none whose bug is marked fixed; no `skip` without a reason; no test outside the §6.2 layout; fast suite still under 60 s. Classify every hit as intentional or leftover.

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
- One `plot_level = 3` interactive session over VNC (S4-Q2), opening the printed viewer URL by hand (S4-Q4), replaying `answers_<channel>_3_v2.txt` with `KCWI_STRICT_INPUT=1`, then one live session answering by hand.
- `cube_method = "warp"` still runs end to end.

### R.4 Documentation hand-off list — Sonnet 5.5 (medium)

The user documentation lives in the separate `KCWI_DRP-Docs` repository and is updated under its own plan (S4-Q5). This task writes `docs/docs_handoff.md` in this repository: a checklist of every user-visible change that the docs plan must cover, each with the PR that introduced it and the source of truth (config comment, changelog entry, docstring). At least: bokeh-plan 5.1b's table (installation without Firefox/geckodriver/selenium/phantomjs; `enable_plot_viewer`, `plot_viewer_port`, `plot_viewer_open_browser = False`, `enable_bokeh` deprecation; the viewer URL and VNC viewing; viewer-port troubleshooting; `StartBokeh` → `StartPlotViewer`); the interactive-mode page (P.7: `?` help, quitting, EOF behaviour, `KCWI_STRICT_INPUT`, answers in the log); `cube_method`, `cube_covariance`, `cube_save_matrix` and the new output files (coverage, correlation, `W`); a short "cube construction" explanation and the invsens re-derivation note; API entries for `cube_matrix` and `kcwi_interactive`; the testing rule from T.10.

Validation: every changelog entry from PR0–PR7 that changes user-visible behaviour, a config key or an output file appears in the list.

👤 Human check: read the list and hand it to the documentation plan.

### R.5 Release — Human, with Keck maintainers

Run the full regression profile on the release candidate first. Version bump: `pyproject.toml` still says `1.3.1` although `v1.3.2` is tagged, so set the version from the tag history, not the file; because the default cube method changes science results, agree with Keck whether this is a minor or major version. Consolidated changelog, release notes listing: Bokeh removal and new viewer; new prompt behaviour; each bug fix with its effect; matrix default and the invsens re-derivation requirement; new keys and products. Gate G7.

---

## 11. Consolidated human checklist

| When | Task | What you do | Done when |
|---|---|---|---|
| Start | S.1 | Drop `stash@{0}`; fetch; confirm `develop` has `v1.3.2`; tag `upgrade-start` | Tag exists on a commit containing `v1.3.2`, tree clean |
| Start | S.2 | Repair the psutil install in the baseline env; confirm it runs Bokeh and opens the browser; create the Python 3.12 dev env later with matching numeric pins | Versions and test baseline recorded |
| Start | S.3 | Run blue and red reference reductions at levels 0, 1, 3; capture answers; archive outputs and screenshots; check replay and logs | `docs/baseline/README.md` written, replay identical, no failed primitives in logs |
| Start | S.4 | Ask maintainers questions 1–14 (10b after M.5a) | Answers in §2 (done 2026-10-22 except 10b) |
| After S.5 | S.5 | Spot-check 5 figure rows and 5 prompt rows | Noted in inventory |
| Track T | T.0 | Approve `docs/testing.md`, especially tolerance and bless policies | Approved in PR0 |
| Track T | T.7 | Assemble `$KCWI_REGRESSION_DATA` from S.3; check the drafted manifest; first run on unmodified code; "environment" bless if needed; keep a second copy of the data | Suite passes on `upgrade-start` |
| Track T | T.8 | Set per-module coverage minimums; accept or reject breakages only regression catches | Written in `docs/testing.md` |
| Every PR | §1 rule | Run fast and regression suites; paste summaries; do any bless yourself, with a reason | Summaries in the PR |
| Track A | A.0.3 | Run prototype criteria in a browser on macOS and Ubuntu | Decision written in bokeh plan |
| Track A | A.3.1 | Interactive MakeInvsens session with new cursor readout | Noted in PR1 |
| Track A | A.4.3 | Review comparison page; interactive level-3 session; sign-off | Sign-off with accepted differences |
| Before P.1 | P.0 | Review every prompt explanation for scientific accuracy | Each section marked approved |
| Track P | P.5 | Approve or reject each bug fix; inspect B1 plots | Decision per fix in PR3b |
| Track P | P.6 | Interactive session with invalid input, `?`, Ctrl-D; read the log; record `answers_*_v2.txt`; replay it again after PR3b | Notes in PR3a; v2 file replays on PR3a and PR3b |
| Before M.2e | M.0 | Choose sky regions and std aperture; confirm thresholds | `regions_*.json` and thresholds committed |
| Track M | M.2e | Blink warp vs. matrix cubes; check coverage; sign off report; `M.2e initial` bless | Sign-off in PR4 |
| Track M | M.3 | Compare one sky-line profile | Noted in PR5 |
| Track M | M.4 | Check std star does not drift; compare invsens | Noted in PR5 |
| Track M | M.5c | Check noise curves; downstream tool smoke test | Noted in PR6 |
| Track M | PR6 | Ask S.4 Q10b; present covariance keys and products for G5b | Decision in §2 |
| Track M | M.6 | Present evidence to Keck for G6 | Decision in §2 |
| Release | R.3 | Clean-machine run, remote interactive session, warp fallback | Checklist complete |
| Release | R.4–R.5 | Hand the docs list to the docs plan; release notes, version | Released |

---

## 12. Risks, known gaps and rollback

| Risk / gap | Mitigation |
|---|---|
| No nod-and-shuffle real data (Q4) | N&S validated on synthetic frames in P.4 and M.2c. 👤 Optional: fetch one N&S dataset from KOA before M.6 and run it; if not done, the release notes say N&S under `matrix` is validated on synthetic data only |
| Matrix default changes science results for all users (Q5) | Gate G6, method guard M.2d, `warp` kept (D5), thresholds fixed before results are seen (M.0) |
| Prompt wording may be a de-facto interface at Keck | S.4 Q6 before P.0 |
| Bug fixes change outputs at default `plot_level` (B1) | Separate PR3b, measured effect per fix, maintainer gate G4 |
| Long-running sequence; upstream `develop` moves | Rebase each track branch on `develop` before its PR; rerun the fast and regression suites after each rebase |
| No CI (Q11, Q15): tests can be skipped or rot | The §1 pre-PR rule, with both test summaries pasted into every PR; R.1 checks test hygiene; T.9's PR template if Keck agrees |
| Regression references go stale or are blessed carelessly | Bless needs a reason, is recorded in the manifest history, is only done by you, and only for the reasons listed in T.7 |
| Library upgrades change floating-point output with no code change | The "environment" bless step in T.7; provenance (environment file) stored with every reference |
| Regression data exist on one machine only | 👤 Keep a second copy of `$KCWI_REGRESSION_DATA` (it cannot be re-downloaded by script, Q12); the manifest's checksums detect corruption |
| Bokeh plan viewer decision branches (WebAgg or fallback) | Only Track A is affected; P and M are independent of the viewer |
| Memory for unbinned (1×1) frames in `build_W` | M.2a builds per source slice; R.2 checks full-frame unbinned memory explicitly |
| A primitive fails but the run "succeeds" (exit status always 0, §6.1) | Tier 1 scans logs for failed primitives (D15); S.3 checks its own logs |
| IERS table updates change velocity corrections and `icubew` with no code change | Auto-download off in tests and runner; IERS versions in provenance (D18) |
| Regression runs too slow to do for every PR, so they get skipped | Quick profile per PR, full profile only where required (D16) |
| Baseline env breaks again (it was already broken by a mixed psutil install) | `env_baseline.txt` recorded in S.2; rebuild instructions in `docs/baseline/README.md` |

Rollback: each PR is independently revertable. Until M.6, the matrix code is opt-in. After M.6, `cube_method = "warp"` restores old cubes without a code change. Track A rollback is as in the bokeh plan.

---

## 13. Bugs found while planning (input to T.11, P.4 and P.5)

Found by reading the code at `dfbcf9d` (unchanged at `v1.3.2`); B11–B13 were found by the S.5 inventory, B14–B15 by T.1's import test and B16 by T.2's attribute scan, all confirmed by reading the code. Each must be confirmed by a failing test before it is fixed. Line numbers will move after Track A. "Crash" below means the primitive raises: the framework logs the exception, that frame's chain of events stops, and the rest of the reduction continues (§6.1). "Fixed in" follows D11 and D12.

| ID | Where | Bug | Effect | Changes outputs? | Fixed in |
|---|---|---|---|---|---|
| B1 | `MakeMasterFlat.py:475-481` | When more than one peak is found and `plot_level < 2`, `ipk = properties['peak_heights'].argmax()` is an index into `peaks`, but is then used as a pixel index (`apk = xvals[ipk]`). It should be `peaks[argmax]`. | The flat ledge position is taken near the start of the trimmed range instead of at the highest peak. **This runs at the default `plot_level = 1`** | Yes, when multiple peaks occur | PR3b (P.5) |
| B1b | Same block | When **no** peak is found and `plot_level < 2`, `argmax()` of an empty array raises `ValueError`. | MakeMasterFlat fails; no master flat is written for that group, and later frames are reduced without it or fail to find it | Yes (crash → result) | PR3b (P.5) |
| B1c | `MakeMasterFlat.py:469-474` | Interactive peak: `isnumeric()` loop accepts any non-negative integer with no bounds check; out of range raises `IndexError`. `isnumeric()` is also true for characters such as `²`, which `int()` then rejects with `ValueError`. No way to accept the automatic choice. | Crash on typo | No for valid input | PR3a (P.4) |
| B2 | `MakeInvsens.py:469-480` | Line-mask entry with only one number: `qstr.split()[1]` raises `IndexError`, but only `ValueError` is caught. | MakeInvsens fails; no invsens for that standard | No for valid input | PR3a (P.4) |
| B3 | `MakeInvsens.py:655-680` | Fit order: `ford = int(qstr)` is assigned before the fit; if `np.polyfit` then fails (e.g. negative order), the error is caught but `ford` keeps the bad value while the fit is the old one. Very large orders are not rejected. | Logged/recorded fit order may not match the fit used | Possibly (logged value) | PR3b (P.5) |
| B4 | `ArcOffsets.py:92-99`, `SolveArcs.py:438-445`, `FitCenter.py:339-346` | "Next? <int>": a number ≤ the current bar, negative, or ≥ NBARS never matches the loop index, so the interactive review silently ends. | Confusing; user thinks they will see bar *n* | No | PR3a (P.4) |
| B5 | All `'Q' in q.upper()` sites | Any input containing `q` quits. | Accidental quit | No | PR3a (P.4) |
| B6 | `MakeInvsens.py:365-380` | An invalid limits entry resets *both* limits to WAVGOOD, discarding limits the user accepted earlier in the loop, without saying so. | Surprising reset | Only after invalid input | PR3a (P.4) |
| B7 | `NandshuffSubtractSky.py:185-195, 269-` | Sky scale 0 is accepted (only negatives re-prompt), giving no sky subtraction with no warning. `nan` and `inf` are accepted too. | Silent no-op sky subtraction | Yes, if 0 entered | PR3b (P.5) |
| B8 | `reduce_kcwi.py:436-469`, `kcwi_rti.py:412-440` | The channel check uses `'kr' in frame.lower()` on the full path, so a directory name containing `kr` or `kb` triggers a false "wrong channel" prompt; the prompt repeats for every file. | Spurious prompts, which also consume lines of a piped answers file | No | PR3b (P.5), because the fix removes prompts that valid runs see |
| B9 | `reduce_kcwi.py:405`, `kcwi_rti.py`, `smart_reduce_kcwi.py`; every `bokeh_plot(p, self.context.bokeh_session)` | With `enable_bokeh = False`, `StartBokeh` never runs, so `context.bokeh_session` is never set and every `bokeh_plot` call reached at `plot_level ≥ 1` raises `AttributeError` (e.g. `SolveArcs.py:531`). The MakeMasterFlat ledge figures (`MakeMasterFlat.py:441-454`) have no `plot_level` guard, so there it fails even at `plot_level = 0`. | Most plotting primitives fail; the reduction is unusable with the viewer off | Yes (crash → result) | PR0 (T.11) |
| B10 | `core/bokeh_plotting.py` `bokeh_clear`; called at `SubtractOverscan.py:143` | `bokeh_clear(None)` dereferences `session.document`. `StartBokeh` sets the session to `None` when the server is unreachable. | SubtractOverscan fails for every frame at `plot_level ≥ 3` without a server | Yes (crash → result) | PR0 (T.11) |
| B11 | `kcwi_rti.py:405-423` (`-f` branch) | After the wrong-channel prompt, answering `Q` empties the local `frames` list, but the code then ingests `args.frames`, so the quit is ignored. (`reduce_kcwi.py` ingests `frames` correctly.) | Quit at the channel prompt does nothing in RTI | No | PR3b (P.5), because the fix changes what a valid answer does |
| B12 | `NandshuffSubtractSky.py:303, 305` (non-standard N&S branch) | `sky[0, objrow0, :] = 0.` and `obj[0, objrow0, :] = 0.` use three indices on 2-D arrays (the intended slice is `[0:objrow0, :]`). | NandshuffSubtractSky fails after the user enters the sky scale, so non-standard N&S frames are never sky-subtracted | Yes (crash → result) | PR3b (P.5); synthetic N&S test only (Q4) |
| B13 | `MakeInvsens.py:404-408, 488-493` | Header lines of an existing `.lmsk` file are read with `strip()` and written back with no newline, so all header lines and the first mask line are joined into one line that starts with `#`. | On the next run with that file, the first line mask is silently dropped as a comment | Yes, on re-use of an edited `.lmsk` file | PR3b (P.5) |
| B14 | `kcwidrp/scripts/kcwi_masksky_ds9.py:2`, `pyproject.toml` | The script imports `regions` (since `2aba1ef`), but `pyproject.toml` declares `pyregion`, which nothing imports any more. | The `kcwi_masksky_ds9` entry point fails with `ModuleNotFoundError` in a clean install | No | PR1 (bokeh 5.1, which already edits the dependency list): replace `pyregion` with `regions` |
| B15 | `kcwidrp/scripts/reduce_server.py:11` | Imports `keckdrpframework.utils.easyHTTP`, which `keckdrpframework` 1.0.0 does not ship. Not an entry point. | The module cannot be imported or run | No | R.1: the maintainer decides whether to delete the script or fix it; the T.1 import test pins it as a strict `xfail` until then |
| B16 | `SolveAIT.py:219-224` | Reads `self.action.args.dich_frac`, but `GetAtlasLines.py:330` sets `args.dichroic_fraction` (as `SolveGeom.py:202` reads it). The `AttributeError` fallback always applies. | `DICHFRAC` in SolveAIT's output header is always 1.0; cube headers (from `geom['dich_frac']`) are unaffected | Yes (one header value) | PR3b (P.5), with the measured effect on S.3 outputs |

B1 is the only one known to affect unattended reductions at default settings. S.3's level-1 logs show whether it occurred on the reference data.
