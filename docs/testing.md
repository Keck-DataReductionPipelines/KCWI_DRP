# Testing the KCWI DRP

Status: **approved by the maintainer, 2026-10-22 (task T.0).** The sections on installing, running, regression data and writing tests were filled in by task T.9. Changes to these rules need the same approval.

This document sets the rules every test in `kcwidrp/tests/` follows. Each rule is written so that a reviewer can check it by reading a test. If a test breaks a rule, either the test is fixed or the rule is changed here first.

There is no CI. Tests are run locally before every pull request (see [Before opening a pull request](#before-opening-a-pull-request)).

## 1. Two suites

| Suite | Command | What it contains | Limits |
|---|---|---|---|
| **Fast** | `pytest` | Unit, property (invariant), config and CLI tests | Whole suite under 60 s on a laptop; no display, no network, no data outside the repository |
| **Regression** | `pytest -m regression` (quick profile), `pytest -m regression --regression-profile=full` | Full reductions of real datasets from `$KCWI_REGRESSION_DATA`, compared with blessed reference outputs | Run by hand; skipped with the reason `KCWI_REGRESSION_DATA not set` when the variable is unset |

`pytest` with no arguments never runs a regression test (`-m "not regression"` is the default in `pyproject.toml`).

## 2. Layout

```
kcwidrp/tests/
  conftest.py      shared fixtures: kcwi_config, synthetic_ccddata, tmp_redux, run_primitive
  harness.py       run_primitive(): runs one primitive without the framework queue
  unit/
    core/          functions in kcwidrp/core
    primitives/    primitives and module-level functions in kcwidrp/primitives
    plotting/      the plotting wrapper (Track A)
    interactive/   the prompt helper and the bare-input() scan (Track P)
    cube_matrix/   the sparse-matrix cube code (Track M)
  properties/      hypothesis-based invariant tests
  config_cli/      config files, event tables, entry points, script arguments
  regression/      runner, comparison, report and bless tooling (marker "regression")
```

A test file lives in the directory that matches the code it tests. A test outside this layout is a review failure (checked in R.1).

## 3. Markers

Registered in `pyproject.toml`; `--strict-markers` rejects any other marker.

| Marker | Meaning |
|---|---|
| `regression` | Needs `$KCWI_REGRESSION_DATA`; runs a full reduction. Only in `regression/` |
| `slow` | A fast-suite test that takes more than 2 s. Each use carries a comment saying why it cannot be faster. Expected to be rare |

## 4. What counts as a fast-suite test

A test in the fast suite:

1. writes only under pytest's `tmp_path` (or the `tmp_redux` fixture, which is inside it);
2. opens no network connection: `astropy` IERS auto-download is switched off in `conftest.py`, and a test that needs IERS or ephemeris data uses what ships with astropy;
3. needs no display: it runs with `DISPLAY` unset and opens no window or browser tab;
4. reads no file outside the repository and the test's temporary directory;
5. calls no real `input()` and no real `time.sleep()` (the `run_primitive` harness patches both; other tests patch them with `monkeypatch`);
6. finishes in under 2 s, or is marked `slow`;
7. passes when run alone, in any order, and in parallel with other tests (no shared module state is left changed).

## 5. Where expected values come from

Every expected value in a test carries a comment naming its source, using one of these labels:

| Label | Meaning | Example |
|---|---|---|
| `# source: independent` | Computed by a different, trusted implementation | `astropy`'s `radial_velocity_correction` for the velocity correction |
| `# source: analytic` | Derived by hand from the mathematics, with the derivation in the comment or the docstring | A known polynomial recovered by a fit |
| `# source: published` | Taken from a publication or standard table, cited | Filippenko (1982) for atmospheric dispersion |
| `# source: characterization @ <commit>` | Captured by running the code at a named commit, usually `upgrade-start` | A table of `bm_ledge_position` outputs |

A characterization value states what the code *did*, not what is *right*. It must never be presented as correct. If a characterization value is known to be wrong, the test is an `xfail` (section 7), not an assertion of the wrong value.

## 6. Tolerances

1. **Integers, flags, masks, booleans, strings, shapes, dtypes and HDU names are compared exactly.**
2. **Floating-point values are compared exactly by default** (`numpy.array_equal`, with `equal_nan=True` where NaN is expected).
3. A floating-point tolerance is allowed only with a comment on the same or preceding line that gives the reason and why the number is the size it is. Example:
   ```python
   # Tolerance: polyfit on 5 points, condition number ~1e4, float64 -> ~1e-12 relative.
   numpy.testing.assert_allclose(fit, expected, rtol=1e-10)
   ```
   "To make the test pass" is not a reason.
4. Tolerances are relative (`rtol`) for quantities that scale, absolute (`atol`) for quantities near zero, and both are stated when both apply.
5. **Regression comparisons start out exact** for every floating-point product. A per-product tolerance may be added to the regression manifest only with a written reason in the manifest entry, approved by the maintainer.

## 7. Known bugs: strict xfail

A bug that exists in the code is pinned by a test that asserts the **correct** behaviour and is marked:

```python
@pytest.mark.xfail(strict=True, reason="B2: one-number line mask raises IndexError")
```

- `reason` starts with the bug ID from the upgrade plan (§13) or a task ID (for example `A.3.2`).
- `strict=True` is mandatory: when the bug is fixed and the marker is forgotten, the test fails.
- The pull request that fixes the bug removes the marker in the same commit.
- An `xfail` without an ID, or whose bug is marked fixed, is a review failure (checked in R.1).
- `skip` is used only for a missing external resource (for example `KCWI_REGRESSION_DATA`) and always carries a reason.

## 8. Independence of tests from the code they test

1. Tests are written from the specification (the plan, the task description, docstrings, published behaviour), not from the implementation's source.
2. Tests are written by a different model, or in a fresh session with no access to the implementation conversation, from the code they test. The pull request says who wrote what.
3. A test that only re-executes the implementation's own formula (for example, computing the expected value with the same expression the code uses) is not a test. It is rejected in review.
4. Each new test file is checked once by deliberately breaking the code it tests; at least one test in the file must fail. The breakage used is noted in the pull request. Run these checks with `PYTHONDONTWRITEBYTECODE=1`, and delete `__pycache__` directories after reverting: a breakage of the same length as the original (for example `tan` → `sin`), reverted within the same second, leaves a stale `.pyc` that Python still loads, so later runs test the broken code.

## 9. Property (invariant) tests

- Live in `properties/` and use `hypothesis`.
- Each invariant is a statement that holds for the physics or mathematics, not a restatement of the code (for example "air-to-vacuum conversion is monotonic increasing", not "the function returns `a*x+b`").
- The default `hypothesis` profile is derandomized and small enough to keep the fast suite under 60 s. The `thorough` profile (`--hypothesis-profile=thorough`) is for manual runs.
- Each invariant was seen to fail once against a deliberately broken implementation (noted in the pull request).

## 10. Regression suite and blessing

Defined in task T.7 of the upgrade plan. The rules that matter for reviewers:

1. A regression run never changes a reference. Only `python -m kcwidrp.tests.regression bless --dataset <id> --reason "<text>"` does, and it refuses to run without `--reason`.
2. **Only the maintainer runs a bless.** A model may prepare the command.
3. A bless is allowed only for a reason listed in the upgrade plan's T.7 bless table (for example `environment`, `bokeh port 4.3`, `P.5 B1`, `matrix default`), quoted in the pull request with the diff summary.
4. Every bless is recorded in the manifest history: old and new commit, date, reason and diff summary.
5. A reduction's success is judged from its log and the prompt record, never from its exit status (the framework always exits 0).
6. Wall-clock timings are recorded in the report and never fail a test.

## 11. Coverage

Coverage is measured (`pytest --cov=kcwidrp`) and reported in PR0 against the baseline (`docs/baseline/coverage_upgrade_start.txt`, then `docs/baseline/coverage_t8.txt`).

**Minimum (set by the maintainer after T.8, 2026-10-24): 70 % line coverage per module** for every module in `kcwidrp/core/` and for `kcwidrp/primitives/kcwi_file_primitives.py`. A pull request may not take one of these modules below 70 %, and a new module there starts at 70 % or more. Exempt:

- `kcwidrp/core/bspline/` (vendored third-party code);
- `kcwidrp/core/kcwi_plotting.py` and `kcwidrp/core/bokeh_plotting.py`, until the Bokeh port (PR1) replaces them; the new plotting modules then follow the rule.

Other modules have no minimum; their coverage is reported, and a pull request should not lower it without saying why.

## Installing the test dependencies

Use Python 3.12. Python 3.13 does not work: `pyproject.toml` pins `numpy==1.26.4`, which has no 3.13 wheels. Create a conda environment or a venv, then install the package **editable** with the test extra, from the repository root:

```
conda create -n kcwidrp python=3.12      # or: python3.12 -m venv .venv && source .venv/bin/activate
conda activate kcwidrp
pip install -e ".[test]"
```

The extra adds `pytest` (8 or newer), `hypothesis` and `pytest-cov`.

Editable matters. A reduction (`reduce_kcwi`) and the tests must run the code in the working tree. A plain `pip install .` copies the package into `site-packages`, and then edits to the repository change nothing. Check which copy Python uses:

```
python -c "import kcwidrp; print(kcwidrp.__file__)"
```

The path must be inside your repository (`.../KCWI_DRP/kcwidrp/__init__.py`), not `site-packages`. To also see the installed version, which `kcwi_fits_writer` writes into the HISTORY of every output (it must equal `version` in `pyproject.toml`):

```
python -c "import importlib.metadata as m; print(m.version('kcwidrp'))"
```

A stale `kcwidrp.egg-info/` directory at the repository root (left by an old `setup.py`-style build) shadows the installed metadata whenever Python runs from the repository root, so the version above, and the version in every output's HISTORY, comes out wrong (an old build there reported 1.2.0 instead of the `pyproject.toml` version). It is untracked. If it exists, delete it: `ls -d kcwidrp.egg-info` should say "No such file". See `docs/baseline/README.md`.

## Running the tests

### The fast suite

From the repository root:

```
pytest
```

This runs `kcwidrp/tests/` without the regression tests (`-m "not regression"` is the default in `pyproject.toml`). It needs no display, no network and no data directory. It must finish in under 60 s (section 1), measured as the time pytest reports at the end of a run in an environment that has run the suite before. On the maintainer's laptop (Apple M4) that is about 32 s. The first run in a new environment takes longer (49 s reported, about a minute of wall time, in a fresh venv) because Python first compiles every module. A suite that gets close to the limit is a problem to fix, not to live with.

The run prints many warnings (about 1100: polynomial-fit `RankWarning`s, deprecations in pipeline code, `SyntaxWarning: invalid escape sequence` from docstrings). They are expected and are not failures; deprecations raised from `kcwidrp` code are deliberately left visible (only third-party ones are filtered in `pyproject.toml`).

Useful variations:

```
pytest -q                                            # short output
pytest kcwidrp/tests/unit/primitives/test_run_CorrectGain.py   # one file
pytest kcwidrp/tests/unit/primitives -k "red"        # a directory, tests whose name matches
pytest -rx                                           # also list every xfail with its reason
pytest --cov=kcwidrp --cov-report=term-missing       # coverage (section 11)
```

For the check of section 8.4 (break the code on purpose, see a test fail, revert), run with the byte-code cache off and clean up afterwards:

```
PYTHONDONTWRITEBYTECODE=1 pytest kcwidrp/tests/unit/primitives/test_run_CorrectGain.py
find kcwidrp -name __pycache__ -not -path "kcwidrp/__pycache__" -exec rm -rf {} +   # after reverting the breakage
```

(The tracked file `kcwidrp/__pycache__/__init__.cpython-36.pyc` is part of the repository, which is why the command above skips that directory. If `git status` shows it changed, `git checkout -- kcwidrp/__pycache__` restores it.)

The property tests use a small, fixed set of examples by default. For a longer, random search:

```
pytest kcwidrp/tests/properties --hypothesis-profile=thorough
```

This takes minutes, not seconds; use it when you add or change an invariant (section 9). It writes nothing to `.hypothesis/`.

### What the strict xfails mean

The summary line of a clean run looks like `710 passed, 8 deselected, 60 xfailed`. An `xfailed` test is a known bug (section 7): the test asserts the correct behaviour, fails today, and is marked `xfail(strict=True, reason="B<n>: ...")`. `pytest -rx` lists them with the reasons. They are expected and are not failures.

When one of them **passes**, pytest reports it as `FAILED ... [XPASS(strict)]`. That means the bug no longer exists in the code you are running, usually because your change fixed it. Do the following:

1. Confirm the fix is real and intended (read the test: it asserts the correct behaviour).
2. Delete the `xfail` marker in the same commit as the fix (section 7).
3. Mark the bug as fixed in the upgrade plan's bug list, as the plan's task requires.

If a bug was fixed by accident, or by a change that was not supposed to touch it, stop and find out why before removing the marker.

## Setting up the regression data directory

The regression suite reads everything from one directory, `$KCWI_REGRESSION_DATA`. It is assembled by hand (there is no download script) and lives outside the repository. Its layout:

```
manifest.json                datasets, SHA-256 checksums, settings, tolerances, provenance, bless history
raw/<dataset>/               the raw KOA frames: lev0/*.fits, files_flatlamp_<channel>.txt, kcwi.sky
reference/<entry>/           the blessed outputs of one entry: redux/, plots/, kcwib.proc or kcwir.proc
answers/                     level-3 answer files, answers_<dataset>_<channel>_3.txt
regions_<channel>.json       sky regions and standard-star aperture (added later, by the matrix work)
```

Other files and directories in the data directory (for example an archive of the S.3 screenshots and logs) are ignored by the suite.

`<dataset>` is the S.3 dataset name (for example `ds01_F22491_pt2_60266`). `<entry>` is `<short name>_<channel>_L<plot level>` (for example `ds01_blue_L1`). `raw/<dataset>` and `reference/<entry>` may be symlinks. To save disk, use symlinks for the raw frames and APFS clones for the references (`cp -c -R`, macOS): a clone shares the data blocks with the original and costs no space until a copy is changed.

Two rules for the path:

- **The path to the raw frames must not contain `kr` or `kb`** (bug B8: the reduction would then ask extra questions). This applies to the real, symlink-resolved location of the frames as well as to `$KCWI_REGRESSION_DATA`.
- Put the directory on a disk with room for one reduction's output (the runner refuses to start an entry with too little free space; see below).

### Worked example

Assume the S.3 run directories are in `~/drp/kcwi/s3/<dataset>/<channel>/L<level>/` (each with `redux/`, `plots/` and `kcwib.proc` or `kcwir.proc`), the raw datasets in `~/drp/kcwi/test_suite/suite/data/<dataset>/`, and the level-3 answer files in `docs/baseline/<dataset>/` of the repository. Run this from the repository root, once per dataset and channel (shown for ds01; repeat for `ds03_Compact_groupHS_60406`, and for `red` / `kcwir.proc`):

```bash
export KCWI_REGRESSION_DATA=~/drp/kcwi/regression_data
S3=~/drp/kcwi/s3
SUITE=~/drp/kcwi/test_suite/suite/data
DS=ds01_F22491_pt2_60266 ; SHORT=ds01 ; CH=blue ; PROC=kcwib.proc

mkdir -p $KCWI_REGRESSION_DATA/{raw,reference,answers}

# raw frames: a symlink to the dataset directory
ln -s $SUITE/$DS $KCWI_REGRESSION_DATA/raw/$DS

# references: clones of the S.3 run directories (L0, L1, L3)
for L in 0 1 3; do
  mkdir $KCWI_REGRESSION_DATA/reference/${SHORT}_${CH}_L$L
  cp -c -R $S3/$DS/$CH/L$L/redux $S3/$DS/$CH/L$L/plots $S3/$DS/$CH/L$L/$PROC \
        $KCWI_REGRESSION_DATA/reference/${SHORT}_${CH}_L$L/
done

# level-3 answer files
cp docs/baseline/$DS/answers_${DS}_${CH}_3.txt $KCWI_REGRESSION_DATA/answers/
```

(`cp -c` is the macOS APFS clone flag; on Linux use `cp -R` or `cp --reflink=auto -R`.) The `logs/` and `kcwi.cfg` of the S.3 runs are not needed in the data directory.

Then draft the manifest and check the raw frames:

```
python -m kcwidrp.tests.regression draft-manifest
python -m kcwidrp.tests.regression verify-checksums
```

`draft-manifest` reads `$KCWI_REGRESSION_DATA` (or `--data-dir`). It lists the datasets under `raw/`, computes the SHA-256 of every frame named in the `files_flatlamp_<channel>.txt` lists, makes one entry per dataset, channel and level, and reads the DRP version from the HISTORY of each reference cube. It refuses to overwrite an existing `manifest.json` without `--force`. `verify-checksums` re-reads every raw file and answer file and prints `ok` or the problems (exit status 1 if any).

### What to check by hand in the drafted manifest

The draft is a starting point. The maintainer reads every entry in `manifest.json`:

- **Entries.** One per dataset, channel and level: `ds01_blue_L0`, `_L1`, `_L3`, and the same for red. A level for which you made no reference still gets an entry; remove it, or leave it and bless it with reason `initial` after a first run (see below).
- **Profile and plot level.** `L1` is `quick`; `L0` and `L3` are `full` (decision D16). `plot_level` agrees with the level. A level-3 entry has an `answers` file.
- **Config overrides.** `enable_plot_viewer` is off and `plot_pause` is 0, so no viewer opens during a run.
- **Tolerances.** The drafted manifest has none: every float comparison is exact. A tolerance may be added only with a written `reason` in the entry, approved by the maintainer (section 6); the manifest loader refuses a tolerance without a reason.
- **Provenance gaps.** The draft prints what it could not find under "Could not determine / check by hand". Typical gaps: a reference with no `git version=` in its HISTORY (the S.3 ds01 runs were a plain `pip install .`, so only the version is recorded), `iers: null` (the IERS table versions are not determined until a bless records them), and references that do not exist yet. Fill in what you know; leave `null` rather than guess.
- **Checksums.** `verify-checksums` must say `ok` for every entry.

**Keep a second copy** of the whole directory on another disk. The data cannot be re-created by script, and the checksums only detect damage, not repair it (upgrade plan, section 12). With symlinked raw frames, copy the target files as well (`cp -R -L`).

## Running, reading and blessing regressions

### Running

Run everything below with the test environment active (`conda activate kcwidrp`, or `source .venv/bin/activate`), so that `pytest` and `python` are the environment's; `python -m pytest` works too. Set the variable in the shell where you run pytest:

```
export KCWI_REGRESSION_DATA=~/drp/kcwi/regression_data
pytest -m regression                          # quick profile: the L1 entries, both channels
pytest -m regression --regression-profile=full   # also the L0 and L3 entries
pytest -m regression -k ds01_blue_L1          # one entry
pytest -m regression -k structure             # one tier, all entries
pytest -m regression --keep-runs              # keep redux/ of passing runs
```

With the variable unset, every regression test is skipped with the reason `KCWI_REGRESSION_DATA not set` (`pytest -m regression -rs` shows it). The test ids are `<entry>-tier<N>_<name>`, so `-k` matches an entry name or a tier name.

Each entry is reduced once, by running `reduce_kcwi` as a subprocess in a fresh directory, with the plot viewer off, IERS downloads off, and stdin taken from `/dev/null` (levels 0 and 1) or from the answers file (level 3). The runner checks the raw checksums first; a mismatch stops that entry before it runs. Entries run one at a time, and **an entry's `redux/` is deleted after all its required tiers pass**, so the disk holds roughly one reduction at a time. Before starting, the runner checks that the free space is at least 1.5 times the size of the reference `redux/` (5 GB if there is no reference yet) and stops that entry if not. Failed runs are kept. Logs, the prompt record and the per-entry report are always kept.

Time: on the maintainer's Apple M4, a level-1 reduction takes about 3 to 5 minutes (193 to 288 s in S.3, `docs/baseline/README.md`), plus the comparison. The quick profile has four entries (two datasets, two channels), and took 12 minutes in the first run (2026-10-24: 132 to 172 s of reduction per entry, the rest comparison). The full-profile time is to be recorded here after its first run. Level-3 entries replay the interactive answers and take about the same as level 1. The runner gives up on an entry after 2 hours (`KCWI_REGRESSION_TIMEOUT`, in seconds, changes this).

Where the output goes: `$KCWI_REGRESSION_OUT/<timestamp>/` (default: `kcwi_regression_out/` under the system temporary directory). The directory is created if missing, and each run adds a new timestamped session directory; nothing is removed automatically, so delete old sessions yourself. A quick-profile session in which everything passed takes about 11 MB (the reductions' `redux/` directories are deleted); a failing entry keeps its `redux/` (about 2.5–3.5 GB per entry). The path of `report.md` is printed at the end. The run prints a few `SyntaxWarning`s from pipeline docstrings; they are expected. To run without pytest, with the same checks:

```
python -m kcwidrp.tests.regression run                      # quick profile
python -m kcwidrp.tests.regression run --profile full --entry ds01_blue_L1 --keep-runs
python -m kcwidrp.tests.regression report                   # rewrite report.md of the latest session from results.json
```

### Reading a report

The session directory holds:

- `report.md`: a summary table (one row per entry, one column per tier) ready to paste into a pull request, and below it the details per entry: run status, wall time, code version, and each tier's findings.
- `results.json`: the same results in machine-readable form (it is what `bless` and `report` read).
- `<entry>/`: the run itself: `reduce.out`, `kcwi.cfg`, `files.txt`, the procfile, `logs/`, `plots/`, `entry_report.md`, `prompt_record.jsonl` (only when the reduction asked any prompts, i.e. level 3), and `redux/` only if it was kept (a failing entry, or `--keep-runs`). `plots/` is always kept.
- Difference images (run minus reference, PNG) for the extensions that differ in tier 3; the report gives each file's path.

The table cells are `pass` (equal to the reference), `**FAIL**`, `info` (reported, never fails), `skip` (not applicable, for example no reference yet) and `-` (not run). The tiers, in plain words:

| Tier | What it checks |
|---|---|
| 0 checksums | The raw frames and the answers file are the ones in the manifest |
| 1 completion | The reduction finished: no failed primitive, no unexpected `ERROR` log line, and the prompt record matches the answers file (the number of answers used, no end-of-file; none at all at levels 0 and 1). The process exit status is shown but not checked, since it is always 0 |
| 2 structure | The same output files, the same HDU names, shapes and dtypes, the same header keywords |
| 3 pixels | Every pixel equal to the reference (integers, flags and masks always exactly; floats exactly unless the manifest says otherwise). The report gives the largest difference, the fraction of pixels that differ and where |
| 4 header values | Numeric header values within tolerance, strings equal |
| 5 metrics | Science metrics (arc line centres and widths, standard-star flux, and so on) against thresholds; none are registered yet, so it passes trivially until the matrix work (M.0) adds them |
| 6 plots | The same PNG file names and pixel sizes |
| 7 timings | Wall time per primitive, for information only; never fails |

Start with the first failing tier: tier 1 says the reduction itself broke (read `reduce.out` and the log), tier 2 that files or headers changed, tier 3 that numbers changed (open the difference image: a smooth offset, a few bad pixels and a shifted edge look different).

### Blessing

A bless makes a run the new reference. Rules (section 10): it is never automatic, **only the maintainer runs it** (a model may prepare the command), and it needs a reason from this list, written as `"<reason>: <text>"`:

`initial` (an entry with no reference yet), `environment`, `T.11 viewer-off`, `bokeh port 4.3`, `prompt text PR3a`, `P.5 <bug ID>` (for example `P.5 B2`), `M.2c header`, `M.2e initial`, `M.3`, `M.4`, `M.5`, `matrix default`, `release version`.

The upgrade plan's T.7 bless table says when each applies and to which datasets. Anything else is refused.

```
python -m kcwidrp.tests.regression bless --entry ds01_blue_L1 \
    --reason "environment: numpy 1.26.4 -> 2.x changes the last bit of invsens"
```

It uses the latest session under `$KCWI_REGRESSION_OUT`; `--run <session or run dir>` picks another. It prints the reason, the old and new code versions and the tier summary, then asks you to type `yes`. It refuses, changing nothing, if:

- `--reason` is missing, is not on the list, or has no text after the colon;
- the run did not complete, or tier 0 or 1 failed (a failed reduction cannot become a reference);
- the run's `redux/` was deleted because the run passed: run again with `--keep-runs`;
- the entry was not in that session.

On success it replaces `reference/<entry>/` (redux, plots, procfile, prompt record), writes the new provenance, and appends the date, reason, old and new provenance and diff summary to the entry's `bless_history` in `manifest.json`. Quote the reason and the diff summary in the pull request.

**The "environment" bless.** After the first run on unmodified code, and after any change to the libraries (numpy, scipy, astropy, the IERS tables), differences that the environment explains are blessed with reason `environment`. A difference the environment does not explain is a bug, not a bless.

**Adding a new dataset entry.** Put the dataset under `raw/<dataset>/` (frames, `files_flatlamp_<channel>.txt`, `kcwi.sky`), add its answers file to `answers/` if it has a level-3 entry, run `draft-manifest --output <new file>` and copy the new entries and datasets into `manifest.json` by hand (never `--force` on a manifest that has bless history: it overwrites it), and check the new entries as above. Run the new entry (`-k <entry>`, with `--keep-runs`): with no reference its comparison tiers are skipped and tiers 0 and 1 still run. If tiers 0 and 1 pass, bless it with `--reason "initial: <what the dataset is, which code produced it>"`.

## Writing a primitive test with `run_primitive`

Put a new test in the directory that matches the code (section 2). A primitive `Foo` gets `kcwidrp/tests/unit/primitives/test_run_Foo.py`; a function in `kcwidrp/core/` gets a file under `unit/core/`; an invariant goes under `properties/`. Follow the rules above it: expected values with a `# source:` label (section 5), exact comparison unless a tolerance has a written reason (section 6).

The `run_primitive` fixture runs one primitive without the framework: it builds the real `Action`, the context and the config, calls `Foo(action, context).apply()` in a temporary directory and returns the primitive's output (an `Arguments` object whose `.ccddata` is the frame). `synthetic_ccddata(camera=..., ampmode=..., binning=...)` builds a frame with the header keywords the primitives read. To chain primitives, pass the previous output as the next `ccddata`, as the framework does.

This example runs for real. It chains `TrimOverscan` and `CorrectGain` on a synthetic frame with a known answer, and pins a known bug with a strict xfail:

```python
import numpy as np
import pytest

from kcwidrp.primitives.CorrectGain import CorrectGain
from kcwidrp.primitives.TrimOverscan import TrimOverscan
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer


def test_correct_gain_doubles_a_flat_trimmed_frame(run_primitive,
                                                   synthetic_ccddata):
    ccd = synthetic_ccddata(camera="BLUE", ampmode="ALL")
    for amp in (1, 2, 3, 4):              # BLUE amplifiers are numbered 1-4
        ccd.header["GAIN%d" % amp] = 2.0
    trimmed = run_primitive(TrimOverscan, ccd)            # first primitive
    trimmed.ccddata.data = np.ones(trimmed.ccddata.data.shape)
    out = run_primitive(CorrectGain, trimmed)             # chained
    # source: analytic - every amplifier region is multiplied by its GAINn,
    # here 2.0 for all four, so a frame of ones becomes a frame of twos.
    assert np.array_equal(out.ccddata.data,
                          np.full(trimmed.ccddata.data.shape, 2.0))
    assert out.ccddata.header["GAINCOR"] is True    # source: CorrectGain docstring


@pytest.mark.xfail(strict=True, reason="B20: kcwi_fits_writer raises "
                   "KeyError for a header without any HISTORY card")
def test_writer_accepts_header_without_history(synthetic_ccddata, tmp_path):
    ccd = synthetic_ccddata()
    del ccd.header["HISTORY"]
    kcwi_fits_writer(ccd, output_file="x.fits", output_dir=str(tmp_path))
    assert (tmp_path / "x.fits").exists()      # source: analytic - it should write
```

Save it as, for example, `kcwidrp/tests/unit/primitives/test_run_CorrectGain_example.py` and run it with `pytest <that file> -rx`: one test passes and one is reported as xfailed. The xfail test asserts the **correct** behaviour; when the bug is fixed the test starts to pass, strict mode turns that into a failure, and the fix removes the marker (section 7).

Notes for real tests:

- **Prompts.** A primitive that calls `input()` (at `plot_level` 3, or an interactive review) needs `answers=[...]`, one string per prompt, in order: `run_primitive(Foo, ccd, config_overrides={"plot_level": 3}, answers=["", "q"])`. The harness fails the test, naming the prompt, if the primitive asks with no answers, asks more often than there are answers, or leaves answers unused. `run_primitive.last.prompts` lists the prompts asked. `config_overrides={"plot_level": 3}` changes the instrument config for that run. See `unit/primitives/test_run_ArcOffsets.py` for a full example.
- **A missing attribute.** A primitive reads attributes that, in a real reduction, earlier primitives or the pipeline set. If the harness does not provide one, the test fails with a message like

  `run_primitive: Foo read action.args.bar_avg, which the harness does not provide (in a reduction it is set by FindBars); pass it with action_extra={'bar_avg': ...}`

  Pass the value yourself: `action_extra={"bar_avg": 5}` for an `action.args` attribute, `context_extra={"arcs": [...]}` for a `context` attribute. The message names the primitive that sets it in a real reduction (the lists are `ARGS_NOT_PROVIDED` and `CONTEXT_NOT_PROVIDED` in `kcwidrp/tests/harness.py`). Give the attribute a value built from the specification, not copied from a real run's output. If the message says "not in the harness's documented lists", add the attribute to the right list in `harness.py` (and the test in `unit/test_harness.py` will tell you if you forget).
- **A frame the harness cannot ingest.** If the header lacks a keyword that the ingest helpers need, the harness raises `HarnessIngestError` naming the argument; set the keyword in the header, or pass the argument with `action_extra`.
- Expected values are computed by hand or by an independent routine, never with the formula the primitive uses (section 8). After writing a test file, break the primitive on purpose once and confirm a test fails (section 8.4), and note the breakage in the pull request.

## Before opening a pull request

Do all of the following, and put the results in the pull request description (`.github/pull_request_template.md` has a heading for each item; it is a proposal and can be dropped or changed by the maintainer).

1. Run `pytest` and paste its summary line into the pull request description.
2. If the pull request changes reduction code, run `pytest -m regression` (quick profile; full profile for pull requests that touch plotting or prompt code, and for releases) and paste the report's summary table. If it does not touch reduction code, say "not reduction code" and why.
3. A regression failure is resolved by a fix, or by a bless with a reason from the allowed list (section 10). Quote each bless in the description: the entry, the reason text and the diff summary. Only the maintainer runs a bless.
4. If the pull request fixes a known bug, the `xfail` marker of its test is removed in the same commit, and the description names the bug ID. No `xfail` is added without a bug or task ID in its reason, and none is left for a bug that is fixed.
5. No new `skip` unless a required external resource is missing, and it carries a reason (section 7).
6. Say who wrote the tests and who wrote the code under test (section 8): a different model, or a fresh session. For each new test file, name the deliberate breakage that made a test fail.
7. Run `pytest --cov=kcwidrp --cov-report=term-missing` and paste the lines for the modules the pull request touches. They must not fall below the minimum in section 11 (70 % per module for `core/` and `kcwi_file_primitives.py`, with the exemptions listed there).
8. Check `git status`: no stray files from the tests (`.hypothesis/`, `__pycache__` changes outside what is tracked, files in the repository root).
