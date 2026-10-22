# Testing the KCWI DRP

Status: **draft (task T.0)**, awaiting maintainer approval. Sections marked *(T.9)* are filled in when the test infrastructure exists.

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
4. Each new test file is checked once by deliberately breaking the code it tests; at least one test in the file must fail. The breakage used is noted in the pull request.

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

Coverage is measured (`pytest --cov=kcwidrp`) and reported in PR0 against the baseline. Minimum line coverage per module for `kcwidrp/core/` and `kcwidrp/primitives/kcwi_file_primitives.py` is set by the maintainer after T.8 (proposed: 70 %) and written here. *(T.8)*

## Installing the test dependencies *(T.9)*

```
pip install -e ".[test]"
```
Python 3.12. Details completed in T.9.

## Running the tests *(T.9)*

## Setting up the regression data directory *(T.9)*

## Reading a regression report and blessing *(T.9)*

## Writing a primitive test with `run_primitive` *(T.9)*

## Before opening a pull request

1. Run `pytest` and paste its summary line into the pull request description.
2. If the pull request changes reduction code, run `pytest -m regression` (quick profile; full profile for pull requests that touch plotting or prompt code, and for releases) and paste the report's summary table.
3. A regression failure is resolved by a fix, or by a bless with a reason from the allowed list (section 10).
