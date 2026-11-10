> Copied into this repository from `~/drp/kcwi/KCWI_docs_plan.md` on 2026-11-09 (upgrade plan D.1). This copy is now the source of truth.

# KCWI_DRP-Docs: audit and update plan

Date: 2026-10-29. Audited `KCWI_DRP-Docs` `main` at `5990deb` (2025-01-31) and the live site at kcwi-drp.readthedocs.io. Checked against:

- the released code, `KCWI_DRP` `origin/master` (1.3.3, still on Bokeh);
- the upgrade, `v3-upgrade` / `upgrade/download-fix` (PR0 and PR1 merged: Matplotlib viewer, no Firefox).

This is the "separate plan" for the docs that the upgrade plan names in S4-Q5. Its input from the upgrade is the R.4 hand-off list (`docs/docs_handoff.md`, not written yet).

Paths below are relative to the docs repo, so `source/` means `docs/source/`.

## 1. Summary

1. **The docs build is broken.** Every RTD build since at least 2026-10-22 fails at `pip install -r docs/requirements.txt`. The `latest` and `stable` pages still being served were built on 2025-02-01, and the next push to `main` will fail the same way.
2. **The API reference has been empty since before that.** The live `api.html` renders 3 of its 46 `automodule` entries: `kcwi_masksky_ds9`, `wb` and `wr`. None of the primitives render, and `check_cals` is missing from `scripts.html`. Autodoc fails to import them, but the build doesn't treat that as an error.
3. **The content describes 1.2.0.** 1.3.0 to 1.3.3 are not mentioned. Some pages are wrong even for 1.3.x (wrong config path, removed `setup.py`, a known issue that is already fixed).
4. **Large gaps.** The data-product table lists 14 of the roughly 35 file suffixes the pipeline writes. Nothing covers FITS extensions, flag values, units, nod-and-shuffle, the flux-calibration processing or most of the CLI. Issues #123, #163 and #199 ask for these docs.
5. **Bokeh-era content goes stale at the v3 release.** The released 1.3.x still uses Bokeh, so the Firefox and `enable_bokeh` content is correct today. It must not be removed from `main` until v3 ships.

## 2. Audit findings

### A. Build and infrastructure

| # | Where | Finding |
|---|---|---|
| A1 | `docs/requirements.txt` | Can't be resolved: `pandas~=1.0.3` conflicts with `keckdrpframework` from git HEAD, which needs `pandas>=2.0.0`. RTD builds 35334219, 35337984 and 35365851 (`develop`, 2026-10-22/23) failed with `ResolutionImpossible`. |
| A2 | `docs/requirements.txt:23` | Installs `kcwidrp @ ...@kcrm_merge`, a 2023 branch, so even a working build documents old code. The file also pins runtime packages the docs build doesn't need: `bokeh~=2.0`, `pyregion`, `pytest~=5.4`, `setuptools~=46`. |
| A3 | `.readthedocs.yml:16-18` | Builds with Python 3.7 on ubuntu-22.04, but `kcwidrp` requires Python >=3.12. This is the probable cause of A4. The file also has a leftover commented-out v1 config block (lines 1-5). |
| A4 | `source/api.rst`, `source/scripts.rst` | The automodule entries for the primitives and `check_cals` render nothing. Import failures only produce warnings because the build isn't run with `-W`. |
| A5 | `source/conf.py` | `release = '1.0'`; copyright says 2023. `sphinx.ext.autodoc` is listed twice. `analytics_id: 'UA-XXXXXXX-1'` is a placeholder, and `display_version` was removed from sphinx_rtd_theme 3.x. `html_theme_path = ["../.."]` and `templates_path = ['_templates']` point to nothing. |
| A6 | RTD versions | Active versions are `1.0.1`, `1.1.0`, `1.2.0`, `v1.2.0` (a duplicate of 1.2.0), `v1.3.1` (whose content still says 1.2), `develop`, `stable` and `latest`. No version has 1.3.x content. |
| A7 | `README.rst:6` | Broken link: `Keck-DataReductionPipeline/KCWI_DRP` is missing an "s". The README is also dated 01-Apr-2021. |
| A8 | `changelog.txt` | A single entry that stops mid-sentence ("...to avoid issue when"). |

### B. Errors in the current content (wrong for 1.3.x)

| # | Where | Finding | Fix |
|---|---|---|---|
| B1 | `quickstart.rst:13`, `running.rst:120` | Config path given as `kcwidrp/config/kcwi.cfg`. | `kcwidrp/configs/kcwi.cfg`; better, point to `--write_config`. |
| B2 | `configuration.rst:24` | Example `reduce_kcwi -r -c ... -f kb0001.fits` uses the red flag with a blue file. | `-b` with `kb…`, or `-r` with `kr…`. |
| B3 | `index.rst:28` | "Mac intel M2 and M3 chips". | Apple Silicon (M-series), as in changelog 1.2.0. |
| B4 | `issues.rst:40-47` | Known issue "Pipeline ingests blue files as red": two instances can't run at once. | Fixed in 1.2.0 (changelog). Remove it. |
| B5 | `issues.rst:49-107` | Firefox workaround: edit `kcwi_plotting.save_plot`, then `python setup.py install`. There has been no `setup.py` since the move to Poetry. `firefox_binary=`/`firefox_options=` were removed in Selenium 4. Partly superseded by PR #218 (1.3.2). | Replace with the current advice for 1.3.x (`conda install -c conda-forge firefox geckodriver`, the snap note). Delete it for v3. |
| B6 | `development.rst:73` | Requires "a pass from the automatic TravisCI tests". | Travis has been removed (T.10). Use the testing rule from `docs/testing.md` §1: fast-suite and regression summaries go in the PR description. |
| B7 | `development.rst:124,130` | Release steps say `python setup.py sdist bdist_wheel` and "version number in `setup.py`". | `python -m build` (or `poetry build`); the version lives in `pyproject.toml`. |
| B8 | `development.rst:41-63` | The branch model doesn't mention `v3-upgrade` (D28). It also says PRs merge "during the next developer meeting". | Update once v3 has merged; in the meantime a short note is enough. |
| B9 | `sky_subtraction.rst:127-129` | Says to enter `skip` "on the line after the raw file". | It's the second column of the same line (`MakeMasterSky.py:131-141`). |
| B10 | `sky_subtraction.rst:40-42` | `kcwi_masksky_ds9` is "contained in the `scripts` directory". | It's an installed command: `kcwi_masksky_ds9 <image>_intf.fits <regions>.reg`. The docstring's requirement that regions be saved in **physical** coordinates is missing from the page. |
| B11 | `running.rst:62-68` vs `sky_subtraction.rst:99-107` | The two pages give different re-reduction steps. One says delete the outputs and the proc entry; the other says set `clobber = True` and delete the proc entry. | `kcwi_pipeline.py` only skips a frame when it is in the proctab and `clobber` is False (issue #48). Write one procedure and link to it from both pages. Note B36: `kcwi.sky` rows aren't matched when frames are given with a path, until PR3b. |
| B12 | `outputs/masters.csv` | Lists `marcs` and `mdarks`. | The code writes `marc` and `mdark`. |
| B13 | `calibration.rst:10,23` | Links to the `kcrm_merge` branch. | Link to `master`. |
| B14 | `data_products.rst:17` | Sends readers to the IDL pipeline's `DATA_PRODUCTS` file for details. | Replace with a full products page (C3). |
| B15 | `index.rst:8-13`, `versions.rst` | Headline is "Release 1.2". `versions.rst` has no 1.2 entry and none for 1.3.x. | Write 1.2 to 1.3.3 entries from `changelog.txt`. |
| B16 | `scripts.rst:34-46` | The `check_cals` usage block omits `-a/--auto`. | Generate usage blocks with `sphinx-argparse`, or copy them from `--help`. |
| B17 | `updating_docs.rst:55` | Email address `keck.hawaii.edi`. | `.edu`. |
| B18 | Typos | "triggerd" (`updating_docs.rst:5`), "owndership" (`development.rst:105`), "Matuszewki" (`index.rst:88`, check the spelling), "a the official" (`index.rst:4`), "noticable" (`issues.rst:19`). | |
| B19 | `index.rst:67-69` | Asks users to "cite the KCWI DRP" but gives no citation. | Add the citation (KCWI instrument paper, plus an ASCL or Zenodo DOI if there is one). |

### C. Missing coverage (independent of version)

| # | Gap | Source of truth |
|---|---|---|
| C1 | **Most of the CLI is undocumented.** `reduce_kcwi` options `-i`, `-d`, `-m`, `-I`, `-w`, `-W`, `-s`, `-v` are missing, and monitor mode is commented out in `running.rst:162`. Decide which options are supported for users, and document those. | `reduce_kcwi.py:35-108` |
| C2 | **Undocumented scripts:** `modhead`, `smart_reduce_kcwi`, `start_kcwi_rti` (or a note that it is internal to Keck), and a `kcwi_masksky_ds9` usage block. | `pyproject.toml [tool.poetry.scripts]` |
| C3 | **Data products.** Missing suffixes include `invsens`, `std`, `ea`, `ocube*`/`scube*`/`dcube*` (N&S), `objf`, `skyf`, `delmap`, `trace`, `_geom.pkl`, `aitspec`, and the `saveintims` intermediates (`fsat`, `crr`, `def`, `gain`, `trim`, `warped`). There are no FITS extension descriptions (`UNCERT`, `FLAGS`, `MASK`, `NOSKYSUB`), no flag bit values (#163, #98), no per-product units or `BUNIT` (#123, #199), no key header keywords (`SKYCOR`, HISTORY provenance), and no proc-table format. | The `kcwi_fits_writer` calls (grep `suffix=`), `kcwi_file_primitives.py` |
| C4 | **Configuration reference.** `clobber` is used on the sky page but never described. `output_directory` and `procfile` are missing. The `[KCWI]`/`[BLUE]`/`[RED]` precedence rules aren't explained. | `kcwidrp/configs/kcwi.cfg` |
| C5 | **Processing overview.** Nothing says which recipe runs for which image type, what order group mode uses, or how many of each calibration are needed. There is no flow diagram. | `kcwi_pipeline.py` event table |
| C6 | **Nod and shuffle**: no page. | `NandshuffSubtractSky.py` |
| C7 | **Flux calibration processing.** `calibration.rst` covers observing a standard only, not how `invsens` is made and applied, or how to use a standard from another night. | `MakeInvsens.py`, `FluxCalibrate.py` |
| C8 | **Install versions.** PyPI's latest is 1.3.1, but 1.3.2 and 1.3.3 are tagged or merged, so `pip install kcwidrp` misses the StackFlats fix. Either publish to PyPI or tell users to install from git. *This is a release gap rather than a docs gap; it needs a maintainer action.* | PyPI |
| C9 | **Developer docs.** No link to `docs/testing.md` and no description of the test layout. | KCWI_DRP `docs/testing.md` |
| C10 | **FAQ items not yet written:** KOA versus original file names; why a frame was skipped (proc table); reducing the red and blue channels together; plot levels. | — |

### D. Stale once v3 is released (keep on `main` until then)

These come from the `Unreleased` changelog section and the R.4 scope in the upgrade plan.

| # | Where | Change |
|---|---|---|
| D1 | `installing.rst:47-56`, `issues.rst:49-107` | No Firefox, geckodriver, selenium or phantomjs. Delete these sections. |
| D2 | `quickstart.rst:15-21`, `configuration.rst:144-171` | `enable_bokeh` is replaced by `enable_plot_viewer`, `plot_viewer_port` (5006, then the next 9) and `plot_viewer_open_browser = False`. Users now open the printed URL themselves. Document the deprecation, viewing over VNC, SSH port forwarding (documented but not tested), and port troubleshooting. |
| D3 | `issues.rst:143-159` (PyQt5) | Check again under Matplotlib: the viewer uses the WebAgg path, so a Qt binding should not be needed. |
| D4 | `api.rst:115` | `StartBokeh` becomes `StartPlotViewer`; add `SendHTTP` and `SolveAIT`; decide whether `Template` stays. Later: `cube_matrix`, `kcwi_interactive`. |
| D5 | New page | Interactive mode (P.7): `?` help, quitting, EOF behaviour, `KCWI_STRICT_INPUT`, answers in the log. |
| D6 | New page + products | Cube construction: `cube_method` (and its default), `cube_covariance`, `cube_save_matrix`, the new products (coverage, correlation, `W`), and the note that invsens must be re-derived. |
| D7 | `versions.rst`, `index.rst` | v3 release notes; the first versions entry written in the new format. |
| D8 | `development.rst` | Dependencies live in `pyproject.toml` only (`requirements.txt` and `environment.yml` are gone); install with `pip install -e ".[test]"`; the testing rule. |

## 3. Plan

Each phase is one PR or a small set of PRs. Phases 0 to 2 target the docs repo's `main`, which serves the released 1.3.x. Phase 3 is written on a separate `v3` branch so that `main` stays correct for released users. Phase 4 merges everything into KCWI_DRP.

### Phase 0: Restore the build (docs repo, `main`)

1. `.readthedocs.yml`: ubuntu-24.04, Python 3.12, and delete the commented block.
2. Replace `docs/requirements.txt` with only what the docs build needs: `sphinx`, `sphinx_rtd_theme`, `sphinx-automodapi` and `sphinx-argparse`, plus `kcwidrp @ git+…KCWI_DRP.git@master`. The package's own dependencies install from its `pyproject.toml`.
3. `conf.py`: `release` read from `importlib.metadata.version("kcwidrp")`; remove the duplicate extension, the placeholder analytics ID and `display_version`; delete the dead `html_theme_path`/`templates_path` lines.
4. Build locally with `sphinx-build -W --keep-going -n`, so that a module that can't be imported fails the build. Fix A4 until every automodule renders.
5. Clean up RTD versions: hide the duplicate `v1.2.0` and deactivate `develop`.

Done when an RTD build of `main` passes and `api.html` lists all primitives.

### Phase 1: Correct errors for 1.3.x (docs repo, `main`)

Fix B1 to B19. Write the 1.2 to 1.3.3 entries in `versions.rst` from `changelog.txt`, and change the `index.rst` headline to the current release. Tag the result as `1.3.3` on RTD.

### Phase 2: Fill version-independent gaps (docs repo, `main`)

In priority order. Issue numbers go in commit messages so the issues can be closed.

1. Data products reference (C3): suffixes, extensions, flag values, units and headers. Closes the docs part of #123, #163 and #199. This needs a maintainer decision on the flag values first, because the code and docstring disagree (#98: saturation is flagged as 8, but the docstring says 1).
2. A single re-reduction procedure (B11) and the `kcwi.sky` format reference (columns, `skip`, `cont`, comments). Mention the B36 path limitation until PR3b ships.
3. CLI and scripts reference generated with `sphinx-argparse` (C1, C2, B16), so it can't drift from the code again.
4. Configuration reference (C4): every key in `kcwi.cfg`, grouped, with "do not change" marked. Add a check that compares the keys in `kcwi.cfg` with the keys documented, matching the R.1 sweep rule in the upgrade plan.
5. Processing overview and diagram (C5); nod-and-shuffle (C6); flux calibration (C7).
6. FAQ additions (C10).
7. C8: raise the PyPI gap with the maintainers. It's outside the docs, but the install page depends on it.

### Phase 3: v3 content (docs repo, `v3` branch, following the upgrade)

1. Branch `v3` from `main` after Phase 2. Point an RTD version at it (hidden) so testers can read it.
2. Write D1, D2, D3 and D4 (the viewer half) now, since PR1 has merged.
3. Write D5 after PR3a/PR3b, and D6 after the M track and the M.6 default-flip decision. The R.4 hand-off list is the checklist for this phase: tick every item in it.
4. On the v3 release (R.5), merge `v3` into `main` and write D7 and D8.

Done when every R.4 item is ticked and `-W` builds pass.

### Phase 4: Merge into KCWI_DRP (final step)

Decisions are recommended in §4.

1. **Prepare the docs repo:** freeze `main` and announce the move in its README.
2. **Import with history:**
   - On a clone, run `git filter-repo --path docs/ --path .readthedocs.yml`. This drops the docs repo's `README.rst`, `LICENSE.rst`, `CODEOWNERS` and `changelog.txt`, which would otherwise collide with KCWI_DRP's.
   - In KCWI_DRP, on a branch from `develop` (the target branch depends on timing, see the note below), run `git merge --allow-unrelated-histories`.
3. **Resolve the clash with the planning docs** in `docs/`, where `drp_upgrade_plan.md`, `testing.md`, `baseline/` and the others live:
   - Move them to `docs/dev/` and update the links (README's "Tests" section, the PR template, cross-links between the plans).
   - Sphinx sources stay in `docs/source/`.
4. **Packaging:**
   - Add a `docs` extra to `pyproject.toml` with the Sphinx packages, and delete `docs/requirements.txt`.
   - `.readthedocs.yaml` at the repo root: `python.install: [{method: pip, path: ., extra_requirements: [docs]}]`. The docs then always build against the code in the same commit.
5. **Bring the content into the repo:**
   - `development.rst` and `updating_docs.rst` describe one repository.
   - `development.rst` links to `testing.md`, or includes it via `myst-parser`.
   - Add "docs build with `-W`" to the pre-PR checklist in `testing.md`.
   - Add a test that compares the config keys in `kcwi.cfg` with the documented keys, in `tests/config_cli/`.
6. **Switch RTD:**
   - Point the `kcwi-drp` project at KCWI_DRP and set `latest` to `develop` and `stable` to the newest tag.
   - Turn on PR preview builds.
   - Check that the old docs versions are still reachable (§4, decision 1).
7. **Archive** `KCWI_DRP-Docs` on GitHub, with a README pointing to the new location.
8. **Close the loop in the upgrade plan:** R.4's hand-off list becomes "edit `docs/source/` in the same PR" for later work, and the plan's decisions record gets a new dated row for this change, following its own convention.

Done when RTD builds from KCWI_DRP, old versions are reachable, and the docs repo is archived.

**Timing changed 2026-11-09 (upgrade plan D32): Phase 4 now runs before Phase 3**, as Track D (§8a of the upgrade plan) and its own PR against `v3-upgrade`. Phase 3's content is then written in `docs/source/` by each upgrade PR, and there is no separate `v3` docs branch. The RTD switch and archiving (steps 6–7) wait for the v3 release.

**Original note on timing.** If Phase 4 is done before Phase 3 instead of after it, the upgrade PRs (PR3a, PR3b, M) can each carry their own docs changes. That removes the R.4 hand-off step and the separate `v3` branch. Merging last, as asked, is safe but means writing v3 docs out of tree.

## 4. Decisions needed

1. **Old RTD versions after the switch.** RTD versions are branches or tags of the connected repo. KCWI_DRP's tags `v1.1.0` to `v1.3.x` contain no Sphinx sources, so when the project is repointed, the 1.0.1 to 1.3.x docs would disappear on the next sync.
   - *Recommended:* a second RTD project (for example `kcwi-drp-legacy`) on the archived docs repo, linked from `versions.rst`.
   - Alternatives: attach static HTML zips to GitHub releases, or accept losing them.
2. **Where the planning docs go:** `docs/dev/` (recommended) or a top-level `dev/`.
3. **Markdown support:** add `myst-parser` so that `testing.md` can be part of the site. Recommended: yes, for developer pages only, and keep user pages in RST.
4. **Flag values (C3):** document what the code does, or what the docstring says, before #98 is fixed.
5. **Supported CLI options (C1):** which of `-d`, `-m`, `-i`, `-I`, `-w`, `-W` and `-s` are for users, and which should be hidden.

## 5. Progress

### 2026-10-29: Phases 0 to 2 done (local branches, not pushed)

The work is on three stacked branches in `~/drp/kcwi/KCWI_DRP-Docs`: `docs/phase0-build`, then `docs/phase1-fixes`, then `docs/phase2-content`.

- **Base branch.** The work starts from `origin/1.3.1`, not `main`. `1.3.1` holds 11 unmerged commits from 2026-03 (Poetry, Python 3.12, the Bokeh FAQ) and is a fast-forward of `main`.
- **Build.** It passes `sphinx-build -W` with Python 3.12 and kcwidrp `master` (1.3.3). `api.html` now shows 50 modules (it showed 3).
- **Not done (needs RTD admin access):** Phase 0 step 5. Hide the duplicate `v1.2.0` version and deactivate `develop`. After merging, point RTD's default branch, or `main`, at the new content.
- **Changes from the plan:**
  - The API page is back in the table of contents. Commit 8aaa9f9 on `1.3.1` had removed it.
  - The command-line reference uses build-time `--help` output instead of sphinx-argparse, because the scripts don't expose their parsers.
  - `conf.py` now fails the build when a `kcwi.cfg` key or a console script is undocumented.
- **Decisions taken by default (§4):**
  - Flag values are documented as the code sets them: 2, 4, 8, 64, 128.
  - The framework CLI options (`-i`, `-d`, `-m`, `-I`, `-w`, `-W`, `-s`) are listed but described as for real-time use only.
- **Please check:** the DRP citation. Only Morrissey et al. 2018 is given; add a DOI or ASCL entry if one exists.

### Code issues found while writing (for the issue triage)

These are true of `master` 1.3.3. Check each against `v3-upgrade`.

1. `kcwi_masksky_ds9` imports `regions`, but 1.3.x depends on `pyregion`, so the script fails in a clean install. Fixed in v3 (changelog).
2. `smart_reduce_kcwi` fails with an ImportError. Fixed in v3 (changelog).
3. `output_directory` is ignored by many primitives, which hard-code `'redux'`: FluxCalibrate, SubtractDark, MakeMasterSky, MakeMasterFlat, SubtractSky, ExtractArcs, MakeMasterBias, SubtractBias.
4. Group mode (`-g`) never processes `DARK` frames (`reduce_kcwi.py` `imtypes`).
5. For nod and shuffle, WavelengthCorrections only corrects `icube`. The `ocubed`, `scubed` and `dcubed` cubes have no radial velocity or air-to-vacuum correction.
6. `LINETHRESH`, `TUKEYALPHA` and `FRACMAX` in `[KCWI]` are always overwritten by `[BLUE]`/`[RED]`. Any other key in the channel sections is ignored.
7. The CLI overrides test truthiness, so `-o 0` and `-t 0` are ignored.
8. The comment above `get_unique_CCD_master_name` and `get_unique_STATEID_master_name` says "Currently NOT USED", but master biases and darks are named with them.
9. PyPI has 1.3.1, while `master` is 1.3.3 (C8).

### 2026-11-09: Phases 0–2 up for merge; §4 decisions 1–3 taken

- Phases 0–2 are in one PR into `main`: [KCWI_DRP-Docs#4](https://github.com/Keck-DataReductionPipelines/KCWI_DRP-Docs/pull/4), branch `docs/phase2-content` (pushed). Rebuilt before opening: `sphinx-build -W --keep-going` passes in a fresh Python 3.12 environment from `docs/requirements.txt` (kcwidrp `master`, which still reports 1.3.3); `api.html` lists 50 modules.
- Decisions (maintainer):
  1. Old versions: a legacy RTD project on the archived docs repository, linked from `versions.rst`.
  2. Planning docs move to `docs/dev/`.
  3. `myst-parser` for developer pages only; user pages stay reStructuredText.
- Next: Phase 4 now runs as Track D of the upgrade plan (§8a, D32), starting with D.1 once #4 is merged.
