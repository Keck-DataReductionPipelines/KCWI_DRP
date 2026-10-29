# Baseline environment and reference reductions

Tasks S.2 and S.3 of [drp_upgrade_plan.md](../drp_upgrade_plan.md). Only text files (`*.md`, `*.txt`, `*.json`) in this directory are committed (see `.gitignore`); FITS outputs, plots and logs from S.3 stay local.

## Baseline environment (S.2, 2026-10-22)

- Conda env `kcwidrp`, Python 3.12.13, macOS arm64 (`Darwin arm64`).
- `kcwidrp` installed editable (`pip install -e .`) from `upgrade-start`. That tag was first `9be28f5` (`v1.3.2`), moved to `6270cba` (1.3.3, the StackFlats hotfix, D21), and then to `af46243` (1.3.4, the MakeMasterFlat BM-ledge hotfix B1/B1b, PR #223) before S.3; the environment itself did not change.
- On 2026-10-23 a plain `pip install .` had replaced the editable install with a copy in `site-packages` (identical to 1.3.3). The editable install was restored the same day (`pip install --no-deps -e ".[test]"`); check with `python -c "import kcwidrp; print(kcwidrp.__file__)"`, which must point into this repository.
- A stale, untracked `kcwidrp.egg-info/` (version 1.2.0, from an old local build) at the repository root shadowed the installed metadata whenever Python ran from the repository root, so `importlib.metadata.version('kcwidrp')`, which `kcwi_fits_writer` writes into every output's HISTORY, reported 1.2.0 there. It was moved out of the repository on 2026-10-22. If it reappears (e.g. after a `setup.py`-style build), delete it before S.3 or any regression run.
- Full package list: [env_baseline.txt](env_baseline.txt). Key versions: numpy 1.26.4, scipy 1.14.1, astropy 6.1.7, scikit-image 0.26.0, astroscrappy 1.1.0, bokeh 2.4.3, matplotlib 3.9.4, keckdrpframework 1.0.0, psutil 5.7.3, pytest 5.4.3.

### Repair done before recording

The env had two psutil installs mixed in one directory: conda's psutil 7.2.2 (pulled in by `ipython`) and pip's 5.7.3 (pinned by `pyproject.toml`). `import psutil` failed (`_psutil_osx has no attribute getpagesize`), so `reduce_kcwi --help` failed. Fixed by `pip uninstall psutil` (twice, removing both), `conda remove --force psutil`, then `pip install "psutil>=5.7.0,<5.8.0"`.

Known, accepted: `pip check` reports `ipython 9.17.1 requires psutil>=7`. The baseline keeps the pin from `upgrade-start`; psutil is only used to list processes and does not affect reduction results.

To rebuild this env from scratch: create a Python 3.12 conda env, `pip install -r docs/baseline/env_baseline.txt` (skipping lines that fail on another platform), then `pip install -e .` from `upgrade-start`.

### Validation

- `python -c "import psutil, kcwidrp, bokeh, matplotlib"` succeeds.
- `reduce_kcwi --help` exits 0.
- Existing tests at `upgrade-start`:
  - `python -m pytest kcwidrp/tests` **fails to start**: pytest 5.4.3 (the runtime pin) is incompatible with Python 3.12 (`TypeError: required field "lineno" missing from alias` in assertion rewriting). This is a pre-existing failure, fixed by T.1 (pytest ≥ 8, D13).
  - `python -m pytest kcwidrp/tests --assert=plain -p no:cacheprovider`: **52 passed** (46 primitive imports, `test_import_pipeline`, three config-parse tests, `test_check_running_process`, and `test_import_kcwi_file_primitives`), 20 warnings, 9.7 s.
- 👤 Still to confirm by hand: the Bokeh browser tab opens from a `reduce_kcwi` run in this env.

## Reference reductions (S.3, 2026-10-23)

Run by the maintainer; checks below were run on the outputs the same day.

### Datasets (plan D24)

From the maintainer's test suite, `~/drp/kcwi/test_suite/suite/` (`manifest.csv` lists every KOA file name with its role):

| Dataset | Night | Channel | Grating | Slicer | Binning | Science | Standard | Sky frame |
|---|---|---|---|---|---|---|---|---|
| `ds01_F22491_pt2_60266` | 2023-11-18 | BLUE | BL, 4500 Å | Large | 2×2 | F22491_pt2, 1320 s | feige34 | yes |
| | | RED | RH1, 6140 Å | Large | 2×2 | F22491_pt2, 300 s | feige34 | yes |
| `ds03_Compact_groupHS_60406` | 2024-04-06 | BLUE | BM, 4775 Å | Medium | 2×2 | Compact groupHS, 200 s | hd93521 | yes |
| | | RED | RL, 7250 Å | Medium | 2×2 | Compact groupHS, 200 s | hd93521 | yes |

Science and standard frames are slow readout (`TUP` / `L2U2`); arcs and continuum bars are `ALL` / `L2U2L1U1`, with no bias in that mode (overscan subtraction only), as is usual for KCWI calibrations.

### How each pass was run

One directory per pass, `~/drp/kcwi/s3/<dataset>/<channel>/L<level>/`, plus `L3_replay/`. In each:

```bash
DATA=~/drp/kcwi/test_suite/suite/data
cp $DATA/$DS/kcwi.sky .
sed "s|^lev0/|$DATA/$DS/lev0/|" $DATA/$DS/files_flatlamp_$CH.txt > files.txt
reduce_kcwi --write_config
sed -i '' "s/^plot_level = .*/plot_level = $L/" kcwi.cfg
reduce_kcwi $FLAG -c kcwi.cfg -l files.txt -g                                   # L0, L1
tee answers_${DS}_${CH}_3.txt | reduce_kcwi $FLAG -c kcwi.cfg -l files.txt -g   # L3
reduce_kcwi $FLAG -c kcwi.cfg -l files.txt -g < ../L3/answers_${DS}_${CH}_3.txt # L3_replay
```

`FLAG` is `-b` or `-r`. Every `kcwi.cfg` differs from the shipped default only in `plot_level` (checked). The Bokeh server was running (`enable_bokeh = True`).

The level-3 answer files are in this directory (`<dataset>/answers_<dataset>_<channel>_3.txt`: 939, 844, 364 and 321 lines). As captured by `tee`, each had one more line: a blank line typed after the reduction had finished (while `tee` was still reading the terminal), which no prompt ever read. The first full regression run found it (2026-10-24: each L3 run recorded one prompt fewer than the file had lines, every recorded answer matching its line), and the line was removed from these files and from the regression data directory's copies, with the manifest checksums updated. The originals are still in each `~/drp/kcwi/s3/<dataset>/<channel>/L3/`. When capturing answers, end `tee` with Ctrl-D without pressing Enter first, or check the line count against the prompt record afterwards. Screenshots of the live-only plots are in each `L3/screenshots/` (60, 13, 9 and 4 files for ds01 blue, ds01 red, ds03 blue, ds03 red).

### Code used

- ds01 (both channels, 09:38–12:06): kcwidrp **1.3.3** from a plain `pip install .` of the repository (identical to `6270cba`). The 1.3.4 change only touches the BM ledge code, which BL and RH1 never reach, so these outputs are what 1.3.4 gives.
- ds03 (both channels, 14:15–15:26): kcwidrp **1.3.4**, editable install of `upgrade/pr0-tests` (HISTORY: `git version=upgrade-start-37-g14f7751`), whose reduction code equals `upgrade-start` (`af46243`) apart from the behaviour-neutral T.11 and B23 changes.

### Checks

- All 16 passes completed. No `Failed executing primitive`, `while invoking` or `EOFError` in any log.
- **Replays identical.** For each of the four setups, `L3_replay` is identical to `L3`: every FITS file in `redux/` (65–70 per run), every HDU, every header card apart from dates, HISTORY and file paths; the same PNG names; the same proctab. The comparison does detect real differences (L1 vs. L3 of ds03 blue differ in invsens and `icubes`, as the non-default answers should make them).
- **Levels 0 and 1 give identical cubes** (`icube*`, 15 files per setup, every HDU) on all four setups, as D16 assumes.
- `plots/` at level 0 holds exactly the two SolveArcs PNGs in every setup.
- **`Picking the highest peak` never appears.** The BM ledge code only runs for ds03 blue; there every pass found **0 peaks** and logged `No BM ledge peak found, skipping BM ledge correction` (the 1.3.4 behaviour, D23). The several-peaks case (B1) is therefore validated on synthetic data only.
- ERROR lines that do not stop the reduction, all in the red channel (5 per pass):
  - `CorrectDefects: Defect list not found for defect_L2U2L1U1_2x2.dat` (3 per pass): no defect list ships for red `L2U2L1U1` at 2×2 (the arc and contbar readout), so those frames are not defect-corrected. Known from T.4.
  - `MakeCube: No reference geometry, cannot make cube!` (2 per pass): red has three arcs; the geometry is solved from the first, and the second and third arcs find no geometry under their own names and skip their arc cubes, which nothing uses. Benign, but logged at ERROR.

### Timings (level 1, Apple M4 MacBook Pro `Mac16,1`, 10 cores (4P + 6E), 16 GB, macOS 26.5.1)

Seconds per primitive, from the DRPF event timestamps (1 s resolution), summed over the calls (number of calls in brackets):

| Setup | GenerateMaps | MakeCube | WavelengthCorrections | CorrectDar | MakeMasterFlat | Whole pass |
|---|---|---|---|---|---|---|
| ds01 blue | 4 (1) | 25 (5) | 5 (3) | 5 (3) | 15 (1) | 218 |
| ds01 red | 4 (3) | 25 (7) | 5 (3) | 3 (3) | 15 (1) | 193 |
| ds03 blue | 4 (1) | 27 (5) | 3 (3) | 4 (3) | 16 (1) | 183 |
| ds03 red | 4 (3) | 32 (7) | 4 (3) | 9 (3) | 15 (1) | 288 |

Level-3 passes took 14–39 minutes of wall time, set by the interactive answers; replays took 3–5 minutes.
