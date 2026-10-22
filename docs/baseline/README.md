# Baseline environment and reference reductions

Tasks S.2 and S.3 of [drp_upgrade_plan.md](../drp_upgrade_plan.md). Only text files (`*.md`, `*.txt`, `*.json`) in this directory are committed (see `.gitignore`); FITS outputs, plots and logs from S.3 stay local.

## Baseline environment (S.2, 2026-10-22)

- Conda env `kcwidrp`, Python 3.12.13, macOS arm64 (`Darwin arm64`).
- `kcwidrp` installed editable (`pip install -e .`) from `upgrade-start` (`9be28f5`, `v1.3.2`).
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

## Reference reductions (S.3)

To be filled in during S.3: dataset IDs (KOA file names), observing dates, config used, command lines, machine (model, CPU, core count), per-primitive timings, whether `Picking the highest peak` appears in the level-1 logs, and anything unusual.
