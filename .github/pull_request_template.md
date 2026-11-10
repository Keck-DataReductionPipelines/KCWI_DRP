## Summary

<!-- What changes and why. Link the upgrade-plan task or issue. -->

## Fast suite

<!-- Paste the last line of `pytest`, e.g. "710 passed, 8 deselected, 60 xfailed in 32s".
     Plotting changed: also the last line of `pytest -m viewer`.
     Regression tooling changed: also the last line of `pytest -m regression_tools`. -->

## Regression suite

<!-- Reduction code changed: paste the summary table of report.md from `pytest -m regression`
     (quick profile; `--regression-profile=full` for plotting or prompt changes and releases).
     Otherwise write "not reduction code" and why. -->

## Documentation

<!-- User-visible behavior, config keys, scripts, outputs or docs/ changed: the docs/source/ pages
     updated, and the last line of `sphinx-build -W --keep-going docs/source docs/_build`
     (after `pip install -e ".[docs]"`). Otherwise write "no user-visible change". -->

## Blesses

<!-- For each bless (maintainer only): entry, reason ("<allowed reason>: <text>") and diff summary.
     Write "none" if there are none. -->

## Bugs fixed / xfails removed

<!-- Bug IDs fixed, with the xfail markers removed in the same commit. "none" if none. -->

## Who wrote the tests

<!-- Which model or person wrote the tests and which the code under test (docs/dev/testing.md, section 8).
     For each new test file, the deliberate breakage that made a test fail. -->

## Checklist

- [ ] `pytest` summary pasted above
- [ ] Regression summary pasted above, or "not reduction code" with the reason
- [ ] Every bless quoted with an allowed reason
- [ ] `xfail` markers removed with their fixes; no `xfail` without a bug or task ID
- [ ] No new `skip` without a reason (missing external resource only)
- [ ] Test authorship stated; breakage check done for each new test file
- [ ] Coverage lines for the touched modules pasted (`pytest --cov=kcwidrp --cov-report=term-missing`)
- [ ] Docs updated and built with `-W`, or "no user-visible change"
- [ ] `git status` shows no stray files
