"""The manual regression suite (T.7): one test per manifest entry and tier.

Run with ``pytest -m regression`` (quick profile) or
``pytest -m regression --regression-profile=full``; select with ``-k``, for
example ``-k ds01`` or ``-k structure``. The work is done by the
``entry_run`` fixture (conftest.py); the tiers are defined in
``kcwidrp.tests.regression.compare``. Skipped, with the reason
``KCWI_REGRESSION_DATA not set``, when the data directory is unset.
"""
import pytest


@pytest.mark.regression
def test_regression(entry_run, tier):
    result = entry_run.tier(tier)
    if result.status == "skipped":
        pytest.skip(result.summary)
    if result.status == "fail":
        lines = [result.summary] + ["  " + d for d in result.details[:20]]
        if len(result.details) > 20:
            lines.append("  ... %d more (see the report)" %
                         (len(result.details) - 20))
        pytest.fail("\n".join(lines), pytrace=False)
