"""Regression-suite hooks (T.1). The runner and fixtures arrive with T.7."""
import os

import pytest

SKIP_REASON = "KCWI_REGRESSION_DATA not set"


def pytest_collection_modifyitems(config, items):
    """Skip every test marked ``regression`` when the data directory is unset."""
    if os.environ.get("KCWI_REGRESSION_DATA"):
        return
    skip = pytest.mark.skip(reason=SKIP_REASON)
    for item in items:
        if item.get_closest_marker("regression") is not None:
            item.add_marker(skip)
