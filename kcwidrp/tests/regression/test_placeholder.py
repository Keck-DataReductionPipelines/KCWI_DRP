"""Placeholder so the regression skip can be validated; T.7 replaces it."""
import os
from pathlib import Path

import pytest


@pytest.mark.regression
def test_regression_data_directory_exists():
    # source: analytic -- the variable names the directory that holds the raw
    # data and blessed references (docs/testing.md section 10).
    data_dir = Path(os.environ["KCWI_REGRESSION_DATA"])
    assert data_dir.is_dir()
