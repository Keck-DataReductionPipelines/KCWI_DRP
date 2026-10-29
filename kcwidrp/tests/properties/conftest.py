"""Hypothesis set-up for the property tests (task T.5).

* Registers the ``kcwidrp-default`` and ``thorough`` profiles (``_profiles``)
  and loads ``kcwidrp-default`` unless ``--hypothesis-profile`` was given.
* Points hypothesis' storage directory (example database, unicode charmap
  cache, constants cache) at a fresh temporary directory, removed at the end
  of the session, so no ``.hypothesis/`` directory is created in the working
  tree (testing.md 4.1).  This must happen when the conftest is loaded: the
  constants cache is written before any fixture runs.
"""
import shutil
import tempfile

from hypothesis import settings
from hypothesis.configuration import set_hypothesis_home_dir

from kcwidrp.tests.properties import _profiles


def pytest_configure(config):
    # pytest_configure is a historic hook: it also runs for this conftest when
    # it is loaded after the hypothesis plugin's own pytest_configure (bare
    # ``pytest`` from the repository root).  Only load our default when the
    # user did not choose a profile.
    if not config.getoption("hypothesis_profile", default=None):
        settings.load_profile(_profiles.DEFAULT)
    home = tempfile.mkdtemp(prefix="kcwidrp-hypothesis-")
    set_hypothesis_home_dir(home)

    def _cleanup():
        set_hypothesis_home_dir(None)
        shutil.rmtree(home, ignore_errors=True)
    config.add_cleanup(_cleanup)
