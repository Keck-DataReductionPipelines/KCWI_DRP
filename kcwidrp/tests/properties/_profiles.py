"""Hypothesis settings profiles for the property tests (task T.5).

Rules: docs/testing.md section 9.  Registered here and loaded by
``kcwidrp/tests/properties/conftest.py``.

Profiles
--------
kcwidrp-default
    Loaded automatically.  ``derandomize=True`` (the same examples on every
    run, so a pass or a failure is reproducible), ``max_examples=25``,
    no example database (nothing is written to ``.hypothesis/``), no
    deadline (a slow laptop must not turn timing into a failure: the suite
    time is watched by the 60 s fast-suite rule instead).
thorough
    For manual runs:
    ``pytest kcwidrp/tests/properties --hypothesis-profile=thorough``.
    Random (not derandomized), ``max_examples=1000``, no database.

The profile names must be registered before hypothesis' pytest plugin reads
``--hypothesis-profile`` (in ``pytest_configure``).  That works when the
properties directory (or a file in it) is named on the command line, which
makes this directory's conftest.py an *initial* conftest.  A bare
``pytest --hypothesis-profile=thorough`` from the repository root fails with
"profile not registered" unless the root conftest imports this module.

Tests that run a whole primitive (several runs per example) use
``scaled(fraction)`` so that both profiles stay proportionate.
"""
from hypothesis import HealthCheck, settings

DEFAULT = "kcwidrp-default"
THOROUGH = "thorough"

_COMMON = dict(
    database=None,          # no example database: nothing written to the repo
    deadline=None,          # timing is not part of any invariant
    print_blob=False,
    # The detector-arithmetic and proctab properties use the factory fixtures
    # (synthetic_ccddata, run_primitive) inside @given.  The factories build
    # a fresh object on every call, so reusing the function-scoped fixture
    # across examples shares no state between examples.
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)

settings.register_profile(DEFAULT, derandomize=True, max_examples=25,
                          **_COMMON)
settings.register_profile(THOROUGH, derandomize=False, max_examples=1000,
                          **_COMMON)


def scaled(fraction):
    """``settings`` with ``fraction`` of the active profile's max_examples.

    Evaluated when the test module is imported, i.e. after the profile has
    been chosen (pytest_configure runs before collection).
    """
    n = max(3, int(round(settings.default.max_examples * fraction)))
    return settings(max_examples=n)
