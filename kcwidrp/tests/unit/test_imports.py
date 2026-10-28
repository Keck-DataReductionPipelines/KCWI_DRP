"""Every module under kcwidrp.{primitives,core,pipelines,scripts} imports.

One test case per module, discovered with pkgutil.walk_packages, so a new
module is picked up without editing this file (T.1). Replaces the generated
test_import_primitives.py and test_import_pipeline.py.
"""
import importlib
import pkgutil

import pytest

PACKAGES = ("kcwidrp.primitives", "kcwidrp.core", "kcwidrp.pipelines",
            "kcwidrp.scripts")

# Modules that fail to import in the baseline environment (python 3.12,
# bokeh 2.4.3, keckdrpframework 1.0.0 at upgrade-start), with the plan ID that
# will fix or remove the failure. strict=True: when the cause is fixed the
# test XPASSes, fails, and this entry must be deleted (docs/testing.md s.7).
KNOWN_IMPORT_FAILURES = {
    # ModuleNotFoundError: No module named 'keckdrpframework.utils.easyHTTP'
    # (absent from keckdrpframework 1.0.0).
    "kcwidrp.scripts.reduce_server": "B15: imports "
        "keckdrpframework.utils.easyHTTP, absent from keckdrpframework 1.0.0",
}


def discover_modules():
    names = []
    for package_name in PACKAGES:
        package = importlib.import_module(package_name)
        names.append(package_name)
        for info in pkgutil.walk_packages(package.__path__,
                                          package_name + "."):
            names.append(info.name)
    return sorted(set(names))


def _params():
    params = []
    for name in discover_modules():
        marks = []
        if name in KNOWN_IMPORT_FAILURES:
            marks.append(pytest.mark.xfail(
                strict=True, raises=ImportError,
                reason=KNOWN_IMPORT_FAILURES[name]))
        params.append(pytest.param(name, id=name, marks=marks))
    return params


@pytest.mark.parametrize("module_name", _params())
def test_module_imports(module_name):
    importlib.import_module(module_name)
