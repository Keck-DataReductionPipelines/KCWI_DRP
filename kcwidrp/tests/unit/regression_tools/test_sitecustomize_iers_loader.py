"""The runner's sitecustomize must leave astropy.utils.iers as a normal import
would, except that conf.auto_download is False.

One snippet is run in a subprocess normally and again as the runner does it
(sitecustomize directory first on PYTHONPATH, KCWI_PROMPT_RECORD set).
"""
import json
import os
import subprocess
import sys

from kcwidrp.tests.regression import runner

RECORDER_DIR = str(runner.SITECUSTOMIZE_DIR)

SNIPPET = r'''
import hashlib, importlib.resources, inspect, json, pkgutil
import astropy.utils.iers as iers
import astropy.utils.data as data
ld = iers.__loader__
try:
    inspect.getsource(iers)
    src = True
except Exception:
    src = False
print(json.dumps({
    "loader": type(ld).__module__ + "." + type(ld).__name__,
    "get_resource_reader": hasattr(ld, "get_resource_reader"),
    "get_data": hasattr(ld, "get_data"),
    "is_package": hasattr(ld, "is_package"),
    "entries": sorted(p.name for p in
                      importlib.resources.files("astropy.utils.iers").iterdir()),
    "sha": hashlib.sha256(
        pkgutil.get_data("astropy.utils.iers", "iers.py")).hexdigest(),
    "getsource": src,
    "auto_download": iers.conf.auto_download,
    "data_loader": type(data.__loader__).__module__ + "."
                   + type(data.__loader__).__name__,
}))
'''


def run_snippet(tmp_path, runner_mode, record=True):
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    env.pop("KCWI_PROMPT_RECORD", None)
    if runner_mode:
        env["PYTHONPATH"] = RECORDER_DIR
        if record:
            env["KCWI_PROMPT_RECORD"] = str(tmp_path / "rec.jsonl")
    out = subprocess.run([sys.executable, "-c", SNIPPET], capture_output=True,
                         text=True, env=env, cwd=str(tmp_path))
    assert out.returncode == 0, out.stderr
    return json.loads(out.stdout.strip().splitlines()[-1])


def test_iers_module_identical_except_auto_download(tmp_path):
    # source: analytic -- the hook may change only conf.auto_download
    normal = run_snippet(tmp_path, False)
    hooked = run_snippet(tmp_path, True)
    assert normal["auto_download"] is True
    assert hooked["auto_download"] is False
    normal.pop("auto_download")
    hooked.pop("auto_download")
    assert hooked == normal


def test_iers_loader_keeps_normal_type_and_methods(tmp_path):
    # source: analytic -- spec: SourceFileLoader with the resource methods
    hooked = run_snippet(tmp_path, True)
    # (the defining module is _frozen_importlib_external on this Python)
    assert hooked["loader"].endswith(".SourceFileLoader")
    assert hooked["loader"] == run_snippet(tmp_path, False)["loader"]
    assert hooked["get_resource_reader"] is True
    assert hooked["get_data"] is True
    assert hooked["is_package"] is True


def test_iers_resources_getdata_and_getsource_work(tmp_path):
    # source: analytic -- spec: resources, pkgutil.get_data, inspect.getsource
    normal = run_snippet(tmp_path, False)
    hooked = run_snippet(tmp_path, True)
    assert "iers.py" in hooked["entries"]
    assert hooked["entries"] == normal["entries"]
    assert hooked["sha"] == normal["sha"]
    assert hooked["getsource"] is True


def test_other_modules_loader_untouched(tmp_path):
    # source: analytic -- spec: astropy.utils.data keeps its normal loader
    assert (run_snippet(tmp_path, True)["data_loader"]
            == run_snippet(tmp_path, False)["data_loader"])


def test_without_prompt_record_auto_download_is_default(tmp_path):
    # source: analytic -- spec: no KCWI_PROMPT_RECORD, the hook does nothing
    assert run_snippet(tmp_path, True, record=False)["auto_download"] is True
