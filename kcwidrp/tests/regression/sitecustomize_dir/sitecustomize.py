"""Input recorder for the regression runner (T.7, decisions D15 and D18).

The runner puts this directory first on the reduction subprocess's PYTHONPATH.
When ``KCWI_PROMPT_RECORD`` is set, this module

* wraps ``builtins.input`` so every call appends one JSON line to that file,
  ``{"prompt": ..., "answer": ...}`` or ``{"prompt": ..., "eof": true}`` (the
  EOFError is then re-raised), flushing each line; and
* sets ``astropy.utils.iers.conf.auto_download = False`` as soon as
  ``astropy.utils.iers`` is imported (a lazy hook: astropy is not imported
  here, so import order and start-up time are unchanged).

With ``KCWI_PROMPT_RECORD`` unset it does nothing. Any other ``sitecustomize``
that the interpreter would have found is executed afterwards.
"""
import builtins
import importlib.abc
import importlib.machinery
import importlib.util
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))


def _install_recorder(path):
    original = builtins.input

    def _write(obj):
        with open(path, "a", encoding="utf-8") as f:
            f.write(json.dumps(obj, ensure_ascii=False) + "\n")
            f.flush()

    def recording_input(prompt=""):
        try:
            answer = original(prompt)
        except EOFError:
            _write({"prompt": str(prompt), "eof": True})
            raise
        _write({"prompt": str(prompt), "answer": answer})
        return answer

    recording_input._kcwi_recorder = True
    builtins.input = recording_input


class _IersLoader(importlib.abc.Loader):
    def __init__(self, inner):
        self._inner = inner

    def create_module(self, spec):
        return self._inner.create_module(spec)

    def exec_module(self, module):
        self._inner.exec_module(module)
        try:
            module.conf.auto_download = False
        except Exception:                                # pragma: no cover
            pass


class _IersHook(importlib.abc.MetaPathFinder):
    target = "astropy.utils.iers"

    def find_spec(self, name, path, target=None):
        if name != self.target:
            return None
        spec = importlib.machinery.PathFinder.find_spec(name, path)
        if spec is None or spec.loader is None:
            return None
        spec.loader = _IersLoader(spec.loader)
        return spec


def _install_iers_hook():
    if "astropy.utils.iers" in sys.modules:
        try:
            sys.modules["astropy.utils.iers"].conf.auto_download = False
        except Exception:                                # pragma: no cover
            pass
        return
    sys.meta_path.insert(0, _IersHook())


def _chain():
    """Run the sitecustomize this one shadows, if there is one."""
    here = os.path.normcase(os.path.realpath(_HERE))
    rest = [p for p in sys.path
            if os.path.normcase(os.path.realpath(p or ".")) != here]
    spec = importlib.machinery.PathFinder.find_spec("sitecustomize", rest)
    if spec is None or spec.loader is None:
        return
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)


_record = os.environ.get("KCWI_PROMPT_RECORD")
if _record:
    _install_recorder(_record)
    _install_iers_hook()
_chain()
