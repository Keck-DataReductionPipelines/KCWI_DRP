"""Config keys are documented both ways (task D.4).

* every key of every section of ``kcwidrp/configs/kcwi.cfg`` appears in
  ``docs/source/configuration.rst``, as ``key`` (double backticks) or as a
  ``key = ...`` line inside a ``.. code-block::`` block;
* every ``name = value`` line inside a ``.. code-block:: python`` block of
  that page names a key that exists in ``kcwi.cfg`` (any section).

Case of keys is kept (``TAPERFRAC``), so the comparison is case-sensitive.
"""
import configparser
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
CFG = REPO / "kcwidrp" / "configs" / "kcwi.cfg"
RST = REPO / "docs" / "source" / "configuration.rst"

DIRECTIVE = re.compile(r"^(\s*)\.\.\s+code-block::\s*(\S*)\s*$")
ASSIGN = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=")


def _cfg_keys():
    parser = configparser.ConfigParser()
    parser.optionxform = str
    with open(CFG) as fh:
        parser.read_file(fh)
    return {k for s in parser.sections() for k in parser[s]}


def _block_keys(language=None):
    """Keys from ``name = value`` lines in code-block directives of the page
    (only those of ``language`` when given)."""
    keys = set()
    lines = RST.read_text().splitlines()
    i = 0
    while i < len(lines):
        m = DIRECTIVE.match(lines[i])
        i += 1
        if not m:
            continue
        base = len(m.group(1))
        lang = m.group(2)
        while i < len(lines):
            line = lines[i]
            if line.strip() and len(line) - len(line.lstrip()) <= base:
                break
            i += 1
            if language is not None and lang != language:
                continue
            if line.lstrip().startswith("#"):
                continue
            a = ASSIGN.match(line)
            if a:
                keys.add(a.group(1))
    return keys


def test_every_cfg_key_is_documented():
    # source: published - keys of kcwidrp/configs/kcwi.cfg must appear in
    # docs/source/configuration.rst
    text = RST.read_text()
    in_blocks = _block_keys()
    missing = sorted(k for k in _cfg_keys()
                     if "``%s``" % k not in text and k not in in_blocks)
    assert not missing, "\n".join(
        "kcwi.cfg key %s is not documented in configuration.rst" % k
        for k in missing)


def test_every_documented_key_exists():
    # source: published - python code-blocks of docs/source/configuration.rst
    # list keys of kcwidrp/configs/kcwi.cfg
    extra = sorted(_block_keys("python") - _cfg_keys())
    assert not extra, "\n".join(
        "configuration.rst documents %s, which is not a kcwi.cfg key" % k
        for k in extra)
