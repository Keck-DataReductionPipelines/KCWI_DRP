"""Config references (plan T.6): the config keys the code reads, the keys the
files define, and the two config files agree.

An AST scan of the package (``kcwidrp/``, not ``kcwidrp/tests/``) finds

* attribute-style reads ``<...>.config.instrument.<key>`` (59 distinct keys at
  upgrade-start; the number is asserted below so a change is noticed);
* ``getattr(<...>.config.instrument, '<key>', ...)`` and ``hasattr(...)``
  with a string literal (a non-literal name fails the test with file:line);
* in the scripts, reads through the local ``kcwi_config`` object:
  ``kcwi_config.<key>`` and ``kcwi_config.BLUE['<key>']`` / ``.RED[...]``.

Case: keckdrpframework's ``ConfigClass.__getattr__`` lower-cases the name
before looking it up in ``properties`` (and ConfigParser lower-cases the keys
of the file), so ``config.instrument.TAPERFRAC`` reads ``taperfrac``. All
comparisons below are made on lower-case names.

Any other use of ``config.instrument`` (for example binding it to a variable)
fails the scan test, because the key reads that follow could not be seen.
"""
import ast
import configparser
import pathlib
import warnings

import pytest
from keckdrpframework.config.framework_config import ConfigClass

import kcwidrp
from kcwidrp.core.kcwi_pkg_resources import get_resource_path

PKG_ROOT = pathlib.Path(kcwidrp.__file__).resolve().parent
CHANNELS = ("BLUE", "RED")

# --- allow-lists -------------------------------------------------------------

# Attribute-style reads of keys that are NOT in kcwi.cfg because a script sets
# them at run time on the instrument config object.
RUNTIME_SET_KEYS = {
    # reduce_kcwi.py / kcwi_rti.py: `kcwi_config.cwd = os.getcwd()`
    "cwd": "set by the scripts from os.getcwd() before the framework starts",
    # reduce_kcwi.py:526-527 from the -w and -W options (kcwi_rti.py same)
    "wait_for_event": "set by the scripts from the -w option",
    "continuous": "set by the scripts from the -W option",
}

# getattr(..., 'key', default) / hasattr(..., 'key') keys that may be absent
# from kcwi.cfg. Empty today: the keys read this way (the plot viewer options;
# the command-line override checks in reduce_kcwi.py and kcwi_rti.py) are all
# in kcwi.cfg. Add an entry, with a reason, only for a
# key whose absence from kcwi.cfg is intended.
GETATTR_MAY_BE_ABSENT = {}

# "Deprecated" allow-list: keys the code still reads for backward
# compatibility although kcwi.cfg no longer defines them (so they are exempt
# from the "must exist in kcwi.cfg" rule). Each entry must still be read by
# the code (a stale entry fails).
DEPRECATED_READ_ONLY = {
    "enable_bokeh": "bokeh plan 5.1: deprecated alias of enable_plot_viewer, "
                    "read by plot_config.viewer_requested() to honour old "
                    "config files and warn",
}

# Keys defined in kcwi.cfg that no code reads. Each needs a reason. The
# CRR_* entries are bug B25 (plan §13); remove them when it is fixed.
UNREAD_CFG_KEYS = {
    "crr_gain": "B25: RemoveCosmicRays.py passes gain=1.0 to detect_cosmics "
                "instead of reading CRR_GAIN",
    "crr_readnoise": "B25: RemoveCosmicRays.py computes the read noise from the "
                     "OSCNRN header keywords instead of reading CRR_READNOISE",
    "crr_pssl": "B25: not passed to detect_cosmics (astroscrappy default is used)",
    "crr_satlevel": "B25: not passed to detect_cosmics (astroscrappy default is "
                    "used)",
    "crr_niter": "B25: not passed to detect_cosmics (astroscrappy default is "
                 "used)",
}

# Keys in kcwi_koarti.cfg that are not in kcwi.cfg.
KOARTI_ONLY_KEYS = {
    "inter": "legacy key ('inter = 1' next to a commented-out "
             "'interactive = 1'); nothing reads it",
}

# Keys defined in both files with different values (compared after the
# ConfigClass conversion of the value). Both sides are pinned.
KOARTI_VALUE_DIFFERENCES = {
    # kcwi_koarti.cfg comment: "Set lower to enable CCR". The RTI file turns
    # cosmic ray rejection off in effect (no frame is that long).
    "crr_minexptime": (60.0, 100000.0),
}


# Keys read inside plot_config.viewer_requested(), which the launch scripts
# call with config.instrument (the scan cannot see inside it).
VIEWER_REQUESTED_KEYS = ("enable_plot_viewer", "enable_bokeh")


# --- the scan ------------------------------------------------------------------

def _package_files():
    return sorted(p for p in PKG_ROOT.rglob("*.py")
                  if "tests" not in p.relative_to(PKG_ROOT).parts)


def _rel(path):
    return str(path.relative_to(PKG_ROOT.parent))


def _is_instrument(node):
    """node is `<...>.config.instrument` (or `config.instrument`)."""
    if not (isinstance(node, ast.Attribute) and node.attr == "instrument"):
        return False
    inner = node.value
    return ((isinstance(inner, ast.Attribute) and inner.attr == "config")
            or (isinstance(inner, ast.Name) and inner.id == "config"))


class Scan:
    def __init__(self):
        self.attr_reads = {}      # lower key -> [file:line]
        self.attr_stores = {}     # lower key -> [file:line]
        self.getattr_keys = {}    # lower key -> [(file:line, has_default)]
        self.nonliteral = []      # file:line
        self.aliases = []         # file:line of other uses of .instrument
        self.script_reads = {}    # script -> {lower key: [file:line]}
        self.channel_reads = {}   # script -> {lower key: [file:line]}
        for path in _package_files():
            with warnings.catch_warnings():
                # unescaped "\*" in docstrings of the code under scan
                warnings.simplefilter("ignore", SyntaxWarning)
                tree = ast.parse(path.read_text(), filename=str(path))
            self._scan(path, tree)

    def _scan(self, path, tree):
        where = lambda n: f"{_rel(path)}:{n.lineno}"
        parent = {}
        for node in ast.walk(tree):
            for child in ast.iter_child_nodes(node):
                parent[child] = node
        is_script = "scripts" in path.relative_to(PKG_ROOT).parts
        for node in ast.walk(tree):
            if isinstance(node, ast.Attribute) and node.attr == "instrument":
                par = parent.get(node)
                if isinstance(par, ast.Attribute) and par.value is node:
                    store = isinstance(par.ctx, ast.Store)
                    target = self.attr_stores if store else self.attr_reads
                    target.setdefault(par.attr.lower(), []).append(where(par))
                elif (isinstance(par, ast.Call) and par.args
                        and par.args[0] is node
                        and isinstance(par.func, ast.Name)
                        and par.func.id in ("getattr", "hasattr")):
                    self._getattr(path, par)
                elif (isinstance(par, ast.Assign) and node in par.targets):
                    pass  # `framework.config.instrument = kcwi_config`
                elif (isinstance(par, ast.Call) and par.args
                        and par.args[0] is node
                        and isinstance(par.func, ast.Name)
                        and par.func.id == "viewer_requested"):
                    # plot_config.viewer_requested(config.instrument, ...)
                    # reads these two keys with getattr(..., None)
                    for key in VIEWER_REQUESTED_KEYS:
                        self.getattr_keys.setdefault(key, []).append(
                            (where(par), True))
                else:
                    self.aliases.append(where(node))
            if is_script:
                self._script_reads(path, node)

    def _getattr(self, path, call):
        loc = f"{_rel(path)}:{call.lineno}"
        name = call.args[1] if len(call.args) > 1 else None
        if isinstance(name, ast.Constant) and isinstance(name.value, str):
            has_default = (call.func.id == "getattr" and len(call.args) > 2)
            self.getattr_keys.setdefault(name.value.lower(), []).append(
                (loc, has_default or call.func.id == "hasattr"))
        else:
            self.nonliteral.append(loc)

    def _script_reads(self, path, node):
        script = path.stem
        loc = f"{_rel(path)}:{getattr(node, 'lineno', 0)}"
        # kcwi_config.<key>
        if (isinstance(node, ast.Attribute)
                and isinstance(node.value, ast.Name)
                and node.value.id == "kcwi_config"
                and node.attr not in CHANNELS
                and isinstance(node.ctx, ast.Load)):
            self.script_reads.setdefault(script, {}).setdefault(
                node.attr.lower(), []).append(loc)
        # kcwi_config.BLUE['key']
        if (isinstance(node, ast.Subscript)
                and isinstance(node.value, ast.Attribute)
                and isinstance(node.value.value, ast.Name)
                and node.value.value.id == "kcwi_config"
                and node.value.attr in CHANNELS
                and isinstance(node.slice, ast.Constant)):
            self.channel_reads.setdefault(script, {}).setdefault(
                str(node.slice.value).lower(), []).append(loc)


@pytest.fixture(scope="module")
def scan():
    return Scan()


# --- config files ------------------------------------------------------------

def _own_keys(relpath):
    path = get_resource_path("kcwidrp", relpath)
    parser = configparser.ConfigParser(default_section="KCWI")
    parser.read(path)
    return set(parser.defaults()), {s: set(parser._sections[s])
                                    for s in parser.sections()}


def _loaded(relpath):
    return ConfigClass(get_resource_path("kcwidrp", relpath),
                       default_section="KCWI")


@pytest.fixture(scope="module")
def kcwi_keys():
    return _own_keys("configs/kcwi.cfg")[0]


@pytest.fixture(scope="module")
def koarti_keys():
    return _own_keys("configs/kcwi_koarti.cfg")[0]


# --- tests ---------------------------------------------------------------------

def test_scan_finds_expected_reads(scan):
    # source: characterization @ dfbcf9d (plan T.6: 59 distinct keys)
    # (58 after bokeh plan 5.1: enable_bokeh is no longer read as an
    # attribute by the scripts; viewer_requested() reads it via getattr)
    assert len(scan.attr_reads) == 58, sorted(scan.attr_reads)
    # source: characterization @ upgrade-start, amended by bokeh 3.1 and 5.1:
    # reduce_kcwi.py and kcwi_rti.py (8 override checks each, same keys), the
    # two keys read through viewer_requested(), and StartPlotViewer.py's two
    # keys (StartBokeh.py's three keys are gone)
    assert sorted(scan.getattr_keys) == [
        "atoff", "enable_bokeh", "enable_plot_viewer", "fracmax", "linelist",
        "linethresh", "midfrac", "plot_viewer_open_browser",
        "plot_viewer_port", "skipsky", "taperfrac", "tukeyalpha"]


def test_no_aliases_of_config_instrument(scan):
    # source: specification: the scan only sees reads written as
    # <...>.config.instrument.<key>; any other use hides reads from it
    assert not scan.aliases, (
        "config.instrument used other than as `.key` or getattr/hasattr "
        f"first argument (keys read through it are invisible): "
        f"{scan.aliases}")


def test_getattr_names_are_literals(scan):
    # source: specification (plan T.6)
    assert not scan.nonliteral, (
        f"getattr/hasattr on config.instrument with a non-literal name: "
        f"{scan.nonliteral}")


def test_config_lookup_is_case_insensitive():
    # ConfigClass.__getattr__ lower-cases the name; the references test
    # relies on it (code reads TAPERFRAC, the file defines TAPERFRAC and the
    # parser stores 'taperfrac').
    config = _loaded("configs/kcwi.cfg")
    # source: keckdrpframework ConfigClass.__getattr__ (properties.get(key.lower()))
    assert config.TAPERFRAC == config.taperfrac == 0.2
    assert "taperfrac" in config.properties
    assert "TAPERFRAC" not in config.properties


def test_attribute_reads_exist_in_kcwi_cfg(scan, kcwi_keys):
    missing = {}
    for key, locs in scan.attr_reads.items():
        if key in kcwi_keys or key in RUNTIME_SET_KEYS:
            continue
        if key in DEPRECATED_READ_ONLY:
            continue
        missing[key] = locs
    # source: specification (plan T.6): each attribute-style key exists in
    # kcwi.cfg (or is on a commented allow-list)
    assert not missing, (
        "config.instrument.<key> read in the code but not defined in "
        f"kcwi.cfg: {missing}")


def test_runtime_set_allow_list_is_exact(scan, kcwi_keys):
    for key in RUNTIME_SET_KEYS:
        assert key in scan.attr_reads, f"stale allow-list entry {key}"
        assert key not in kcwi_keys, f"{key} is in kcwi.cfg: drop the entry"
    # each runtime key is really assigned by a script
    stored = set(scan.attr_stores)
    stored_by_script = {"wait_for_event", "continuous"} & stored
    assert stored_by_script == {"wait_for_event", "continuous"}


def test_getattr_keys_exist_or_are_listed(scan, kcwi_keys):
    missing = {}
    for key, uses in scan.getattr_keys.items():
        if key in kcwi_keys:
            continue
        if key in GETATTR_MAY_BE_ABSENT or key in DEPRECATED_READ_ONLY:
            assert all(has_default for _, has_default in uses), (
                f"{key}: listed as may-be-absent but read without a "
                f"default at {uses}")
            continue
        missing[key] = uses
    # source: specification (plan T.6)
    assert not missing, (
        f"getattr/hasattr keys absent from kcwi.cfg and not on an "
        f"allow-list: {missing}")


def test_allow_lists_are_not_stale(scan, kcwi_keys):
    read = set(scan.attr_reads) | set(scan.getattr_keys)
    for name, table in (("GETATTR_MAY_BE_ABSENT", GETATTR_MAY_BE_ABSENT),
                        ("DEPRECATED_READ_ONLY", DEPRECATED_READ_ONLY)):
        for key in table:
            assert key in read, f"{name}: {key} is no longer read"
            assert key not in kcwi_keys or name == "DEPRECATED_READ_ONLY"
    for key in UNREAD_CFG_KEYS:
        assert key in kcwi_keys, f"UNREAD_CFG_KEYS: {key} not in kcwi.cfg"
        assert key not in read, f"UNREAD_CFG_KEYS: {key} is read now"


def test_every_kcwi_cfg_key_is_read(scan):
    main_keys, sections = _own_keys("configs/kcwi.cfg")
    read = set(scan.attr_reads) | set(scan.getattr_keys)
    for reads in scan.script_reads.values():
        read |= set(reads)
    for reads in scan.channel_reads.values():
        read |= set(reads)
    unread = sorted(k for k in main_keys
                    if k not in read and k not in UNREAD_CFG_KEYS)
    # source: specification (plan T.6)
    assert not unread, (
        f"keys in kcwi.cfg that no code reads (remove them, or add them to "
        f"UNREAD_CFG_KEYS with a reason): {unread}")
    # the channel sections are read only by key name through kcwi_config.BLUE
    # / .RED in the scripts
    channel_read = set()
    for reads in scan.channel_reads.values():
        channel_read |= set(reads)
    for sec, keys in sections.items():
        extra = sorted(keys - channel_read)
        assert not extra, (f"[{sec}] keys no script reads through "
                           f"kcwi_config.{sec}[...]: {extra}")


def test_koarti_keys_exist_in_kcwi_cfg(kcwi_keys, koarti_keys):
    extra = sorted(k for k in koarti_keys
                   if k not in kcwi_keys and k not in KOARTI_ONLY_KEYS)
    # source: specification (plan T.6)
    assert not extra, (f"keys in kcwi_koarti.cfg but not in kcwi.cfg: "
                       f"{extra}")
    for key in KOARTI_ONLY_KEYS:
        assert key in koarti_keys and key not in kcwi_keys, (
            f"stale KOARTI_ONLY_KEYS entry {key}")


def test_koarti_value_differences_are_listed(kcwi_keys, koarti_keys):
    main = _loaded("configs/kcwi.cfg").properties
    rti = _loaded("configs/kcwi_koarti.cfg").properties
    differing = {k: (main[k], rti[k]) for k in sorted(koarti_keys & kcwi_keys)
                 if main[k] != rti[k]}
    # source: characterization @ upgrade-start (the one difference)
    assert differing == KOARTI_VALUE_DIFFERENCES, (
        f"values that differ between kcwi.cfg and kcwi_koarti.cfg: "
        f"{differing}; listed: {KOARTI_VALUE_DIFFERENCES}")


# Which config file each launch script loads when -c is not given.
SCRIPT_DEFAULT_CONFIG = {
    "reduce_kcwi": "configs/kcwi.cfg",
    "kcwi_rti": "configs/kcwi_koarti.cfg",
}


def _channel_params():
    params = []
    for script, relpath in SCRIPT_DEFAULT_CONFIG.items():
        marks = []
        if script == "kcwi_rti":
            marks.append(pytest.mark.xfail(
                strict=True,
                reason="B22: kcwi_koarti.cfg, the default config of "
                       "kcwi_rti, has no [BLUE]/[RED] sections, so "
                       "kcwi_config.BLUE['linethresh'] (kcwi_rti.py, with "
                       "-b or -r, which are required) is None[...] and "
                       "raises TypeError"))
        params.append(pytest.param(script, relpath, id=script, marks=marks))
    return params


@pytest.mark.parametrize("script, relpath", _channel_params())
def test_default_config_has_channel_keys_the_script_reads(
        scan, script, relpath):
    _, sections = _own_keys(relpath)
    wanted = set(scan.channel_reads[script])
    problems = []
    for sec in CHANNELS:
        missing = sorted(wanted - sections.get(sec, set()))
        if missing:
            problems.append(f"[{sec}] in {relpath} lacks {missing}")
    # source: specification: reduce_kcwi.py / kcwi_rti.py main() read these
    # keys from kcwi_config.BLUE / .RED for -b / -r, one of which is required
    assert not problems, "; ".join(problems)
