"""plot_config.viewer_requested: the enable switch of bokeh plan 5.1.

Specification (docs/bokeh_port_plan.md, "The enable switch" and the
"Interfaces fixed for parallel implementation" subsection):

* ``enable_plot_viewer`` present: it decides; ``enable_bokeh`` is then ignored.
* Only ``enable_bokeh`` present: honoured.
* Neither present: viewer off.
* ``enable_bokeh`` present (in any combination) -> a deprecation message is
  logged, and at most once per process.  ``enable_bokeh`` absent -> no message.
* The module imports no Matplotlib.

Configs are real ``keckdrpframework`` ``ConfigClass`` objects, as in a
reduction (``framework.config.instrument``).  Read from its source: an absent
key reads as ``None`` (``ConfigClass.__getattr__``), ``True``/``False`` in a
cfg file are parsed to bools (``_getType``), and a quoted value such as
``'False'`` stays a string.

Each test loads a private copy of the module (``fresh_plot_config``) so the
"once per process" state never depends on test order.
"""
import importlib.util
import itertools
import logging
import subprocess
import sys
import textwrap

import pytest

MODULE = "kcwidrp.core.plot_config"
ABSENT = "absent"


def make_instrument(tmp_path, **values):
    """A real ConfigClass; keys with value ABSENT are left out of the file."""
    from keckdrpframework.config.framework_config import ConfigClass
    lines = ["[KCWI]", "plot_level = 1"]
    strings = {}
    for key, val in values.items():
        if val is ABSENT:
            continue
        if isinstance(val, str):
            strings[key] = val          # delivered as a str, set below
        else:
            lines.append("%s = %s" % (key, val))
    path = tmp_path / "kcwi_test.cfg"
    path.write_text("\n".join(lines) + "\n")
    config = ConfigClass(str(path), default_section="KCWI")
    for key, val in strings.items():
        config.properties[key] = val
    return config


class Capture(logging.Handler):
    def __init__(self):
        super().__init__(logging.DEBUG)
        self.records = []

    def emit(self, record):
        self.records.append(record)

    def deprecations(self):
        return [r for r in self.records
                if "enable_bokeh" in r.getMessage()]


@pytest.fixture
def fresh_plot_config():
    """Load a private copy of the module, without touching sys.modules."""
    spec = importlib.util.find_spec(MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def logger_capture():
    logger = logging.getLogger("plot_config_test_%d" % id(object()))
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    handler = Capture()
    logger.addHandler(handler)
    return logger, handler


# The key is (enable_plot_viewer, enable_bokeh); values True, False or ABSENT.
def expected_viewer(epv, eb):
    if epv is not ABSENT:
        return epv          # source: 5.1 "present: it decides"
    if eb is not ABSENT:
        return eb           # source: 5.1 "only enable_bokeh: honoured"
    return False            # source: 5.1 "neither present: viewer off"


CELLS = list(itertools.product([ABSENT, True, False], repeat=2))
CELL_IDS = ["epv_%s-eb_%s" % (a, b) for a, b in CELLS]


@pytest.mark.parametrize("epv,eb", CELLS, ids=CELL_IDS)
def test_enable_table(fresh_plot_config, logger_capture, tmp_path, epv, eb):
    logger, handler = logger_capture
    cfg = make_instrument(tmp_path, enable_plot_viewer=epv, enable_bokeh=eb)
    result = fresh_plot_config.viewer_requested(cfg, logger)
    assert result is expected_viewer(epv, eb)
    # source: 5.1 -- "the deprecation message appears exactly when
    # enable_bokeh is present"
    assert len(handler.deprecations()) == (0 if eb is ABSENT else 1)


@pytest.mark.parametrize("epv,eb", CELLS, ids=CELL_IDS)
def test_enable_table_string_values(fresh_plot_config, logger_capture,
                                    tmp_path, epv, eb):
    # A quoted value reaches the code as the string "True" / "False".  The
    # string "False" must not be treated as truthy.
    def as_str(v):
        return v if v is ABSENT else str(v)
    logger, handler = logger_capture
    cfg = make_instrument(tmp_path, enable_plot_viewer=as_str(epv),
                          enable_bokeh=as_str(eb))
    assert isinstance(cfg.enable_plot_viewer, (str, type(None)))
    result = fresh_plot_config.viewer_requested(cfg, logger)
    assert result is expected_viewer(epv, eb)
    assert len(handler.deprecations()) == (0 if eb is ABSENT else 1)


def test_neither_key_is_off_and_silent(fresh_plot_config, logger_capture,
                                       tmp_path):
    logger, handler = logger_capture
    cfg = make_instrument(tmp_path)
    assert fresh_plot_config.viewer_requested(cfg, logger) is False
    assert handler.records == []


def test_enable_bokeh_ignored_when_new_key_present(fresh_plot_config,
                                                   logger_capture, tmp_path):
    logger, handler = logger_capture
    cfg = make_instrument(tmp_path, enable_plot_viewer=False,
                          enable_bokeh=True)
    assert fresh_plot_config.viewer_requested(cfg, logger) is False
    cfg = make_instrument(tmp_path, enable_plot_viewer=True,
                          enable_bokeh=False)
    assert fresh_plot_config.viewer_requested(cfg, logger) is True


def test_deprecation_logged_once_in_process(fresh_plot_config, logger_capture,
                                            tmp_path):
    logger, handler = logger_capture
    cfg = make_instrument(tmp_path, enable_bokeh=True)
    for _ in range(5):
        assert fresh_plot_config.viewer_requested(cfg, logger) is True
    assert len(handler.deprecations()) == 1
    # a call without the old key in between does not matter either
    clean = make_instrument(tmp_path, enable_plot_viewer=True)
    fresh_plot_config.viewer_requested(clean, logger)
    fresh_plot_config.viewer_requested(cfg, logger)
    assert len(handler.deprecations()) == 1


def test_deprecation_message_names_replacement(fresh_plot_config,
                                               logger_capture, tmp_path):
    logger, handler = logger_capture
    cfg = make_instrument(tmp_path, enable_bokeh=True)
    fresh_plot_config.viewer_requested(cfg, logger)
    text = handler.deprecations()[0].getMessage()
    # source: 5.1 -- "a deprecation message" about enable_bokeh, which
    # enable_plot_viewer replaces.
    assert "deprecat" in text.lower()
    assert "enable_plot_viewer" in text


def test_logger_defaults_to_kcwi_logger(fresh_plot_config, tmp_path):
    handler = Capture()
    kcwi = logging.getLogger("KCWI")
    kcwi.addHandler(handler)
    old_level = kcwi.level
    kcwi.setLevel(logging.DEBUG)
    try:
        cfg = make_instrument(tmp_path, enable_bokeh=True)
        fresh_plot_config.viewer_requested(cfg)     # logger=None
        fresh_plot_config.viewer_requested(cfg)
    finally:
        kcwi.removeHandler(handler)
        kcwi.setLevel(old_level)
    # source: 3.2 interface -- logger=None means the package logger 'KCWI'
    # (assumption: the plan says "logger=None" without naming the fallback);
    # and still only once.
    assert len(handler.deprecations()) == 1


_ONCE_SCRIPT = textwrap.dedent('''
    import logging, sys, tempfile, os
    from keckdrpframework.config.framework_config import ConfigClass
    from kcwidrp.core.plot_config import viewer_requested

    class H(logging.Handler):
        n = 0
        def emit(self, record):
            if "enable_bokeh" in record.getMessage():
                H.n += 1

    logger = logging.getLogger(sys.argv[2])
    logger.setLevel(logging.DEBUG)
    logger.addHandler(H())
    d = tempfile.mkdtemp()
    def cfg(text):
        p = os.path.join(d, "c.cfg")
        open(p, "w").write("[KCWI]\\n" + text)
        return ConfigClass(p, default_section="KCWI")
    clean = cfg("enable_plot_viewer = True\\n")
    old = cfg("enable_bokeh = True\\n")
    use = None if sys.argv[1] == "none" else logger
    viewer_requested(clean, use)   # no old key: must not use up the "once"
    assert H.n == 0
    for _ in range(4):
        viewer_requested(old, use)
    print("COUNT", H.n)
''')


@pytest.mark.parametrize("mode,logger_name",
                         [("explicit", "plot_config_once_probe"),
                          ("none", "KCWI")])
def test_deprecation_once_per_process_fresh_interpreter(mode, logger_name):
    # A fresh interpreter has no state from other tests (order independence).
    proc = subprocess.run([sys.executable, "-c", _ONCE_SCRIPT, mode,
                           logger_name],
                          capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, proc.stderr
    assert "COUNT 1" in proc.stdout


def test_each_process_logs_once():
    # "once per process", not once per machine: two interpreters, one each.
    for _ in range(2):
        proc = subprocess.run([sys.executable, "-c", _ONCE_SCRIPT, "explicit",
                               "plot_config_once_probe2"],
                              capture_output=True, text=True, timeout=120)
        assert "COUNT 1" in proc.stdout, proc.stderr


def test_module_does_not_import_matplotlib():
    code = ("import sys, kcwidrp.core.plot_config; "
            "print([m for m in ('matplotlib', 'matplotlib.pyplot', 'bokeh', "
            "'selenium') if m in sys.modules])")
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True,
                          text=True, timeout=120)
    assert proc.returncode == 0, proc.stderr
    # source: 3.2/5.1 interface -- "No Matplotlib import"; the old stack must
    # not come in through this module either.
    assert proc.stdout.strip() == "[]"


def test_returns_a_real_bool(fresh_plot_config, tmp_path):
    for epv in (True, False):
        cfg = make_instrument(tmp_path, enable_plot_viewer=epv)
        assert type(fresh_plot_config.viewer_requested(cfg)) is bool
