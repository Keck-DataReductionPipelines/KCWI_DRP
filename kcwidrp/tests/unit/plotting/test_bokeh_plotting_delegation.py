"""bokeh_plotting delegates to the plot viewer (bokeh plan 3.1).

Specification: ``bokeh_plot(plot, session)`` -> ``session.show(plot)``;
``bokeh_clear(session)`` -> ``session.clear()``; both are no-ops for ``None``
and never raise.  ``init_plot_session`` and ``check_running_process`` are kept.
"""
import logging
import subprocess
import sys
import time
import uuid

import pytest


class FakeViewer:
    def __init__(self, show_exc=None, clear_exc=None):
        self.shown = []
        self.cleared = 0
        self.show_exc = show_exc
        self.clear_exc = clear_exc

    def show(self, obj):
        self.shown.append(obj)
        if self.show_exc:
            raise self.show_exc

    def clear(self):
        self.cleared += 1
        if self.clear_exc:
            raise self.clear_exc


class Collect(logging.Handler):
    def __init__(self):
        super().__init__(logging.DEBUG)
        self.records = []

    def emit(self, record):
        self.records.append(record)


@pytest.fixture
def warnings_seen():
    """WARNING+ records from 'KCWI' and from the root logger."""
    handler = Collect()
    loggers = [logging.getLogger("KCWI"), logging.getLogger()]
    old = [lg.level for lg in loggers]
    for lg in loggers:
        lg.addHandler(handler)
        lg.setLevel(logging.DEBUG)
    yield lambda: [r for r in handler.records if r.levelno >= logging.WARNING]
    for lg, lvl in zip(loggers, old):
        lg.removeHandler(handler)
        lg.setLevel(lvl)


@pytest.fixture
def a_plot():
    from kcwidrp.core.kcwi_plotting import figure
    p = figure(title="delegation", x_axis_label="x", y_axis_label="y",
               plot_width=300, plot_height=200)
    p.line([0, 1, 2], [0, 1, 4])
    return p


def test_plot_with_no_session_returns_none(a_plot):
    from kcwidrp.core.bokeh_plotting import bokeh_plot
    assert bokeh_plot(a_plot, None) is None


def test_clear_with_no_session_returns_none():
    from kcwidrp.core.bokeh_plotting import bokeh_clear
    assert bokeh_clear(None) is None


def test_plot_delegates_to_show_once(a_plot):
    from kcwidrp.core.bokeh_plotting import bokeh_plot
    viewer = FakeViewer()
    assert bokeh_plot(a_plot, viewer) is None
    assert len(viewer.shown) == 1
    assert viewer.shown[0] is a_plot
    assert viewer.cleared == 0


def test_clear_delegates_to_clear_once():
    from kcwidrp.core.bokeh_plotting import bokeh_clear
    viewer = FakeViewer()
    assert bokeh_clear(viewer) is None
    assert viewer.cleared == 1
    assert viewer.shown == []


def test_plot_accepts_a_grid(a_plot):
    from kcwidrp.core.bokeh_plotting import bokeh_plot
    from kcwidrp.core.kcwi_plotting import gridplot
    grid = gridplot([a_plot, a_plot], ncols=2, plot_width=300,
                    plot_height=200)
    viewer = FakeViewer()
    bokeh_plot(grid, viewer)
    assert viewer.shown == [grid]


def test_failing_show_is_swallowed_and_warned(a_plot, warnings_seen):
    from kcwidrp.core.bokeh_plotting import bokeh_plot
    viewer = FakeViewer(show_exc=RuntimeError("viewer exploded"))
    bokeh_plot(a_plot, viewer)                  # must not raise
    assert len(viewer.shown) == 1
    assert len(warnings_seen()) >= 1


def test_failing_clear_is_swallowed_and_warned(warnings_seen):
    from kcwidrp.core.bokeh_plotting import bokeh_clear
    viewer = FakeViewer(clear_exc=OSError("socket gone"))
    bokeh_clear(viewer)                         # must not raise
    assert viewer.cleared == 1
    assert len(warnings_seen()) >= 1


def test_base_style_exceptions_from_show_do_not_escape(a_plot):
    # "never raise": any ordinary Exception subclass, including odd ones.
    from kcwidrp.core.bokeh_plotting import bokeh_plot
    for exc in (ValueError("v"), KeyError("k"), TypeError("t"),
                AttributeError("a")):
        bokeh_plot(a_plot, FakeViewer(show_exc=exc))


def test_init_plot_session_sets_none_when_missing():
    from types import SimpleNamespace
    from kcwidrp.core.bokeh_plotting import init_plot_session
    context = SimpleNamespace()
    init_plot_session(context)
    # source: existing docstring (B9) -- "None if Bokeh is not running"
    assert context.bokeh_session is None


def test_init_plot_session_keeps_existing_session():
    from types import SimpleNamespace
    from kcwidrp.core.bokeh_plotting import init_plot_session
    viewer = FakeViewer()
    context = SimpleNamespace(bokeh_session=viewer)
    init_plot_session(context)
    assert context.bokeh_session is viewer


def test_check_running_process_finds_and_misses():
    from kcwidrp.core.bokeh_plotting import check_running_process
    token = "kcwi_probe_%s" % uuid.uuid4().hex
    proc = subprocess.Popen([sys.executable, "-c",
                             "import time; time.sleep(60)", token])
    try:
        deadline = time.time() + 10
        found = False
        while time.time() < deadline and not found:
            found = check_running_process(token)
            time.sleep(0.05)
        assert found is True
    finally:
        proc.kill()
        proc.wait()
    # source: analytic -- a random token that is in no command line
    assert check_running_process("kcwi_absent_%s" % uuid.uuid4().hex) is False
