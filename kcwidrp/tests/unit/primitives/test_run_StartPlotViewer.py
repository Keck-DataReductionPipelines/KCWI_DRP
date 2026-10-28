"""StartPlotViewer run through ``run_primitive`` (bokeh plan 3.1).

The primitive starts a PlotViewer on ``plot_viewer_port`` (or one of the
next nine), stores it as ``context.bokeh_session`` and its URL as
``context.plot_viewer_url``, and prints the URL on its own line; when no port
is free it warns and leaves ``bokeh_session`` None (never raises).
"""
import logging
import random
import socket

import pytest

from kcwidrp.primitives.StartPlotViewer import StartPlotViewer


@pytest.fixture
def base_port():
    return random.randint(30000, 40000) // 10 * 10


def _stop(ctx):
    v = getattr(ctx, "bokeh_session", None)
    if v is not None:
        v.stop()


def test_url_printed_on_one_line_and_stored(run_primitive, synthetic_ccddata,
                                            base_port, capsys, caplog):
    caplog.set_level(logging.INFO)
    run_primitive(StartPlotViewer, synthetic_ccddata(),
                  config_overrides={"plot_viewer_port": base_port})
    ctx = run_primitive.last.context
    try:
        url = ctx.plot_viewer_url
        # source: specification - plan 0.3 decision: http://127.0.0.1:<port>/
        # with the port in plot_viewer_port .. plot_viewer_port + 9
        assert url.startswith("http://127.0.0.1:")
        assert base_port <= int(url.rsplit(":", 1)[1].strip("/")) \
            <= base_port + 9
        assert ctx.bokeh_session.url == url
        # source: specification - plan 0.3: "printed to the terminal on its
        # own line" (a label around the URL on that line is allowed)
        out = capsys.readouterr().out
        assert len([ln for ln in out.splitlines() if url in ln]) == 1
        # source: specification - plan 0.3: "The chosen URL is logged at INFO"
        assert any(url in r.getMessage() and r.levelno == logging.INFO
                   for r in caplog.records)
    finally:
        _stop(ctx)


def test_no_free_port_leaves_session_none(run_primitive, synthetic_ccddata,
                                          base_port):
    held = []
    try:
        for port in range(base_port, base_port + 10):
            s = socket.socket()
            s.bind(("127.0.0.1", port))
            s.listen()
            held.append(s)
        run_primitive(StartPlotViewer, synthetic_ccddata(),
                      config_overrides={"plot_viewer_port": base_port})
        ctx = run_primitive.last.context
        # source: specification - plan 3.1: a failed start warns and the
        # reduction continues without a viewer (R10)
        assert ctx.bokeh_session is None
        assert ctx.plot_viewer_url is None
    finally:
        for s in held:
            s.close()
