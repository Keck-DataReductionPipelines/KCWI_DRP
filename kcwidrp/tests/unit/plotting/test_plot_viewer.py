"""PlotViewer (bokeh plan 0.3 decision, 3.1 and "Interfaces fixed").

Written from the specification, not from the implementation.  A browser tab
is emulated by a WebSocket client that speaks mpl.js's protocol, as in
docs/prototype/viewer_0.3/check_proto.py.

Assumptions about what the spec leaves open (see the task report):
* the page is served at ``viewer.url`` and the WebSocket at ``<url>ws``
  (the prototype's layout, which is mpl's WebAgg embedding example);
* the first full frame a freshly connected tab receives is a PNG whose
  IHDR gives the figure's pixel size.

All tests use high ports from a random base so they do not collide with a
real viewer (5006) or with each other.
"""
import asyncio
import json
import logging
import os
import random
import socket
import struct
import subprocess
import sys
import textwrap
import time

import numpy as np
import pytest

LOW, HIGH = 25000, 55000

# Real servers and WebSocket round trips (~12 s for the file): kept out of
# the default run for the 60 s budget; run with `pytest -m viewer`.
pytestmark = pytest.mark.viewer


# ----------------------------------------------------------------- helpers
def _bindable(port):
    for host in ("127.0.0.1", "0.0.0.0"):
        s = socket.socket()
        try:
            s.bind((host, port))
        except OSError:
            return False
        finally:
            s.close()
    return True


def free_base(span=10):
    """A base port with `span` consecutive free ports above it."""
    for _ in range(500):
        base = random.randint(LOW, HIGH)
        if all(_bindable(base + i) for i in range(span)):
            return base
    raise RuntimeError("no free port block found")


@pytest.fixture
def base_port():
    return free_base()


@pytest.fixture
def holders():
    """Hold ports with listening sockets: hold(port, host='127.0.0.1')."""
    socks = []

    def hold(port, host="127.0.0.1"):
        s = socket.socket()
        s.bind((host, port))
        s.listen()
        socks.append(s)
        return s
    yield hold
    for s in socks:
        s.close()


@pytest.fixture
def make_viewer(base_port):
    from kcwidrp.core.plot_viewer import PlotViewer
    made = []

    def make(port=None, **kw):
        v = PlotViewer(port=base_port if port is None else port, **kw)
        made.append(v)
        return v
    yield make
    for v in made:
        try:
            v.stop()
        except Exception:
            pass


class Collect(logging.Handler):
    def __init__(self):
        super().__init__(logging.DEBUG)
        self.records = []

    def emit(self, record):
        self.records.append(record)

    def messages(self, minlevel=logging.DEBUG):
        return [r.getMessage() for r in self.records
                if r.levelno >= minlevel]


@pytest.fixture
def kcwi_log():
    handler = Collect()
    logger = logging.getLogger("KCWI")
    old = logger.level
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)
    yield handler
    logger.removeHandler(handler)
    logger.setLevel(old)


def line_plot(title="plot", w=640, h=400, y=None, n=200):
    from kcwidrp.core.kcwi_plotting import figure
    x = np.arange(n, dtype=float)
    y = np.sin(x / 15.0) if y is None else y
    p = figure(title=title, x_axis_label="x", y_axis_label="y",
               plot_width=w, plot_height=h)
    p.line(x, y)
    return p


def grid_plot(title="grid", ncols=2, nplots=4, w=300, h=200):
    from kcwidrp.core.kcwi_plotting import figure, gridplot
    plots = []
    for i in range(nplots):
        p = figure(title="%s %d" % (title, i), plot_width=w, plot_height=h)
        p.line(np.arange(50.0), np.cos(np.arange(50.0) / (5 + i)))
        plots.append(p)
    return gridplot(plots, ncols=ncols, plot_width=w, plot_height=h)


class Tab:
    """A fake browser tab speaking mpl.js's WebSocket protocol."""

    def __init__(self, viewer):
        self.url = viewer.url.replace("http://", "ws://").rstrip("/") + "/ws"
        self.events = []
        self.frames = []            # (w, h) of every binary PNG frame
        self.blobs = []
        self.conn = None
        OPEN_TABS.append(self)

    async def connect(self):
        from tornado.websocket import websocket_connect
        self.conn = await websocket_connect(self.url)
        for m in ({"type": "supports_binary", "value": True},
                  {"type": "send_image_mode"},
                  {"type": "set_device_pixel_ratio", "device_pixel_ratio": 1},
                  {"type": "refresh"}):
            await self.conn.write_message(json.dumps(m))

    QUIET = 0.15    # s of silence that ends a pump early

    async def pump(self, seconds):
        """Handle messages for up to ``seconds``.  Returns early once the
        exchange has settled: something arrived, no requested frame is
        outstanding, and nothing more came for QUIET seconds."""
        end = time.time() + seconds
        got_any = awaiting_frame = False
        while time.time() < end:
            wait = end - time.time()
            if got_any and not awaiting_frame:
                wait = min(wait, self.QUIET)
            try:
                msg = await asyncio.wait_for(self.conn.read_message(),
                                             max(0.01, wait))
            except asyncio.TimeoutError:
                if got_any and not awaiting_frame:
                    break
                continue
            if msg is None:
                break
            got_any = True
            if isinstance(msg, bytes):
                w, h = struct.unpack(">II", msg[16:24])
                self.frames.append((w, h))
                self.blobs.append(msg)
                awaiting_frame = False
                await self.conn.write_message(json.dumps({"type": "ack"}))
            else:
                ev = json.loads(msg)
                self.events.append(ev)
                if ev["type"] == "draw":
                    awaiting_frame = True
                    await self.conn.write_message(json.dumps({"type": "draw"}))

    async def send(self, **m):
        await self.conn.write_message(json.dumps(m))

    def close(self):
        # Closing needs a running loop; run() closes every tab it opened.
        pass

    def labels(self):
        return [e["label"] for e in self.events
                if e["type"] == "figure_label"]

    def last(self, typ):
        ev = [e for e in self.events if e["type"] == typ]
        return ev[-1] if ev else None


OPEN_TABS = []


def run(coro):
    loop = asyncio.new_event_loop()
    try:
        return loop.run_until_complete(coro)
    finally:
        for tab in OPEN_TABS:
            if tab.conn is not None:
                tab.conn.close()
        OPEN_TABS.clear()
        loop.run_until_complete(asyncio.sleep(0.05))
        loop.close()


def observe(viewer, seconds=1.5):
    """Open a tab on the viewer, listen, return the Tab."""
    tab = Tab(viewer)

    async def go():
        await tab.connect()
        await tab.pump(seconds)
    run(go())
    return tab


def pixels(blob):
    import io
    import matplotlib.image as mimg
    return mimg.imread(io.BytesIO(blob), format="png")


# ------------------------------------------------------------- start / ports
def test_start_returns_true_and_exposes_url_and_port(make_viewer, base_port):
    v = make_viewer()
    assert v.start() is True
    assert v.port == base_port                  # source: spec, first free port
    assert v.url == "http://127.0.0.1:%d/" % base_port
    assert v.dead is False


def test_default_port_is_5006():
    from kcwidrp.core.plot_viewer import PlotViewer
    # source: 0.3 -- "plot_viewer_port (default 5006)", open_browser default
    # False (S4-Q4).  Signature only: do not start, 5006 may be in use.
    import inspect
    assert inspect.signature(PlotViewer).parameters["port"].default == 5006
    assert inspect.signature(
        PlotViewer).parameters["open_browser"].default is False


def test_binds_loopback_only(make_viewer):
    import psutil                               # a package dependency
    v = make_viewer()
    assert v.start() is True
    proc = psutil.Process()
    conns = getattr(proc, "net_connections", None) or proc.connections
    listening = [c for c in conns(kind="inet")
                 if c.status == psutil.CONN_LISTEN and c.laddr.port == v.port]
    assert listening, "no listening socket found for the viewer port"
    # source: 0.3 decision -- "127.0.0.1 only"
    assert {c.laddr.ip for c in listening} == {"127.0.0.1"}


def test_serves_http_on_the_url(make_viewer):
    import urllib.request
    v = make_viewer()
    v.start()
    with urllib.request.urlopen(v.url, timeout=10) as r:
        assert r.status == 200
        assert b"<html" in r.read().lower()


def test_next_free_port_when_first_is_taken(make_viewer, base_port, holders):
    holders(base_port)
    v = make_viewer()
    assert v.start() is True
    assert v.port == base_port + 1


def test_skips_several_taken_ports(make_viewer, base_port, holders):
    for i in (0, 1, 2):
        holders(base_port + i)
    v = make_viewer()
    assert v.start() is True
    assert v.port == base_port + 3


def test_wildcard_listener_counts_as_taken(make_viewer, base_port, holders):
    # source: 0.3 decision -- a 127.0.0.1 bind can silently shadow a
    # 0.0.0.0 listener on macOS; such a port must still be skipped.
    holders(base_port, "0.0.0.0")
    v = make_viewer()
    assert v.start() is True
    assert v.port == base_port + 1


def test_last_port_of_range_is_used(make_viewer, base_port, holders):
    for i in range(9):
        holders(base_port + i)
    v = make_viewer()
    assert v.start() is True
    assert v.port == base_port + 9              # source: "up to port+9"


def test_all_ten_taken_start_false_warns(make_viewer, base_port, holders,
                                         kcwi_log):
    for i in range(10):
        # includes wildcard listeners, which must count as taken too
        holders(base_port + i, "0.0.0.0" if i % 3 == 0 else "127.0.0.1")
    v = make_viewer()
    assert v.start() is False
    warns = kcwi_log.messages(logging.WARNING)
    assert len(warns) >= 1


def test_show_after_failed_start_is_harmless(make_viewer, base_port, holders):
    for i in range(10):
        holders(base_port + i)
    v = make_viewer()
    assert v.start() is False
    for _ in range(5):
        v.show(line_plot())             # must not raise or block
        v.clear()
    v.stop()


def test_eleventh_port_is_never_used(make_viewer, base_port, holders):
    for i in range(10):
        holders(base_port + i)
    holders(base_port + 10)
    v = make_viewer()
    assert v.start() is False


def test_two_viewers_get_two_ports(make_viewer):
    a, b = make_viewer(), make_viewer()
    assert a.start() is True and b.start() is True
    assert a.port != b.port
    assert abs(a.port - b.port) == 1


def test_url_logged(make_viewer, kcwi_log):
    v = make_viewer()
    v.start()
    # source: 0.3/3.1 -- the URL is logged.  The INFO line and the printed
    # line come from StartPlotViewer (test_run_StartPlotViewer.py).
    assert any(v.url in r.getMessage() for r in kcwi_log.records)


# ---------------------------------------------------------- show never fails
def test_show_before_start_does_not_raise(make_viewer):
    v = make_viewer()
    v.show(line_plot())
    v.clear()


def test_show_garbage_does_not_raise(make_viewer):
    v = make_viewer()
    v.start()
    for junk in (None, object(), 42, "text", [1, 2, 3], {}):
        v.show(junk)
    v.show(line_plot())


def test_show_after_stop_does_not_raise(make_viewer):
    v = make_viewer()
    v.start()
    v.stop()
    v.show(line_plot())
    v.clear()
    v.stop()                                    # stop twice is harmless


def test_show_does_not_block_with_no_tab(make_viewer):
    v = make_viewer()
    v.start()
    p = line_plot(n=200000)                     # big enough to notice a render
    t0 = time.time()
    for _ in range(30):
        v.show(p)
    elapsed = time.time() - t0
    # source: 3.1 -- "show() returns without waiting for the render".  With
    # 30 submissions of 200k points an inline render would take far longer;
    # the generous bound only guards against blocking on rendering/IO.
    assert elapsed < 10.0


def test_show_does_not_block_while_a_tab_is_connected(make_viewer):
    v = make_viewer()
    v.start()
    tab = Tab(v)
    times = []

    async def go():
        await tab.connect()
        await tab.pump(0.5)
        p = line_plot(n=100000)
        for _ in range(10):
            t0 = time.time()
            v.show(p)
            times.append(time.time() - t0)
        await tab.pump(1.0)
    run(go())
    tab.close()
    # Tolerance: show() only deep-copies and queues; 1 s per call is a
    # ceiling that no queueing implementation approaches, not a benchmark.
    assert max(times) < 1.0


def test_no_browser_connected_no_exception(make_viewer):
    v = make_viewer()
    v.start()
    for i in range(20):
        v.show(line_plot("nobody %d" % i))
    time.sleep(0.3)
    assert v.dead is False


# ------------------------------------------------------- what a tab receives
def test_tab_receives_frame_of_plot_pixel_size(make_viewer):
    from kcwidrp.core.kcwi_plotting import pixel_size
    v = make_viewer()
    v.start()
    p = line_plot("size check", w=640, h=400)
    v.show(p)
    tab = observe(v)
    assert tab.frames, "no frame arrived"
    # source: 1.1 -- plot_width x plot_height are pixels at 100 dpi; the
    # PNG frame has exactly the plot's pixel size.
    assert tab.frames[0] == (640, 400)
    assert tab.frames[0] == tuple(pixel_size(p))


def test_figure_label_equals_plot_label(make_viewer):
    from kcwidrp.core.kcwi_plotting import plot_label
    v = make_viewer()
    v.start()
    p = line_plot("Label A")
    v.show(p)
    tab = observe(v)
    assert plot_label(p) == "Label A"           # source: spec, title text
    assert tab.labels()
    assert tab.labels()[-1] == "Label A"


def test_grid_label_is_first_figure_title_and_size_is_grid_size(make_viewer):
    v = make_viewer()
    v.start()
    g = grid_plot("cell", ncols=2, nplots=4, w=300, h=200)
    v.show(g)
    tab = observe(v)
    # source: 1.1 -- grid pixel size is ncols*plot_width x nrows*plot_height
    # (2*300 x 2*200); label = first figure's title (Interfaces, plot_label).
    assert tab.frames[0] == (600, 400)
    assert tab.labels()[-1] == "cell 0"


def test_each_submission_shown_at_its_own_size(make_viewer):
    v = make_viewer()
    v.start()
    cases = [(line_plot("single", w=640, h=400), (640, 400)),
             (grid_plot("g4", 2, 4, 300, 200), (600, 400)),
             (grid_plot("g2", 1, 2, 250, 300), (250, 600)),
             (line_plot("again", w=500, h=350), (500, 350))]
    tab = Tab(v)

    async def go():
        await tab.connect()
        await tab.pump(0.5)
        got = []
        for obj, size in cases:
            tab.frames.clear()
            v.show(obj)
            await tab.pump(1.2)
            got.append((size, tab.frames[0] if tab.frames else None))
        return got
    got = run(go())
    tab.close()
    for size, frame in got:
        assert frame == size
    # a resize event told the page about the newest size
    assert [int(round(s)) for s in tab.last("resize")["size"]] == [500, 350]


def test_latest_submission_wins(make_viewer):
    v = make_viewer()
    v.start()
    tab = Tab(v)

    async def go():
        await tab.connect()
        await tab.pump(0.5)
        for i in range(30):
            v.show(line_plot("burst %d" % i))
        await tab.pump(2.0)
    run(go())
    tab.close()
    assert tab.labels()[-1] == "burst 29"


def test_new_plot_clears_zoom_history(make_viewer):
    v = make_viewer()
    v.start()
    v.show(line_plot("zoom me", w=640, h=400))
    tab = Tab(v)

    async def go():
        await tab.connect()
        await tab.pump(1.0)
        # zoom to a box in the middle of the figure (inside any reasonable
        # axes area); browser coordinates, origin top-left
        await tab.send(type="toolbar_button", name="zoom")
        await tab.pump(0.3)
        await tab.send(type="button_press", x=250, y=150, button=0,
                       buttons=1, modifiers=[])
        await tab.send(type="motion_notify", x=400, y=250, button=0,
                       buttons=1, modifiers=[])
        await tab.send(type="button_release", x=400, y=250, button=0,
                       buttons=0, modifiers=[])
        await tab.pump(1.0)
        zoomed = tab.last("history_buttons")
        v.show(line_plot("fresh", w=640, h=400))
        await tab.pump(1.5)
        return zoomed, tab.last("history_buttons")
    zoomed, after = run(go())
    tab.close()
    assert zoomed is not None and zoomed["Back"] is True, \
        "precondition: the zoom gesture must create history"
    # source: 0.3 decision -- "toolbar history reset on a new plot"
    assert after["Back"] is False


def test_cursor_readout_reports_xy(make_viewer):
    v = make_viewer()
    v.start()
    v.show(line_plot("cursor", w=640, h=400))
    tab = Tab(v)

    async def go():
        await tab.connect()
        await tab.pump(1.0)
        await tab.send(type="motion_notify", x=320, y=200, button=0,
                       buttons=0, modifiers=[])
        await tab.pump(0.6)
    run(go())
    tab.close()
    msgs = [e["message"] for e in tab.events if e["type"] == "message"]
    # source: 0.3 decision -- "cursor readout ((x, y) = (...))"
    assert any("(x, y) = (" in m or "x=" in m for m in msgs), msgs


def test_reopened_tab_and_second_tab_show_current_plot(make_viewer):
    v = make_viewer()
    v.start()
    v.show(line_plot("first"))
    first = observe(v)
    assert "first" in first.labels()
    v.show(line_plot("current one"))
    tab_a, tab_b = Tab(v), Tab(v)

    async def go():
        await tab_a.connect()
        await tab_b.connect()
        await asyncio.gather(tab_a.pump(1.5), tab_b.pump(1.5))
    run(go())
    tab_a.close()
    tab_b.close()
    for tab in (tab_a, tab_b):
        assert tab.labels()[-1] == "current one"
        assert (640, 400) in tab.frames


def test_connect_before_first_show_then_show(make_viewer):
    v = make_viewer()
    v.start()
    tab = Tab(v)

    async def go():
        await tab.connect()
        await tab.pump(0.5)
        v.show(line_plot("late"))
        await tab.pump(1.5)
    run(go())
    tab.close()
    assert tab.labels()[-1] == "late"
    assert (640, 400) in tab.frames


def test_clear_then_new_plot_still_shown(make_viewer):
    v = make_viewer()
    v.start()
    v.show(line_plot("before clear"))
    tab = Tab(v)

    async def go():
        await tab.connect()
        await tab.pump(1.0)
        v.clear()
        await tab.pump(0.5)
        v.show(line_plot("after clear"))
        await tab.pump(1.5)
    run(go())
    tab.close()
    assert tab.labels()[-1] == "after clear"


# ------------------------------------------------- snapshot at show() time
def _first_frame_pixels(viewer, seconds=1.5):
    tab = observe(viewer, seconds)
    assert tab.blobs
    return pixels(tab.blobs[0])


def test_in_place_mutation_of_input_after_show_is_not_displayed(make_viewer,
                                                                base_port):
    v = make_viewer()
    v.start()
    data = np.zeros(200)
    p = line_plot("mutate", y=data)
    v.show(p)
    data[:] = 1e6                               # the worker reuses its array
    shown = _first_frame_pixels(v)

    ref = make_viewer(port=base_port + 5)
    ref.start()
    ref.show(line_plot("mutate", y=np.zeros(200)))
    expected = _first_frame_pixels(ref)

    # source: 0.3 -- "the displayed plot shows the values as they were at
    # submission".  A reference viewer is shown the same plot built from
    # never-mutated zeros; the rendered pixels must be identical (same
    # renderer, same input: exact equality).
    assert shown.shape == expected.shape
    assert np.array_equal(shown, expected)


def test_mutating_the_plot_object_after_show_is_not_displayed(make_viewer,
                                                              base_port):
    v = make_viewer()
    v.start()
    x = np.arange(200.0)
    p = line_plot("changed after", y=np.sin(x / 15.0))
    v.show(p)
    p.line(x, 5 * np.ones(200), color="red")    # extra glyph after hand-off
    shown = _first_frame_pixels(v)

    ref = make_viewer(port=base_port + 5)
    ref.start()
    ref.show(line_plot("changed after", y=np.sin(x / 15.0)))
    expected = _first_frame_pixels(ref)
    # source: 3.1 -- "show() deep-copies the recorded description on the
    # caller's thread": later changes to the object are not displayed.
    assert np.array_equal(shown, expected)


def test_sanity_different_data_gives_different_pixels(make_viewer, base_port):
    # Guards the two tests above against comparing blank frames.
    v = make_viewer()
    v.start()
    v.show(line_plot("s", y=np.zeros(200)))
    a = _first_frame_pixels(v)
    w = make_viewer(port=base_port + 5)
    w.start()
    w.show(line_plot("s", y=np.arange(200.0)))
    b = _first_frame_pixels(w)
    assert not np.array_equal(a, b)


# -------------------------------------------------------------- the browser
@pytest.fixture
def browser_calls(monkeypatch):
    import webbrowser
    calls = []
    monkeypatch.setattr(webbrowser, "open",
                        lambda url, *a, **k: calls.append(url) or False)
    return calls


def test_browser_not_opened_by_default(make_viewer, browser_calls):
    v = make_viewer()                           # open_browser False
    v.start()
    v.show(line_plot())
    time.sleep(0.3)
    assert browser_calls == []


def test_browser_not_opened_when_false(make_viewer, browser_calls):
    v = make_viewer(open_browser=False)
    v.start()
    for _ in range(3):
        v.show(line_plot())
    time.sleep(0.3)
    assert browser_calls == []


def test_browser_opened_at_most_once_when_true(make_viewer, browser_calls):
    v = make_viewer(open_browser=True)
    v.start()
    for i in range(10):
        v.show(line_plot("b %d" % i))
    time.sleep(0.5)
    assert len(browser_calls) <= 1              # source: 0.3, "at most once"
    if browser_calls:
        assert browser_calls[0].startswith(v.url.rstrip("/"))


def test_browser_failure_is_not_an_error(make_viewer, monkeypatch):
    import webbrowser

    def boom(*a, **k):
        raise webbrowser.Error("no runnable browser")
    monkeypatch.setattr(webbrowser, "open", boom)
    v = make_viewer(open_browser=True)
    assert v.start() is True                    # headless host: just go on
    v.show(line_plot())
    assert v.dead is False


# ------------------------------------------------------------ fork safety
_FORK_SCRIPT = textwrap.dedent('''
    import multiprocessing, os, sys, time
    from kcwidrp.core.plot_viewer import PlotViewer
    base = int(sys.argv[1])
    v = PlotViewer(port=base)
    assert v.start()
    method = sys.argv[2]
    if method == "pool":
        # as MakeCube: a default Pool() on a non-main thread of a process
        # that has a live viewer; children left running at exit
        pool = multiprocessing.get_context("fork").Pool(2)
        pool.map_async(time.sleep, [60, 60])
        time.sleep(0.5)
        pids = [p.pid for p in multiprocessing.active_children()]
    else:
        ctx = multiprocessing.get_context("fork")
        procs = [ctx.Process(target=time.sleep, args=(60,)) for _ in range(2)]
        for p in procs:
            p.start()
        pids = [p.pid for p in procs]
    print("PORT", v.port, flush=True)
    print("PIDS", *pids, flush=True)
    sys.stdin.readline()
    os._exit(0)             # as the framework: no atexit, no cleanup
''')


def _port_listening(port):
    s = socket.socket()
    s.settimeout(1.0)
    try:
        return s.connect_ex(("127.0.0.1", port)) == 0
    finally:
        s.close()


@pytest.mark.parametrize("method", ["process", "pool"])
def test_forked_children_do_not_hold_the_viewer_port(base_port, method,
                                                     tmp_path):
    env = dict(os.environ)
    env.pop("DISPLAY", None)
    errfile = tmp_path / "stderr.txt"       # a pipe would be held open by
    with open(errfile, "w") as err:         # the surviving children
        proc = subprocess.Popen(
            [sys.executable, "-c", _FORK_SCRIPT, str(base_port), method],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=err,
            text=True, env=env)
    pids = []
    try:
        port = None
        for _ in range(20):             # the viewer prints its URL first
            line = proc.stdout.readline()
            if not line:
                break
            if line.startswith("PORT"):
                port = int(line.split()[1])
            if line.startswith("PIDS"):
                pids = [int(p) for p in line.split()[1:]]
                break
        assert port is not None and pids, errfile.read_text()
        assert _port_listening(port)            # viewer is serving
        proc.stdin.write("\n")                  # parent leaves via os._exit
        proc.stdin.flush()
        proc.wait(30)
        for pid in pids:                        # children still alive
            os.kill(pid, 0)
        deadline = time.time() + 5
        while time.time() < deadline and _port_listening(port):
            time.sleep(0.1)
        # source: 0.3 decision -- "after the parent exits, no child holds
        # the port": the listening socket was closed in the forked children.
        assert not _port_listening(port)
        s = socket.socket()
        try:
            s.bind(("127.0.0.1", port))         # and it can be bound again
        finally:
            s.close()
    finally:
        for pid in pids:
            try:
                os.kill(pid, 9)
            except OSError:
                pass
        if proc.poll() is None:
            proc.kill()
        proc.wait()


def _square(i):
    return i * i


def test_fork_while_serving_a_tab_does_not_hang(make_viewer):
    import multiprocessing
    v = make_viewer()
    v.start()
    tab = Tab(v)

    async def go():
        await tab.connect()
        await tab.pump(0.5)
        ctx = multiprocessing.get_context("fork")
        for i in range(3):
            v.show(line_plot("pool %d" % i, n=20000))
            with ctx.Pool(2) as pool:
                assert pool.map(_square, range(8)) == [i * i
                                                      for i in range(8)]
            await tab.pump(0.2)
    run(go())
    tab.close()
    assert tab.labels()[-1] == "pool 2"
