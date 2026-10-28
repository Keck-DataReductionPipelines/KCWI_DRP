"""Viewer: optional vertical line at the cursor.

Written from the specification.  The page has a checkbox (id ``vline``);
after every draw each tab receives ``{"type": "axes_extents", "extents":
[[left, top, right, bottom], ...]}`` in CSS pixels from the figure's top-left
corner, one entry per distinct axes area.
"""
import asyncio
import json
import random
import re
import socket
import time
import urllib.request

import numpy as np
import pytest

# These tests drive a live viewer; kept out of the default run (docs/testing.md
# section 1).  Run with `pytest -m viewer`.
pytestmark = pytest.mark.viewer


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


@pytest.fixture
def viewer():
    from kcwidrp.core.plot_viewer import PlotViewer
    for _ in range(500):
        port = random.randint(41000, 49000)
        if _bindable(port):
            break
    v = PlotViewer(port=port, open_browser=False)
    assert v.start() is True
    yield v
    v.stop()


class Tab:
    """A fake browser tab speaking mpl.js's WebSocket protocol (pattern of
    test_plot_viewer.py)."""

    def __init__(self, viewer):
        self.url = viewer.url.replace("http://", "ws://").rstrip("/") + "/ws"
        self.events = []
        self.conn = None

    async def connect(self):
        from tornado.websocket import websocket_connect
        self.conn = await websocket_connect(self.url)
        for m in ({"type": "supports_binary", "value": True},
                  {"type": "send_image_mode"},
                  {"type": "set_device_pixel_ratio", "device_pixel_ratio": 1},
                  {"type": "refresh"}):
            await self.conn.write_message(json.dumps(m))

    async def pump(self, seconds, quiet=0.3):
        end = time.time() + seconds
        got_any = awaiting_frame = False
        while time.time() < end:
            wait = end - time.time()
            if got_any and not awaiting_frame:
                wait = min(wait, quiet)
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
                awaiting_frame = False
                await self.conn.write_message(json.dumps({"type": "ack"}))
            else:
                ev = json.loads(msg)
                self.events.append(ev)
                if ev["type"] == "draw":
                    awaiting_frame = True
                    await self.conn.write_message(json.dumps({"type": "draw"}))

    def extents(self):
        return [e["extents"] for e in self.events
                if e["type"] == "axes_extents"]


def observe(viewer, obj, seconds=3.0, then=None):
    """Show ``obj``, open a tab, listen; optionally show ``then`` and listen
    again.  Returns the Tab."""
    viewer.show(obj)
    tab = Tab(viewer)
    loop = asyncio.new_event_loop()
    try:
        async def go():
            await tab.connect()
            await tab.pump(seconds)
            if then is not None:
                viewer.show(then)
                await tab.pump(seconds)
        loop.run_until_complete(go())
        tab.conn.close()
        loop.run_until_complete(asyncio.sleep(0.05))
    finally:
        loop.close()
    return tab


def single_plot(w=640, h=400):
    from kcwidrp.core.kcwi_plotting import figure
    x = np.arange(100.0)
    p = figure(title="one", x_axis_label="x", y_axis_label="y",
               plot_width=w, plot_height=h)
    p.line(x, np.sin(x / 9.0))
    return p


def grid_2x2():
    from kcwidrp.core.kcwi_plotting import figure, gridplot
    figs = []
    for i in range(4):
        p = figure(title="g%d" % i, plot_width=500, plot_height=300)
        p.line(np.arange(50.0), np.cos(np.arange(50.0) / (5 + i)))
        figs.append(p)
    return gridplot(figs, ncols=2, plot_width=500, plot_height=300,
                    toolbar_location=None)


def twin_plot():
    from kcwidrp.core.kcwi_plotting import figure, Range1d, LinearAxis
    p = figure(title="twin", x_axis_label="Wavelength (A)",
               y_axis_label="Flux", plot_width=800, plot_height=500)
    p.x_range = Range1d(start=4000., end=4100.)
    p.extra_x_ranges = {'pix': Range1d(start=10., end=60.)}
    p.add_layout(LinearAxis(x_range_name='pix', axis_label='CCD Y pix'),
                 'above')
    p.line(np.linspace(4000., 4100., 51), np.linspace(0, 1, 51))
    return p


def check_inside(extents, w, h):
    for left, top, right, bottom in extents:
        assert 0 <= left < right <= w
        assert 0 <= top < bottom <= h


# ------------------------------------------------------------------- page
def fetch_page(viewer):
    with urllib.request.urlopen(viewer.url, timeout=10) as r:
        return r.read().decode()


def test_page_has_vline_checkbox_and_label(viewer):
    html = fetch_page(viewer)
    # source: published - viewer cursor-line spec (maintainer request, 2026-10-28) -- checkbox id "vline", label text as given
    assert 'id="vline"' in html
    assert "Show vertical line at cursor" in html
    assert 'type="checkbox"' in html


def test_page_checkbox_is_unchecked_in_html(viewer):
    html = fetch_page(viewer)
    tag = re.search(r'<input[^>]*id="vline"[^>]*>', html).group(0)
    # source: published - viewer cursor-line spec (maintainer request, 2026-10-28) -- unchecked in the HTML (state comes from localStorage)
    assert "checked" not in tag


# ------------------------------------------------------------- fake tabs
def test_single_figure_sends_extents(viewer):
    tab = observe(viewer, single_plot(640, 400))
    sent = tab.extents()
    assert sent, "no axes_extents message received"
    # source: published - viewer cursor-line spec (maintainer request, 2026-10-28) -- one axes area, inside the figure's pixel size
    assert len(sent[-1]) == 1
    check_inside(sent[-1], 640, 400)


def test_grid_extents_in_four_quadrants(viewer):
    tab = observe(viewer, grid_2x2())
    sent = tab.extents()
    assert sent, "no axes_extents message received"
    ext = sent[-1]
    # source: published - viewer cursor-line spec (maintainer request, 2026-10-28) -- 2x2 grid of 500x300 cells: 1000x600 figure
    assert len(ext) == 4
    check_inside(ext, 1000, 600)
    found = set()
    for left, top, right, bottom in ext:
        col, row = int(left // 500), int(top // 300)
        assert right <= (col + 1) * 500 and bottom <= (row + 1) * 300
        found.add((col, row))
    assert found == {(0, 0), (1, 0), (0, 1), (1, 1)}


def test_twin_axis_gives_one_entry(viewer):
    tab = observe(viewer, twin_plot())
    sent = tab.extents()
    assert sent, "no axes_extents message received"
    # source: published - viewer cursor-line spec (maintainer request, 2026-10-28) -- the twin x-axis shares its parent's area
    assert len(sent[-1]) == 1
    check_inside(sent[-1], 800, 500)


def test_new_plot_sends_new_extents(viewer):
    tab = observe(viewer, single_plot(640, 400), then=grid_2x2())
    sent = tab.extents()
    assert len(sent) >= 2
    first, last = sent[0], sent[-1]
    # source: published - viewer cursor-line spec (maintainer request, 2026-10-28) -- first plot has one area, the 2x2 grid that follows has four
    assert len(first) == 1
    assert len(last) == 4
    assert first != last


# ----------------------------------------------------------------- direct
def test_axes_extents_on_rendered_figure():
    from kcwidrp.core.kcwi_plotting import render
    from kcwidrp.core.plot_viewer import axes_extents
    fig = render(single_plot(640, 400))
    fig.canvas.draw()
    ext = axes_extents(fig)
    # source: published - viewer cursor-line spec (maintainer request, 2026-10-28) -- one entry inside the 640x400 figure, left<right, top<bottom
    assert len(ext) == 1
    check_inside(ext, 640, 400)
