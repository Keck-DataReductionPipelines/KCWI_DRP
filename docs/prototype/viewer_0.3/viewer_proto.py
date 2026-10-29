"""Bokeh-plan task 0.3: throwaway prototype of the live plot viewer.

Reproduces the framework's threading: the main thread idles (with a SIGINT
handler that only sets a flag, like keckdrpframework), a daemon worker thread
builds plot descriptions and blocks in input() or time.sleep(), and a
Tornado/WebAgg server on a third (daemon) thread shows whatever the worker
last submitted.  The process leaves through os._exit, as the framework does.

Design under test:
- one long-lived Figure per viewer, cleared and re-rendered per submission,
  with the toolbar's navigation history reset and a 'resize' sent to every
  connected tab;
- the description is deep-copied on the caller's thread; all Matplotlib
  objects live on the viewer thread; hand-off via IOLoop.add_callback;
- rendering is lazy: a submission only stores the description and asks the
  connected tabs to redraw; with no tab connected nothing is drawn;
- ports base..base+9 on 127.0.0.1, probing with connect() first (macOS lets a
  127.0.0.1 bind shadow a wildcard listener);
- os.register_at_fork closes the listening sockets in forked children.

Usage (demo modes, driven by check_proto.py):
    python viewer_proto.py --port 15006 --mode input|sleep|burst|pool|sizes
"""
import argparse
import asyncio
import copy
import json
import logging
import mimetypes
import os
import signal
import socket
import sys
import threading
import time
import webbrowser
from pathlib import Path

import numpy as np
import tornado.httpserver
import tornado.ioloop
import tornado.netutil
import tornado.web
import tornado.websocket
from matplotlib.backends.backend_webagg_core import (
    FigureCanvasWebAggCore, FigureManagerWebAgg, NavigationToolbar2WebAgg)
from matplotlib.figure import Figure

log = logging.getLogger("proto")

PAGE = """<!DOCTYPE html>
<html lang="en"><head>
<link rel="stylesheet" href="_static/css/page.css" type="text/css">
<link rel="stylesheet" href="_static/css/boilerplate.css" type="text/css">
<link rel="stylesheet" href="_static/css/fbm.css" type="text/css">
<link rel="stylesheet" href="_static/css/mpl.css" type="text/css">
<script src="_static/js/mpl_tornado.js"></script>
<script src="mpl.js"></script>
<script>
document.addEventListener("DOMContentLoaded", function () {
  var ws_type = mpl.get_websocket_type();
  var uri = (window.location.protocol === "https:" ? "wss://" : "ws://")
            + window.location.host + "/ws";
  var websocket = new ws_type(uri);
  websocket.addEventListener("close", function () {
    document.getElementById("status").textContent =
      "Pipeline has stopped or the viewer is not reachable.";
  });
  new mpl.figure(1, websocket, mpl_ondownload,
                 document.getElementById("figure"));
});
</script>
<title>%(title)s</title></head>
<body><div id="status"></div><div id="figure" style="margin: 10px;"></div>
</body></html>
"""


class _Manager(FigureManagerWebAgg):
    _toolbar2_class = NavigationToolbar2WebAgg


class _Canvas(FigureCanvasWebAggCore):
    manager_class = _Manager


def render_description(fig, desc):
    """Draw a prototype description into a cleared Figure."""
    w, h = desc["size"]
    fig.set_size_inches(w / fig.dpi, h / fig.dpi, forward=False)
    nrows, ncols = desc.get("grid", (1, 1))
    for i, panel in enumerate(desc["panels"]):
        ax = fig.add_subplot(nrows, ncols, i + 1)
        ax.set_title(panel["title"], loc="left", fontweight="bold")
        for x, y in panel["lines"]:
            ax.plot(x, y)
        if panel.get("twin"):
            ax2 = ax.twiny()
            ax2.set_xlabel("CCD Y pix")
            ax2.plot(*panel["twin"], color="red")
    fig.suptitle(desc.get("label", ""))


class PlotViewer:
    def __init__(self, port=5006, open_browser=False, title="KCWI DRP"):
        self.base_port = port
        self.open_browser = open_browser
        self.title = title
        self.port = None
        self.url = None
        self.dead = False
        self._loop = None
        self._sockets = []
        self._desc = None
        self._browser_tried = False
        self.n_renders = 0

    # ---------------------------------------------------- caller's thread
    def start(self):
        sockets = self._bind()
        if sockets is None:
            log.warning("Plot viewer: ports %d-%d are all in use; running "
                        "without a viewer", self.base_port, self.base_port + 9)
            self.dead = True
            return False
        self._sockets = sockets
        os.register_at_fork(after_in_child=self._close_in_child)
        ready = threading.Event()
        t = threading.Thread(target=self._run, args=(ready,),
                             name="plot-viewer", daemon=True)
        t.start()
        ready.wait(10)
        self.url = "http://127.0.0.1:%d/" % self.port
        log.info("Plot viewer: %s", self.url)
        print("Plot viewer: %s" % self.url, flush=True)
        if self.open_browser and not self._browser_tried:
            self._browser_tried = True
            try:
                webbrowser.open(self.url)
            except Exception as e:      # headless host
                log.debug("could not open a browser: %s", e)
        return True

    def show(self, desc):
        """Never blocks, never raises."""
        if self.dead or self._loop is None:
            return
        try:
            snap = copy.deepcopy(desc)
            self._loop.add_callback(self._submit, snap)
        except Exception as e:
            log.warning("Plot viewer stopped working (%s); plots continue "
                        "to be saved", e)
            self.dead = True

    def stop(self):
        if self._loop is not None:
            self._loop.add_callback(self._loop.stop)

    def _bind(self):
        for port in range(self.base_port, self.base_port + 10):
            probe = socket.socket()
            probe.settimeout(0.2)
            try:
                if probe.connect_ex(("127.0.0.1", port)) == 0:
                    continue            # someone is listening there
            finally:
                probe.close()
            try:
                socks = tornado.netutil.bind_sockets(port, "127.0.0.1")
            except OSError:
                continue
            for s in socks:
                s.set_inheritable(False)
            self.port = port
            return socks
        return None

    def _close_in_child(self):
        for s in self._sockets:
            try:
                s.close()
            except Exception:
                pass

    # ---------------------------------------------------- viewer thread
    def _run(self, ready):
        asyncio.set_event_loop(asyncio.new_event_loop())
        self._loop = tornado.ioloop.IOLoop.current()
        self.figure = Figure()
        self.canvas = _Canvas(self.figure)
        self.manager = _Manager(self.canvas, 1)
        app = tornado.web.Application([
            (r"/", _Page, dict(viewer=self)),
            (r"/mpl.js", _MplJs),
            (r"/_static/(.*)", tornado.web.StaticFileHandler,
             dict(path=FigureManagerWebAgg.get_static_file_path())),
            (r"/download.([a-z0-9.]+)", _Download, dict(viewer=self)),
            (r"/ws", _Socket, dict(viewer=self)),
        ])
        server = tornado.httpserver.HTTPServer(app)
        server.add_sockets(self._sockets)
        ready.set()
        try:
            self._loop.start()
        except Exception:
            log.exception("plot viewer loop ended")
            self.dead = True

    def _submit(self, desc):
        self._desc = desc
        if self.manager.web_sockets:
            self._draw_current()

    def _draw_current(self):
        if self._desc is None:
            return
        self.figure.clear()
        render_description(self.figure, self._desc)
        self.n_renders += 1
        tb = self.manager.toolbar
        if tb is not None:
            tb.update()                 # drop the previous plot's zoom history
        self.canvas._force_full = True
        self.manager.resize(*self.figure.bbox.size)
        self.figure.set_label(self._desc.get("label", ""))  # sent on refresh
        self.manager.set_window_title(self._desc.get("label", ""))
        self.canvas.draw_idle()


class _Page(tornado.web.RequestHandler):
    def initialize(self, viewer):
        self.viewer = viewer

    def get(self):
        self.write(PAGE % dict(title=self.viewer.title))


class _MplJs(tornado.web.RequestHandler):
    def get(self):
        self.set_header("Content-Type", "application/javascript")
        self.write(FigureManagerWebAgg.get_javascript())


class _Download(tornado.web.RequestHandler):
    def initialize(self, viewer):
        self.viewer = viewer

    def get(self, fmt):
        import io
        buf = io.BytesIO()
        self.viewer.figure.savefig(buf, format=fmt)
        self.set_header("Content-Type",
                        mimetypes.types_map.get("." + fmt, "binary"))
        self.write(buf.getvalue())


class _Socket(tornado.websocket.WebSocketHandler):
    supports_binary = True

    def initialize(self, viewer):
        self.viewer = viewer

    def open(self):
        first = not self.viewer.manager.web_sockets
        self.viewer.manager.add_web_socket(self)
        if first:
            self.viewer._draw_current()
        self.set_nodelay(True)

    def on_close(self):
        self.viewer.manager.remove_web_socket(self)

    def on_message(self, message):
        message = json.loads(message)
        if message["type"] == "supports_binary":
            self.supports_binary = message["value"]
        else:
            self.viewer.manager.handle_json(message)

    def send_json(self, content):
        self.write_message(json.dumps(content))

    def send_binary(self, blob):
        if self.supports_binary:
            self.write_message(blob, binary=True)
        else:
            import base64
            self.write_message("data:image/png;base64,"
                               + base64.b64encode(blob).decode("ascii"))


# ------------------------------------------------------------------ demo
def single(label, n=500, w=1000, h=600, data=None):
    x = np.arange(n, dtype=float)
    y = np.sin(x / 20.0) if data is None else data
    return dict(label=label, size=(w, h),
                panels=[dict(title=label, lines=[(x, y)])])


def grid(label, namps):
    ncols = 2 if namps > 2 else 1
    nrows = (namps + ncols - 1) // ncols
    x = np.arange(200.0)
    return dict(label=label, size=(500 * ncols, 300 * nrows),
                grid=(nrows, ncols),
                panels=[dict(title="amp %d" % i,
                             lines=[(x, np.cos(x / (10 + i)))])
                        for i in range(namps)])


def twin(label):
    x = np.linspace(4000, 5000, 300)
    return dict(label=label, size=(1000, 600),
                panels=[dict(title=label, lines=[(x, np.exp(-(x - 4500) ** 2
                                                            / 200.0))],
                             twin=(np.arange(300.0), np.ones(300)))])


def _square(i):
    return i * i


def worker(viewer, mode, flag):
    from multiprocessing import get_context
    if mode == "input":
        viewer.show(single("before input"))
        print("WAITING-INPUT", flush=True)
        line = sys.stdin.readline()
        print("GOT %r" % line, flush=True)
        viewer.show(single("after input"))
    elif mode == "sleep":
        viewer.show(single("sleeping"))
        print("SLEEPING", flush=True)
        while not flag["stop"]:
            time.sleep(0.05)
    elif mode == "pool":
        method = os.environ.get("PROTO_START", "fork")
        ctx = get_context(method)
        for i in range(50):
            viewer.show(single("pool %d" % i, n=20000))
            p = ctx.Pool()
            r = p.map(_square, range(64))
            p.close()
            p.join()
            assert r[-1] == 63 * 63
        print("POOL-OK %s" % method, flush=True)
        # leave a pool's children alive at exit, as a crash would
        p = ctx.Pool(2)
        p.map_async(time.sleep, [30, 30])
        print("POOL-LEFT", flush=True)
        while not flag["stop"]:
            time.sleep(0.05)
    elif mode == "label":
        viewer.show(single(os.environ.get("PROTO_LABEL", "x")))
        print("SHOWN", flush=True)
        while not flag["stop"]:
            time.sleep(0.05)
    if flag["stop"]:
        print("EXITING", flush=True)
        os._exit(0)
    print("DONE", flush=True)
    os._exit(0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=15006)
    ap.add_argument("--mode", default="sleep")
    ap.add_argument("--open-browser", action="store_true")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                        format="%(levelname)s %(message)s")
    viewer = PlotViewer(a.port, open_browser=a.open_browser)
    viewer.start()
    flag = {"stop": False}

    # As keckdrpframework: SIGINT only sets a flag, the worker exits.
    def handler(sig, frame):
        flag["stop"] = True
    signal.signal(signal.SIGINT, handler)
    t = threading.Thread(target=worker, args=(viewer, a.mode, flag),
                         daemon=True)
    t.start()
    while True:                         # main thread idle, as the framework
        time.sleep(0.2)


if __name__ == "__main__":
    main()
