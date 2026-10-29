"""Automated run of the bokeh-plan 0.3 criteria against viewer_proto.py.

A browser is emulated by a WebSocket client speaking mpl.js's protocol
(supports_binary, send_image_mode, refresh, answering each 'draw' event with
a 'draw' request).  Criteria that need eyes on a real browser are reported
as HUMAN.
"""
import asyncio
import gc
import json
import os
import resource
import shutil
import signal
import socket
import struct
import subprocess
import sys
import tempfile
import threading
import time

import numpy as np
from tornado.websocket import websocket_connect

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import viewer_proto as vp  # noqa: E402

PY = sys.executable
results = []


def record(name, ok, detail=""):
    results.append((name, ok, detail))
    print("%-4s %s %s" % ("PASS" if ok is True else ok or "FAIL", name,
                          detail), flush=True)


class Client:
    """A fake browser tab."""

    def __init__(self, port):
        self.port = port
        self.events = []
        self.frames = []
        self.conn = None

    async def connect(self):
        self.conn = await websocket_connect("ws://127.0.0.1:%d/ws" % self.port)
        for m in ({"type": "supports_binary", "value": True},
                  {"type": "send_image_mode"},
                  {"type": "set_device_pixel_ratio", "device_pixel_ratio": 1},
                  {"type": "refresh"}):
            await self.conn.write_message(json.dumps(m))

    async def pump(self, seconds):
        end = time.time() + seconds
        while time.time() < end:
            try:
                msg = await asyncio.wait_for(self.conn.read_message(),
                                             max(0.01, end - time.time()))
            except asyncio.TimeoutError:
                break
            if msg is None:
                break
            if isinstance(msg, bytes):
                w, h = struct.unpack(">II", msg[16:24])
                self.frames.append((w, h))
                await self.conn.write_message(json.dumps({"type": "ack"}))
            else:
                ev = json.loads(msg)
                self.events.append(ev)
                if ev["type"] == "draw":
                    await self.conn.write_message(json.dumps({"type": "draw"}))

    async def send(self, **m):
        await self.conn.write_message(json.dumps(m))

    def labels(self):
        return [e["label"] for e in self.events
                if e["type"] == "figure_label"]

    def last(self, typ):
        ev = [e for e in self.events if e["type"] == typ]
        return ev[-1] if ev else None


def run(coro):
    return asyncio.new_event_loop().run_until_complete(coro)


def spawn(mode, port, env=None, cwd=None):
    e = dict(os.environ, **(env or {}))
    e.pop("DISPLAY", None)
    return subprocess.Popen([PY, os.path.join(HERE, "viewer_proto.py"),
                             "--port", str(port), "--mode", mode],
                            stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, env=e,
                            cwd=cwd)


def wait_line(proc, needle, timeout=30):
    end = time.time() + timeout
    lines = []
    while time.time() < end:
        line = proc.stdout.readline()
        if not line:
            break
        lines.append(line.rstrip())
        if needle in line:
            return line, lines
    return None, lines


def port_of(line):
    return int(line.rsplit(":", 1)[1].strip().rstrip("/"))


def can_bind(port):
    s = socket.socket()
    s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    try:
        s.bind(("127.0.0.1", port))
        return True
    except OSError:
        return False
    finally:
        s.close()


def listeners(port):
    out = subprocess.run(["lsof", "-nP", "-iTCP:%d" % port, "-sTCP:LISTEN"],
                         capture_output=True, text=True).stdout
    return [l for l in out.splitlines()[1:] if l.strip()]


# ------------------------------------------------------------- criteria
def c_input():
    p = spawn("input", 15106)
    url, _ = wait_line(p, "Plot viewer:")
    wait_line(p, "WAITING-INPUT")
    c = Client(port_of(url))

    async def go():
        await c.connect()
        await c.pump(2.0)
    run(go())
    ok = "before input" in c.labels() and c.frames
    record("plot appears while worker blocked in input()", bool(ok),
           "labels=%s frames=%s" % (c.labels(), c.frames[:1]))
    p.stdin.write("\n")
    p.stdin.flush()
    wait_line(p, "DONE")
    p.wait(10)


def c_burst_inprocess():
    v = vp.PlotViewer(15206)
    v.start()
    c = Client(v.port)

    async def go():
        await c.connect()
        await c.pump(1.0)
        rss0 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        for i in range(200):
            v.show(vp.single("burst %d" % i))
            await c.pump(0.1)
        await c.pump(2.0)
        rss1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return rss0, rss1
    rss0, rss1 = run(go())
    from matplotlib.figure import Figure
    nfig = sum(isinstance(o, Figure) for o in gc.get_objects())
    ok = (c.labels()[-1] == "burst 199" and nfig == 1
          and not v.dead and len(v.figure.axes) == 1)
    record("200 submissions at 0.1 s", ok,
           "renders=%d frames=%d figures=%d maxrss %.0f->%.0f MB" % (
               v.n_renders, len(c.frames), nfig, rss0 / 2**20,
               rss1 / 2**20))
    return v, c


def c_cursor_zoom(v):
    c = Client(v.port)
    v.show(vp.single("zoom me"))

    async def go():
        await c.connect()
        await c.pump(1.0)
        ax = v.figure.axes[0]
        lim0 = ax.get_xlim()
        # cursor readout: move over the middle of the axes
        bb = ax.bbox
        x, y = bb.x0 + bb.width / 2, v.figure.bbox.height - (bb.y0
                                                              + bb.height / 2)
        await c.send(type="motion_notify", x=x, y=y, button=0, buttons=0,
                     modifiers=[])
        await c.pump(0.5)
        msg = c.last("message")
        # zoom: toolbar button, drag a box
        await c.send(type="toolbar_button", name="zoom")
        await c.pump(0.3)
        await c.send(type="button_press", x=x - 100, y=y - 50, button=0,
                     buttons=1, modifiers=[])
        await c.send(type="motion_notify", x=x + 100, y=y + 50, button=0,
                     buttons=1, modifiers=[])
        await c.send(type="button_release", x=x + 100, y=y + 50, button=0,
                     buttons=0, modifiers=[])
        await c.pump(1.0)
        lim1 = ax.get_xlim()
        hist = c.last("history_buttons")
        # new submission of a different size: history reset, resize sent
        v.show(vp.grid("grid", 2))
        await c.pump(1.5)
        hist2 = c.last("history_buttons")
        size = c.last("resize")
        return msg, lim0, lim1, hist, hist2, size
    msg, lim0, lim1, hist, hist2, size = run(go())
    record("cursor x/y readout", bool(msg and "(x, y) = (" in msg["message"]),
           repr(msg and msg["message"]))
    record("zoom changes limits", lim1 != lim0 and
           bool(lim1[1] - lim1[0] < lim0[1] - lim0[0]), "%s -> %s" % (lim0, lim1))
    record("pan", "HUMAN", "same event path as zoom; check in a browser")
    record("new plot clears zoom history",
           bool(hist and hist["Back"]) and not hist2["Back"],
           "before=%s after=%s" % (hist, hist2))
    record("resize sent for 500x600 grid",
           [int(round(s)) for s in size["size"]] == [500, 600], str(size))


def c_sizes(v):
    c = Client(v.port)

    async def go():
        await c.connect()
        await c.pump(0.5)
        out = []
        for d in (vp.single("single"), vp.grid("g4", 4), vp.grid("g2", 2),
                  vp.twin("twin"), vp.single("single again")):
            c.frames.clear()
            v.show(d)
            await c.pump(1.0)
            out.append((d["label"], d["size"], c.frames[0] if c.frames
                        else None, len(v.figure.axes)))
        return out
    out = run(go())
    ok = all(f is not None and tuple(f) == tuple(s) for _, s, f, _ in out)
    record("sizes and layouts follow each submission", ok, str(out))


def c_mutation(v):
    data = np.zeros(500)
    v.show(vp.single("mut", data=data))
    data[:] = 99.0
    c = Client(v.port)

    async def go():
        await c.connect()
        await c.pump(1.0)
    run(go())
    y = v.figure.axes[0].lines[0].get_ydata()
    record("in-place mutation after submit not shown",
           bool(np.all(y == 0)), "max y=%s" % y.max())


def c_reconnect(v):
    v.show(vp.single("current one"))
    a, b = Client(v.port), Client(v.port)

    async def go():
        await a.connect()
        await a.pump(0.8)
        a.conn.close()
        await asyncio.sleep(0.3)
        a2 = Client(v.port)
        await a2.connect()
        await b.connect()
        await asyncio.gather(a2.pump(1.0), b.pump(1.0))
        return a2
    a2 = run(go())
    ok = ("current one" in a2.labels() and "current one" in b.labels()
          and a2.frames and b.frames)
    record("reopened tab and second tab show current plot", bool(ok),
           "%s %s" % (a2.labels(), b.labels()))


def c_two_instances():
    p1 = spawn("label", 15306, env={"PROTO_LABEL": "instance one"})
    u1, _ = wait_line(p1, "Plot viewer:")
    p2 = spawn("label", 15306, env={"PROTO_LABEL": "instance two"})
    u2, _ = wait_line(p2, "Plot viewer:")
    wait_line(p1, "SHOWN")
    wait_line(p2, "SHOWN")
    c1, c2 = Client(port_of(u1)), Client(port_of(u2))

    async def go():
        await c1.connect()
        await c2.connect()
        await asyncio.gather(c1.pump(1.5), c2.pump(1.5))
    run(go())
    ok = (port_of(u1) != port_of(u2) and set(c1.labels()) == {"instance one"}
          and set(c2.labels()) == {"instance two"})
    record("two instances, two ports, own plots", ok,
           "%s %s %s %s" % (u1.strip(), u2.strip(), c1.labels(),
                            c2.labels()))
    for p in (p1, p2):
        p.send_signal(signal.SIGINT)
        p.wait(10)


def c_all_ports_taken():
    held = []
    for port in range(15406, 15416):
        s = socket.socket()
        s.bind(("0.0.0.0" if port % 2 else "127.0.0.1", port))
        s.listen()
        held.append(s)
    p = spawn("label", 15406)
    line, lines = wait_line(p, "SHOWN", 20)
    warn = [l for l in lines if "WARNING" in l and "in use" in l]
    rc = None
    p.send_signal(signal.SIGINT)
    rc = p.wait(10)
    record("all ten ports taken: warning, submissions harmless",
           bool(line and warn) and rc == 0, "%s rc=%s" % (warn, rc))
    for s in held:
        s.close()


def c_no_browser():
    import webbrowser
    calls = []
    orig = webbrowser.open
    webbrowser.open = lambda url: calls.append(url) or False
    try:
        env_display = os.environ.pop("DISPLAY", None)
        v = vp.PlotViewer(15506, open_browser=True)
        v.start()
        v.start  # second start not part of the design; tried once only
        for i in range(20):
            v.show(vp.single("nobody %d" % i))
        time.sleep(1.0)
        ok = (v.n_renders == 0 and len(calls) == 1 and not v.dead)
        record("no tab connected, DISPLAY unset: no renders, one browser "
               "attempt", ok, "renders=%d opens=%d" % (v.n_renders,
                                                        len(calls)))
        v.stop()
        if env_display:
            os.environ["DISPLAY"] = env_display
    finally:
        webbrowser.open = orig


def c_ctrl_c():
    tmp = tempfile.mkdtemp()
    env = {"TMPDIR": tmp}
    p = spawn("sleep", 15606, env=env, cwd=tmp)
    url, _ = wait_line(p, "Plot viewer:")
    wait_line(p, "SLEEPING")
    port = port_of(url)
    c = Client(port)

    async def go():
        await c.connect()
        await c.pump(0.5)
    run(go())
    t0 = time.time()
    p.send_signal(signal.SIGINT)
    try:
        rc = p.wait(10)
    except subprocess.TimeoutExpired:
        rc = "HANG"
    dt = time.time() - t0
    left = os.listdir(tmp)
    ok = (rc == 0 and dt < 3 and not left and can_bind(port)
          and not listeners(port))
    record("Ctrl-C exits via os._exit; nothing left; port rebindable", ok,
           "rc=%s %.2fs files=%s" % (rc, dt, left))
    shutil.rmtree(tmp)


def c_pool(method):
    p = spawn("pool", 15706 if method == "fork" else 15716,
              env={"PROTO_START": method})
    url, _ = wait_line(p, "Plot viewer:")
    port = port_of(url)
    c = Client(port)
    stop = threading.Event()

    def browse():
        async def go():
            await c.connect()
            while not stop.is_set():
                await c.pump(0.2)
        run(go())
    t = threading.Thread(target=browse, daemon=True)
    t.start()
    t0 = time.time()
    ok_line, lines = wait_line(p, "POOL-OK", 300)
    dt = time.time() - t0
    wait_line(p, "POOL-LEFT", 30)
    stop.set()
    t.join(3)
    renders = len(c.frames)
    p.send_signal(signal.SIGINT)
    rc = p.wait(10)
    time.sleep(0.5)
    held = listeners(port)
    record("50 x %s Pool() while serving (renders seen: %d)" % (method,
                                                                renders),
           bool(ok_line) and rc == 0, "%.1fs" % dt)
    record("%s pool children do not hold the port after exit" % method,
           not held and can_bind(port), "listeners=%s" % held)
    subprocess.run(["pkill", "-f", "viewer_proto.py --port %d" % port])


if __name__ == "__main__":
    which = sys.argv[1:] or ["all"]
    c_input()
    v, _ = c_burst_inprocess()
    c_cursor_zoom(v)
    c_sizes(v)
    c_mutation(v)
    c_reconnect(v)
    c_two_instances()
    c_all_ports_taken()
    c_no_browser()
    c_ctrl_c()
    c_pool("fork")
    c_pool("spawn")
    print()
    print("platform:", sys.platform, os.uname().release)
    for name, ok, detail in results:
        print("%-5s %s" % ({True: "PASS", False: "FAIL"}.get(ok, ok), name))
