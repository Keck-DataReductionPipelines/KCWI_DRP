"""Live plot viewer for the KCWI DRP (bokeh-plan task 3.1).

A small web server, built on Matplotlib's ``backend_webagg_core`` with our own
Tornado application, that shows the most recent plot submitted by the
pipeline.  Decisions behind the design are recorded in
``docs/bokeh_port_plan.md`` ("0.3 decision"):

- One daemon thread runs a Tornado ``IOLoop``.  Every Matplotlib object of the
  viewer (one long-lived ``Figure``, its canvas and its manager) is created,
  drawn and destroyed on that thread only.
- :meth:`PlotViewer.show` takes an independent copy of the recorded plot on
  the caller's thread (``kcwi_plotting.snapshot``) and hands the copy to the
  loop with ``IOLoop.add_callback``, the loop's only thread-safe entry point.
  It never waits for the render and never raises.
- Rendering is lazy: a submission only stores the copy; it is drawn when a
  browser tab is connected.  Each new plot clears the one figure, re-renders
  it (``kcwi_plotting.render``), resets the toolbar's zoom history and sends
  a resize to every tab.
- The server listens on ``127.0.0.1`` only, on the first free port from
  ``port`` to ``port + 9``.  A port is skipped when a ``connect()`` to it
  succeeds, because on macOS a ``127.0.0.1`` bind silently shadows a wildcard
  listener.
- The framework leaves through ``os._exit``, so nothing here relies on
  clean-up at exit: daemon thread only, no child processes, no temporary or
  lock files, no ``atexit`` hook and no signal handlers.  The listening
  sockets are not inheritable and are closed in forked children
  (``os.register_at_fork``), so ``multiprocessing`` pool workers never
  service or hold the viewer port.

No ``pyplot`` is imported.  The plotting wrapper is imported lazily, when a
plot is first handled, so this module imports without it.
"""
import asyncio
import base64
import io
import json
import logging
import mimetypes
import os
import socket
import threading
import weakref
import webbrowser

logger = logging.getLogger('KCWI')

#: Number of ports tried, starting at the requested one.
PORT_RANGE = 10
#: Pixel size of the empty figure shown before the first plot and on clear().
EMPTY_SIZE = (600, 600)
#: Dots per inch of the viewer figure; the wrapper maps pixels at 100 dpi.
DPI = 100
#: Seconds start() waits for the viewer thread to come up.
START_TIMEOUT = 10.0
#: Seconds stop() waits for the viewer thread to finish.
STOP_TIMEOUT = 5.0

_ALLOWED_HOSTS = ('127.0.0.1', 'localhost')

# Viewers whose listening sockets must be closed in a forked child.  One
# module-level fork hook serves all of them, so starting many viewers (tests)
# does not pile up hooks that keep the viewers alive.
_LIVE_VIEWERS = weakref.WeakSet()
_FORK_HOOK_REGISTERED = False


def _plotting():
    """The plotting wrapper, imported on first use (tests may replace it)."""
    from kcwidrp.core import kcwi_plotting
    return kcwi_plotting


def _close_sockets_in_child():
    """``os.register_at_fork`` hook: runs in a forked child only."""
    for viewer in list(_LIVE_VIEWERS):
        try:
            viewer._close_sockets_after_fork()
        except Exception:
            pass


def _register_fork_hook():
    global _FORK_HOOK_REGISTERED
    if not _FORK_HOOK_REGISTERED and hasattr(os, 'register_at_fork'):
        os.register_at_fork(after_in_child=_close_sockets_in_child)
        _FORK_HOOK_REGISTERED = True


PAGE = """<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8">
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
      "The pipeline has stopped, or the viewer is not reachable. "
      + "Reload this page when the pipeline is running again.";
  });
  new mpl.figure(1, websocket, mpl_ondownload,
                 document.getElementById("figure"));
});
</script>
<title>%(title)s</title></head>
<body>
<div style="margin: 10px; font-family: sans-serif;">
<b>%(title)s</b> &mdash; shows the latest plot of this pipeline process.
</div>
<div id="status" style="margin: 10px; font-family: sans-serif;
     color: #b00;"></div>
<div id="figure" style="margin: 10px;"></div>
</body></html>
"""


class PlotViewer:
    """Browser view of the latest plot of one pipeline process.

    Parameters
    ----------
    port : int
        First port tried; up to ``port + 9`` are tried in turn.
    open_browser : bool
        Open a browser tab on the URL once, after a successful start.

    Attributes
    ----------
    url : str or None
        ``http://127.0.0.1:<port>/`` once started.
    port : int or None
        The port actually bound.
    dead : bool
        True once the viewer cannot show plots (failed start, viewer thread
        gone, or stopped); :meth:`show` and :meth:`clear` then do nothing.
    """

    def __init__(self, port=5006, open_browser=False):
        self.base_port = int(port)
        self.open_browser = bool(open_browser)
        self.port = None
        self.url = None
        self.dead = False
        self.n_renders = 0          # renders done, for tests
        self._stopped = False
        self._warned = False
        self._browser_tried = False
        self._thread = None
        self._loop = None
        self._server = None
        self._sockets = []
        self._pending = None        # latest snapshot, or None for empty
        self._state_lock = threading.Lock()   # guards start/stop only
        self._warn_lock = threading.Lock()    # one dead-viewer warning
        # created on the viewer thread
        self.figure = None
        self.canvas = None
        self.manager = None

    # ------------------------------------------------------- caller's thread
    def start(self):
        """Bind a port and start the viewer thread.

        Returns True when the viewer is serving.  Never raises: on failure a
        warning is logged, ``dead`` is set and False is returned.
        """
        with self._state_lock:
            if self._thread is not None and not self.dead:
                return True
            if self.dead:
                return False
            try:
                ok = self._start()
            except Exception as e:
                logger.warning("Plot viewer could not start (%s); plots are "
                               "still saved", e)
                logger.debug("Plot viewer start failure", exc_info=True)
                ok = False
            if not ok:
                self._close_listening_sockets()
                self.dead = True
                return False
        if self.open_browser and not self._browser_tried:
            self._browser_tried = True
            try:
                webbrowser.open(self.url)
            except Exception as e:      # headless host, no browser
                logger.debug("Plot viewer: could not open a browser: %s", e)
        return True

    def _start(self):
        sockets = self._bind()
        if not sockets:
            logger.warning("Plot viewer: ports %d-%d on 127.0.0.1 are all "
                           "in use; running without a viewer (plots are "
                           "still saved)", self.base_port,
                           self.base_port + PORT_RANGE - 1)
            return False
        self._sockets = sockets
        _register_fork_hook()
        _LIVE_VIEWERS.add(self)
        ready = threading.Event()
        failure = []
        self._thread = threading.Thread(target=self._run,
                                        args=(ready, failure),
                                        name='kcwi-plot-viewer', daemon=True)
        self._thread.start()
        if not ready.wait(START_TIMEOUT):
            logger.warning("Plot viewer thread did not start within %.0f s; "
                           "running without a viewer", START_TIMEOUT)
            return False
        if failure:
            logger.warning("Plot viewer could not start (%s); running "
                           "without a viewer", failure[0])
            return False
        self.url = "http://127.0.0.1:%d/" % self.port
        logger.debug("Plot viewer serving at %s", self.url)
        return True

    def show(self, obj):
        """Show a recorded plot (``PlotFigure`` or ``PlotGrid``).

        Takes a snapshot on the caller's thread and returns at once.
        Never raises.
        """
        if not self._usable():
            return
        try:
            snap = _plotting().snapshot(obj)
        except Exception as e:
            logger.warning("Plot viewer: could not copy the plot (%s); "
                           "it is not shown", e)
            logger.debug("Plot viewer snapshot failure", exc_info=True)
            return
        self._hand_off(snap)

    def clear(self):
        """Show an empty figure.  Never blocks, never raises."""
        if not self._usable():
            return
        self._hand_off(None)

    def stop(self):
        """Stop the server and the viewer thread (tests, orderly restart).

        Not needed at exit: the thread is a daemon and holds no resources
        outside the process.  Never raises.
        """
        with self._state_lock:
            if self._stopped:
                return
            self._stopped = True
            self.dead = True
            try:
                if self._loop is not None and self._thread is not None \
                        and self._thread.is_alive():
                    self._loop.add_callback(self._shutdown)
                    self._thread.join(STOP_TIMEOUT)
            except Exception as e:
                logger.debug("Plot viewer stop: %s", e)
            self._close_listening_sockets()
            _LIVE_VIEWERS.discard(self)

    def _usable(self):
        if self.dead or self._loop is None:
            return False
        if self._thread is None or not self._thread.is_alive():
            self._mark_dead("the viewer thread is no longer running")
            return False
        return True

    def _hand_off(self, snap):
        try:
            self._loop.add_callback(self._submit, snap)
        except Exception as e:
            self._mark_dead(e)

    def _mark_dead(self, reason):
        """Set ``dead``; warn once unless the viewer was stopped on purpose.

        Called from the caller's thread or, when the loop ends, from the
        viewer thread.
        """
        self.dead = True
        with self._warn_lock:
            if self._stopped or self._warned:
                return
            self._warned = True
        logger.warning("Plot viewer stopped working (%s); the reduction "
                       "continues and plots are still saved", reason)

    def _bind(self):
        import tornado.netutil
        for port in range(self.base_port, self.base_port + PORT_RANGE):
            probe = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            probe.settimeout(0.2)
            try:
                if probe.connect_ex(('127.0.0.1', port)) == 0:
                    continue            # something is listening there
            except OSError:
                pass
            finally:
                probe.close()
            try:
                socks = tornado.netutil.bind_sockets(port, '127.0.0.1')
            except OSError:
                continue
            for s in socks:
                s.set_inheritable(False)
            self.port = port
            return socks
        return None

    def _close_listening_sockets(self):
        for s in self._sockets:
            try:
                s.close()
            except Exception:
                pass

    def _close_sockets_after_fork(self):
        """In a forked child: drop the viewer's sockets (no threads run)."""
        self._close_listening_sockets()
        manager = self.manager
        if manager is None:
            return
        for ws in list(getattr(manager, 'web_sockets', ())):
            try:
                ws.ws_connection.stream.socket.close()
            except Exception:
                pass

    # ---------------------------------------------------------- viewer thread
    def _run(self, ready, failure):
        try:
            import tornado.httpserver
            import tornado.ioloop
            from matplotlib.figure import Figure

            asyncio.set_event_loop(asyncio.new_event_loop())
            self._loop = tornado.ioloop.IOLoop.current()
            self.figure = Figure(figsize=(EMPTY_SIZE[0] / DPI,
                                          EMPTY_SIZE[1] / DPI), dpi=DPI)
            canvas_class, manager_class = _webagg_classes()
            self.canvas = canvas_class(self.figure)
            self.manager = manager_class(self.canvas, 1)
            self._server = tornado.httpserver.HTTPServer(_make_app(self))
            self._server.add_sockets(self._sockets)
        except Exception as e:
            logger.debug("Plot viewer thread failed to start", exc_info=True)
            failure.append(e)
            ready.set()
            return
        ready.set()
        try:
            self._loop.start()
        except Exception:
            logger.debug("Plot viewer loop ended", exc_info=True)
        finally:
            if not self._stopped:
                self._mark_dead("the viewer thread has ended")
            try:
                # Let cancelled writes to closed tabs finish before closing.
                aio = asyncio.get_event_loop()
                pending = asyncio.all_tasks(aio)
                for task in pending:
                    task.cancel()
                if pending:
                    aio.run_until_complete(
                        asyncio.gather(*pending, return_exceptions=True))
                self._loop.close(all_fds=False)
            except Exception:
                pass

    def _submit(self, snap):
        self._pending = snap
        if self.manager.web_sockets:
            self._draw_current()

    def _draw_current(self):
        """Render the latest submission into the one figure and push it."""
        fig = self.figure
        snap = self._pending
        label = ''
        try:
            fig.clear()
            if snap is None:
                self._set_pixel_size(EMPTY_SIZE)
            else:
                wrapper = _plotting()
                wrapper.render(snap, fig)
                size = None
                if hasattr(wrapper, 'pixel_size'):
                    size = wrapper.pixel_size(snap)
                if size:
                    self._set_pixel_size(size)
                label = wrapper.plot_label(snap) or ''
            self.n_renders += 1
        except Exception as e:
            logger.warning("Plot viewer: could not draw the plot (%s); the "
                           "reduction continues", e)
            logger.debug("Plot viewer render failure", exc_info=True)
            try:
                fig.clear()
                self._set_pixel_size(EMPTY_SIZE)
            except Exception:
                pass
            label = ''
        try:
            toolbar = self.manager.toolbar
            if toolbar is not None:
                toolbar.update()        # drop the previous plot's zoom history
            self.canvas._force_full = True
            self.manager.resize(*fig.bbox.size)
            fig.set_label(str(label))   # sent to a tab that reconnects
            self.manager.set_window_title(str(label))
            self.canvas.draw_idle()
        except Exception as e:
            logger.debug("Plot viewer: could not update the tabs: %s", e)

    def _set_pixel_size(self, size):
        """Size the figure in CSS pixels whatever the tab's pixel ratio."""
        width, height = size
        dpi = getattr(self.figure, '_original_dpi', None) or DPI
        self.figure.set_size_inches(width / dpi, height / dpi, forward=False)

    def _shutdown(self):
        try:
            for ws in list(self.manager.web_sockets):
                try:
                    ws.close()
                except Exception:
                    pass
            if self._server is not None:
                self._server.stop()
            self.figure.clear()
        except Exception as e:
            logger.debug("Plot viewer shutdown: %s", e)
        self._loop.stop()


# ------------------------------------------------------------- web plumbing
def _webagg_classes():
    from matplotlib.backends.backend_webagg_core import (
        FigureCanvasWebAggCore, FigureManagerWebAgg, NavigationToolbar2WebAgg)

    class _Manager(FigureManagerWebAgg):
        _toolbar2_class = NavigationToolbar2WebAgg

    class _Canvas(FigureCanvasWebAggCore):
        manager_class = _Manager

    return _Canvas, _Manager


def _make_app(viewer):
    import tornado.web
    import tornado.websocket
    from matplotlib.backends.backend_webagg_core import FigureManagerWebAgg

    class _LocalOnly:
        """Refuse requests whose Host is not this machine (DNS rebinding)."""

        def prepare(self):
            host = self.request.host_name or ''
            if host.strip('[]') not in _ALLOWED_HOSTS:
                raise tornado.web.HTTPError(403)

    class _Page(_LocalOnly, tornado.web.RequestHandler):
        def get(self):
            title = "KCWI DRP plot viewer (port %d)" % viewer.port
            self.write(PAGE % dict(title=title))

    class _MplJs(_LocalOnly, tornado.web.RequestHandler):
        def get(self):
            self.set_header('Content-Type', 'application/javascript')
            self.write(FigureManagerWebAgg.get_javascript())

    class _Static(_LocalOnly, tornado.web.StaticFileHandler):
        pass

    class _Download(_LocalOnly, tornado.web.RequestHandler):
        def get(self, fmt):
            if fmt not in viewer.canvas.get_supported_filetypes():
                raise tornado.web.HTTPError(404)
            buf = io.BytesIO()
            viewer.figure.savefig(buf, format=fmt)
            self.set_header('Content-Type',
                            mimetypes.types_map.get('.' + fmt, 'binary'))
            self.write(buf.getvalue())

    class _Socket(_LocalOnly, tornado.websocket.WebSocketHandler):
        supports_binary = True

        def open(self):
            first = not viewer.manager.web_sockets
            viewer.manager.add_web_socket(self)
            if first:
                viewer._draw_current()
            self.set_nodelay(True)

        def on_close(self):
            try:
                viewer.manager.remove_web_socket(self)
            except Exception:
                pass

        def on_message(self, message):
            try:
                message = json.loads(message)
                if message['type'] == 'supports_binary':
                    self.supports_binary = message['value']
                else:
                    viewer.manager.handle_json(message)
            except Exception as e:
                logger.debug("Plot viewer: bad message from a tab: %s", e)

        # A tab that has just gone away must not stop the others updating.
        def send_json(self, content):
            try:
                self.write_message(json.dumps(content))
            except tornado.websocket.WebSocketClosedError:
                pass

        def send_binary(self, blob):
            try:
                if self.supports_binary:
                    self.write_message(blob, binary=True)
                else:
                    self.write_message(
                        'data:image/png;base64,'
                        + base64.b64encode(blob).decode('ascii'))
            except tornado.websocket.WebSocketClosedError:
                pass

    return tornado.web.Application([
        (r'/', _Page),
        (r'/mpl.js', _MplJs),
        (r'/_static/(.*)', _Static,
         dict(path=FigureManagerWebAgg.get_static_file_path())),
        (r'/download.([a-z0-9.]+)', _Download),
        (r'/ws', _Socket),
    ])
