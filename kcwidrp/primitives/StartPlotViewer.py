from keckdrpframework.primitives.base_primitive import BasePrimitive

from kcwidrp.core.plot_viewer import PlotViewer


class StartPlotViewer(BasePrimitive):
    """
    Start the live plot viewer for the KCWI DRP.

    Queued by the launch scripts when the viewer is enabled in the
    configuration (event ``start_plot_viewer``, or its alias
    ``start_bokeh``).  The viewer is a small web server on ``127.0.0.1``;
    open the URL it prints in a browser to see the plots of this pipeline
    process.  Configuration (``kcwi.cfg``):

    - ``plot_viewer_port`` (default 5006): first port tried; the next nine
      are tried if it is taken.
    - ``plot_viewer_open_browser`` (default False): open a browser tab.

    If the viewer cannot start, a warning is logged and the reduction
    continues without it; saved plots are not affected.  On success the
    viewer is stored as ``context.bokeh_session`` and its URL as
    ``context.plot_viewer_url``.
    """

    def __init__(self, action, context):
        """
        Constructor
        """
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        existing = getattr(self.context, 'bokeh_session', None)
        if isinstance(existing, PlotViewer) and not existing.dead:
            self.logger.info("Plot viewer already running: %s",
                             existing.url)
            return self.action.args

        self.context.bokeh_session = None
        self.context.plot_viewer_url = None
        try:
            port = int(getattr(self.config.instrument, 'plot_viewer_port',
                               5006))
            open_browser = bool(getattr(self.config.instrument,
                                        'plot_viewer_open_browser', False))
            viewer = PlotViewer(port=port, open_browser=open_browser)
            started = viewer.start()
        except Exception as e:
            self.logger.warning("Could not start the plot viewer, live "
                                "plots disabled: %s", e)
            return self.action.args
        if not started:
            self.logger.warning("Plot viewer not started, live plots "
                                "disabled; plots are still saved")
            return self.action.args

        self.context.bokeh_session = viewer
        self.context.plot_viewer_url = viewer.url
        self.logger.info("Plot viewer: %s", viewer.url)
        rule = "=" * 72
        print("\n%s\n  Plot viewer: %s  (open this in a browser)\n%s\n"
              % (rule, viewer.url, rule), flush=True)
        return self.action.args
