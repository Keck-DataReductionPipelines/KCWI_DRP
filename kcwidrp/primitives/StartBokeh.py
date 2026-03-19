from keckdrpframework.primitives.base_primitive import BasePrimitive

from bokeh.client import pull_session
from bokeh.plotting.figure import figure
from bokeh.layouts import column

from kcwidrp.core.kcwi_plotting import configure_plot_driver


class StartBokeh(BasePrimitive):
    """
    Start the bokeh server for the KCWI DRP.

    This is enabled through the ``enable_bokeh`` instrument configuration
    parameter in the `kcwi.cfg` file.  Set to ``True`` to enable plots.

    """

    def __init__(self, action, context):
        """
        Constructor
        """
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):

        # session = pull_session(session_id='kcwi', url='http://localhost:5006')
        session = pull_session()
        self.logger.info("Enabling BOKEH plots")
        p = figure()
        c = column(children=[p])
        session.document.clear()
        session.document.add_root(c)
        self.context.bokeh_session = session
        session.show(c)

        firefox_compat = getattr(self.config.instrument, 'plot_firefox_compat', False)
        prewarm = getattr(self.config.instrument, 'plot_prewarm_firefox', False)
        if firefox_compat and prewarm:
            configure_plot_driver(firefox_compat=firefox_compat, prewarm=prewarm)

        return True
