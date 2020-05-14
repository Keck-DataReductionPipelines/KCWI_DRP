from keckdrpframework.primitives.base_primitive import BasePrimitive

from bokeh.client import pull_session
from bokeh.plotting.figure import figure
from bokeh.layouts import column


class StartBokeh(BasePrimitive):

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

        return True
