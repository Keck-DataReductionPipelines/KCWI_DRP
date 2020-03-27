from keckdrpframework.primitives.base_primitive import BasePrimitive

from bokeh.client import push_session
from bokeh.io import curdoc
from bokeh.plotting.figure import figure
from bokeh.layouts import column



class StartBokeh(BasePrimitive):

    def __init__(self, action, context):
        '''
        Constructor
        '''
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):

        self.logger.info("Enabling BOKEH plots")

        self.context.bokeh_session = push_session(curdoc())
        p = figure()
        c = column(children=[p])
        curdoc().add_root(c)
        self.context.bokeh_session.show(c)

