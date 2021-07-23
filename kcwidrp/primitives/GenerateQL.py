from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.kcwi_plotting import save_plot
from bokeh.plotting import figure


import numpy as np

class GenerateQL(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        data = self.action.args.ccddata
        medianed = np.median(data, axis=0) + 1
        medianed_list = medianed.tolist()
        p = figure()
        p = self.remove_plot_features(p)
        w = np.shape(medianed)[1]
        h = np.shape(medianed)[0]
        print(f"X: {w}, Y: {h}")
        p.image([medianed_list], x=0, y=0, dw=w, dh=h,
                             level="image")
        
        save_plot(p, filename="test.png")
    
    def remove_plot_features(self, p):
        """Removes axis and toolbar from an input figure object

        Parameters
        ----------
        p : bokeh.plotting.figure
            Figure object to clean up

        Returns
        -------
        bokeh.plotting.figure
            Figure object with axis and toolbar attributes set to none
        """
        
        p.xgrid.visible = False
        p.ygrid.visible = False
        p.toolbar.logo = None
        p.toolbar_location = None
        p.xaxis.visible = None
        p.yaxis.visible = None
        p.xgrid.grid_line_color = None
        p.ygrid.grid_line_color = None
        return p
