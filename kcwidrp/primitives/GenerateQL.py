"""
Generates 2D Quick Look images from action.args.ccddata

author: mbrodheim
"""
from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.kcwi_plotting import save_plot
from bokeh.plotting import figure

import numpy as np
from pathlib import Path

class GenerateQL(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

        # Dict with the methods that should be used to generate each 2d image
        self.ql_methods = {
            'median' : self._median,
            'mean' : self._mean,
            'whitelight' : self._whitelight
        }

    def _perform(self):

        data = self.action.args.ccddata
        
        for method in self.ql_methods:
            self.logger.info(f"Generating {method} quicklook frame")
            # Collapse the image into 2D with the given method
            flattened = self.ql_methods[method](data)
            # Bias the image to get rid of negative values (so we can get a png)
            flattened_biased = flattened + np.min(flattened)
            # Convert to python list object for bokeh
            flattened_list = flattened_biased.tolist()

            ql_name = "ql_" + Path(self.action.args.name).stem + \
                        "_" + method + ".png"
            # Image dimensions
            w = np.shape(flattened)[1]
            h = np.shape(flattened)[0]
            
            p = figure(plot_width=w, plot_height=h)
            p = self._remove_plot_features(p)
            p.image([flattened_list], x=0, y=0, dw=w, dh=h,
                                level="image")
            save_plot(p, filename=ql_name)
            
        
        return self.action.args
    
    def _remove_plot_features(self, p):
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


    def _whitelight(self, arr):
        """Calculates sum of all arr elements normalized by median

        Parameters
        ----------
        arr : np.ndarray
            3D array with no NaN values

        Returns
        -------
        np.ndarray
            2D array summed along z axis and normalized over median
        """
        arr = arr + np.min(arr)
        sum = np.sum(arr, axis=0)
        return sum / np.median(sum)


    def _median(self, arr):
        """Calculates median of all arr elements

        Parameters
        ----------
        arr : np.ndarray
            3D array with no NaN values

        Returns
        -------
        np.ndarray
            2D array meaned along z axis
        """

        return np.median(arr, axis=0)
    
    def _mean(self, arr):
        """Calculates mean of all arr elements

        Parameters
        ----------
        arr : np.ndarray
            3D array with no NaN values

        Returns
        -------
        np.ndarray
            2D array meaned along z axis
        """

        return np.mean(arr, axis=0)