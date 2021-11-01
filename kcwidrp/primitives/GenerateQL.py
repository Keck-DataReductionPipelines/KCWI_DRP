"""
Generates 2D Quick Look images from action.args.ccddata

author: mbrodheim
"""
from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.kcwi_plotting import save_plot
from bokeh.plotting import figure
from bokeh.palettes import gray

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
        """
        Generates a quicklook frame using the method in rti.ql_method
        """

        # Indices get rid of DAR padding pixels
        top = -self.config.rti.DAR_top
        bottom = self.config.rti.DAR_bottom
        left = self.config.rti.DAR_left
        right = -self.config.rti.DAR_right
        scale_factor = self.config.rti.ql_scale

        data = self.action.args.ccddata.data[:,bottom:top,left:right]
        method = self.config.rti.ql_method
        
        self.logger.info(f"Generating {method} quicklook frame")
        # Collapse the image into 2D with the given method
        flattened = self.ql_methods[method](data)
        # Stretch the contrast
        min = np.min(flattened)
        max = np.max(flattened)
        zeroed = flattened - min # Set lowest pixel value to 0
        normalized = zeroed / (max - min) # Make every pixel fall between 0-1
        scaled = np.array(normalized * 255, dtype=np.uint8) # cast to uint8

        # Convert to python list object for bokeh
        flattened_list = scaled.tolist()

        ql_name = "ql_" + Path(self.action.args.name).stem + "_" + method

        # Image dimensions
        w = scaled.shape[1] * scale_factor
        h = scaled.shape[0] * scale_factor
        
        p = figure(plot_width=w, plot_height=h)
        p = self._remove_plot_features(p)
        p.image([flattened_list], x=0, y=0,
                            level="image", dh=1, dw=1,
                            palette=gray(256))
        save_plot(p, filename=ql_name + ".png", width=w*10, height=h*10)
            
        
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