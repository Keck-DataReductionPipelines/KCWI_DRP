import numpy as np
from bokeh.io import export_png

import os
import logging

logger = logging.getLogger('KCWI')


def get_plot_lims(data, padding=0.05):
    """Get plot limits using data range plus padding fraction"""
    dmin = np.nanmin(data)
    dmax = np.nanmax(data)
    ddel = dmax - dmin
    return dmin - padding * ddel, dmax + padding * ddel


def oplot_slices(fig, yrange):
    """Overplot slices vertically on plot"""
    for ix in range(1, 24):
        sx = ix * 5 - 0.5
        fig.line([sx, sx], yrange, color='black', line_dash='dashdot')


def set_plot_lims(fig, xlim=None, ylim=None):
    """Set bokeh figure plot ranges"""
    if xlim:
        fig.x_range.start = xlim[0]
        fig.x_range.end = xlim[1]
    if ylim:
        fig.y_range.start = ylim[0]
        fig.y_range.end = ylim[1]


def save_plot(fig, filename=None):
    if filename is None:
        fnam = os.path.join('plots', 'kcwi_drp_plot.png')
    else:
        fnam = os.path.join('plots', filename)
    export_png(fig, filename=fnam)

    logger.info(">>> Saving to %s" % fnam)
