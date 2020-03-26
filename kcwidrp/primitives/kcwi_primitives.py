
from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments
from keckdrpframework.primitives.base_img import BaseImg
from .kcwi_file_primitives import *
from keckdrpframework.core.bokeh_plotting import bokeh_plot
import ccdproc
from astropy.io import fits as pf
from astropy.nddata import VarianceUncertainty
from astropy.coordinates import SkyCoord
from astropy import units as u
import os

import matplotlib.pyplot as pl
import numpy as np
import scipy as sp
from scipy.signal.windows import boxcar
from scipy.optimize import curve_fit
from bokeh.plotting import figure, show
from bokeh.models import Range1d, LinearAxis
from bokeh.layouts import gridplot
from bokeh.models.markers import X
from bokeh.io import export_png
from bokeh.util.logconfig import basicConfig, bokeh_logger as bl
import logging

from scipy.signal import find_peaks
from skimage import transform as tf
from scipy.interpolate import interpolate
from scipy.signal import find_peaks
from scipy import signal
from scipy.ndimage import gaussian_filter1d
from scipy.stats import sigmaclip
import time
from multiprocessing import Pool
import pkg_resources
import pickle

import pandas as pd

import ref_index

try:
    import _lacosmicx
except ImportError:
    print("Please install lacosmicx from github.com/cmccully/lacosmicx.")
    quit()


def pascal_shift(coef=None, x0=None):
    """Shift coefficients to a new reference value (X0)

    This should probably go somewhere else, but will be needed here.
    """
    if not coef:
        print("Error, no coefficients for pascal_shift.")
        return None
    if not x0:
        print("Warning, no reference value (x0) supplied")
        return coef
    if len(coef) == 7:
        usecoeff = list(reversed(coef))
        fincoeff = [0.] * 7
    else:
        if len(coef) > 7:
            print("Warning - this routine only handles up to 7 coefficients.")
            usecoeff = list(reversed(coef[0:7]))
            fincoeff = [0.] * len(coef)
        else:
            usecoeff = [0.] * 7
            fincoeff = usecoeff
            for ic, c in enumerate(coef):
                usecoeff[len(coef)-(ic+1)] = coef[ic]
    # get reference values
    x01 = x0
    x02 = x0**2
    x03 = x0**3
    x04 = x0**4
    x05 = x0**5
    x06 = x0**6
    # use Pascal's Triangle to shift coefficients
    fincoeff[0] = usecoeff[0] - usecoeff[1] * x01 + usecoeff[2] * x02 \
        - usecoeff[3] * x03 + usecoeff[4] * x04 - usecoeff[5] * x05 \
        + usecoeff[6] * x06

    fincoeff[1] = usecoeff[1] - 2.0 * usecoeff[2] * x01 \
        + 3.0 * usecoeff[3] * x02 - 4.0 * usecoeff[4] * x03 \
        + 5.0 * usecoeff[5] * x04 - 6.0 * usecoeff[6] * x05

    fincoeff[2] = usecoeff[2] - 3.0 * usecoeff[3] * x01 \
        + 6.0 * usecoeff[4] * x02 - 10.0 * usecoeff[5] * x03 \
        + 15.0 * usecoeff[6] * x04

    fincoeff[3] = usecoeff[3] - 4.0 * usecoeff[4] * x01 \
        + 10.0 * usecoeff[5] * x02 - 20.0 * usecoeff[6] * x03

    fincoeff[4] = usecoeff[4] - 5.0 * usecoeff[5] * x01 \
        + 15.0 * usecoeff[6] * x02

    fincoeff[5] = usecoeff[5] - 6.0 * usecoeff[6] * x01

    fincoeff[6] = usecoeff[6]
    # Trim if needed
    if len(coef) < 7:
        fincoeff = fincoeff[0:len(coef)]
    # Reverse for python
    return list(reversed(fincoeff))
    # END: def pascal_shift()


def gaus(x, a, mu, sigma):
    """Gaussian fitting function"""
    return a * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))













































