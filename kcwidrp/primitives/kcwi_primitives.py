
from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments
from keckdrpframework.primitives.base_img import BaseImg
from .kcwi_file_primitives import *
from keckdrpframework.core.bokeh_plotting import bokeh_plot
import ccdproc
from astropy.io import fits as pf
from astropy.nddata import VarianceUncertainty
import os

import matplotlib.pyplot as pl
import numpy as np
import scipy as sp
from scipy.signal.windows import boxcar
from scipy.optimize import curve_fit
from bokeh.plotting import figure, show
from bokeh.models import Range1d
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

import pandas as pd

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


def get_line_window(y, c, thresh=0., verbose=False):
    """Find a window that includes the fwhm of the line"""
    nx = len(y)
    # check edges
    if c < 2 or c > nx - 2:
        if verbose:
            print("input center too close to edge")
        return None, None, 0
    # get initial values
    x0 = c - 2
    x1 = c + 2
    mx = np.nanmax(y[x0:x1+1])
    count = 5
    # check low side
    if x0 - 1 < 0:
        if verbose:
            print("max check: low edge hit")
        return None, None, 0
    while y[x0-1] > mx:
        x0 -= 1
        count += 1
        if x0 - 1 < 0:
            if verbose:
                print("Max check: low edge hit")
            return None, None, 0

    # check high side
    if x1 + 1 > nx:
        if verbose:
            print("max check: high edge hit")
        return None, None, 0
    while y[x1+1] > mx:
        x1 += 1
        count += 1
        if x1 + 1 > nx:
            if verbose:
                print("Max check: high edge hit")
            return None, None, 0
    # adjust starting window to center on max
    cmx = x0 + y[x0:x1+1].argmax()
    x0 = cmx - 2
    x1 = cmx + 2
    mx = np.nanmax(y[x0:x1 + 1])
    # make sure max is high enough
    if mx < thresh:
        return None, None, 0
    #
    # expand until we get to half max
    hmx = mx * 0.5
    #
    # Low index side
    prev = mx
    while y[x0] > hmx:
        if y[x0] > mx or x0 <= 0 or y[x0] > prev:
            if verbose:
                if y[x0] > mx:
                    print("hafmax check: low index err - missed max")
                if x0 <= 0:
                    print("hafmax check: low index err - at edge")
                if y[x0] > prev:
                    print("hafmax check: low index err - wiggly")
            return None, None, 0
        prev = y[x0]
        x0 -= 1
        count += 1
    # High index side
    prev = mx
    while y[x1] > hmx:
        if y[x1] > mx or x1 >= nx or y[x1] > prev:
            if verbose:
                if y[x1] > mx:
                    print("hafmax check: high index err - missed max")
                if x1 >= nx:
                    print("hafmax check: high index err - at edge")
                if y[x1] > prev:
                    print("hafmax check: high index err - wiggly")
            return None, None, 0
        prev = y[x1]
        x1 += 1
        count += 1
    # where did we end up?
    if c < x0 or x1 < c:
        if verbose:
            print("initial position outside final window")
        return None, None, 0

    return x0, x1, count
    # END: get_line_window()


def findpeaks(x, y, wid, sth, ath, pkg=None, verbose=False):
    """Find peaks in spectrum"""
    # derivative
    grad = np.gradient(y)
    # smooth derivative
    win = boxcar(wid)
    d = sp.signal.convolve(grad, win, mode='same') / sum(win)
    # size
    nx = len(x)
    # set up windowing
    if not pkg:
        pkg = wid
    hgrp = int(pkg/2)
    hgt = []
    pks = []
    sgs = []
    # loop over spectrum
    # limits to avoid edges given pkg
    for i in np.arange(pkg, (nx - pkg)):
        # find zero crossings
        if np.sign(d[i]) > np.sign(d[i+1]):
            # pass slope threshhold?
            if (d[i] - d[i+1]) > sth * y[i]:
                # pass amplitude threshhold?
                if y[i] > ath or y[i+1] > ath:
                    # get subvectors around peak in window
                    xx = x[(i-hgrp):(i+hgrp+1)]
                    yy = y[(i-hgrp):(i+hgrp+1)]
                    if len(yy) > 3:
                        try:
                            res, _ = curve_fit(gaus, xx, yy,
                                               p0=[100., x[i], 1.])
                            r = abs(x - res[1])
                            t = r.argmin()
                            if abs(i - t) > pkg:
                                if verbose:
                                    print(i, t, x[i], res[1], x[t])
                            else:
                                hgt.append(res[0])
                                pks.append(res[1])
                                sgs.append(abs(res[2]))
                        except RuntimeError:
                            continue
    # clean by sigmas
    cvals = []
    cpks = []
    sgmd = None
    if len(pks) > 0:
        cln_sgs, low, upp = sigmaclip(sgs, low=3., high=3.)
        for i in range(len(pks)):
            if low < sgs[i] < upp:
                cpks.append(pks[i])
                cvals.append(hgt[i])
        # sgmn = cln_sgs.mean()
        sgmd = float(np.nanmedian(cln_sgs))
    else:
        print("No peaks found!")
    return cpks, sgmd, cvals
    # END: findpeaks()


class SubtractOverscan(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        # image sections for each amp
        bsec, dsec, tsec, direc = self.action.args.map_ccd
        namps = len(bsec)
        # polynomial fit order
        if namps == 4:
            porder = 2
        else:
            porder = 7
        # header keyword to update
        key = 'OSCANSUB'
        keycom = 'Overscan subtracted?'
        # is it performed?
        performed = False
        # loop over amps
        # session = self.context.bokeh_session

        for ia in range(namps):
            # get gain
            gain = self.context.data_set.get_info_column(self.action.args.name,
                                                         'GAIN%d' % (ia + 1))
            # check if we have enough data to fit
            if (bsec[ia][3] - bsec[ia][2]) > self.config.instrument.minoscanpix:
                # pull out an overscan vector
                x0 = bsec[ia][2] + self.config.instrument.oscanbuf
                x1 = bsec[ia][3] - self.config.instrument.oscanbuf
                y0 = bsec[ia][0]
                y1 = bsec[ia][1] + 1
                osvec = np.nanmedian(
                    self.action.args.ccddata.data[y0:y1, x0:x1], axis=1)
                nsam = x1 - x0
                xx = np.arange(len(osvec), dtype=np.float)
                # fit it, avoiding first 50 px
                if direc[ia]:
                    # forward read skips first 50 px
                    oscoef = np.polyfit(xx[50:], osvec[50:], porder)
                else:
                    # reverse read skips last 50 px
                    oscoef = np.polyfit(xx[:-50], osvec[:-50], porder)
                # generate fitted overscan vector for full range
                osfit = np.polyval(oscoef, xx)
                # calculate residuals
                resid = (osvec - osfit) * math.sqrt(nsam) * gain / 1.414
                sdrs = float("%.3f" % np.std(resid))
                self.logger.info("Amp%d Read noise from oscan in e-: %.3f" %
                                 ((ia + 1), sdrs))
                self.action.args.ccddata.header['OSCNRN%d' % (ia + 1)] = \
                    (sdrs, "amp%d RN in e- from oscan" % (ia + 1))

                if self.context.config.plot_level >= 1:
                    x = np.arange(len(osvec))
                    p = figure(title=self.action.args.plotlabel +
                               ', Overscan amp %d' % (ia+1),
                               x_axis_label='x', y_axis_label='counts',
                               plot_width=self.config.instrument.plot_width,
                               plot_height=self.config.instrument.plot_height)
                    p.line(x, osvec)
                    p.line(x, osfit, line_color='red', line_width=3)
                    bokeh_plot(p)
                    if self.context.config.plot_level >= 2:
                        input("Next? <cr>: ")
                    else:
                        time.sleep(self.context.config.instrument.plot_pause)
                # subtract it
                for ix in range(dsec[ia][2], dsec[ia][3] + 1):
                    self.action.args.ccddata.data[y0:y1, ix] = \
                        self.action.args.ccddata.data[y0:y1, ix] - osfit
                performed = True
            else:
                self.logger.info("not enough overscan px to fit amp %d")

        self.action.args.ccddata.header[key] = (performed, keycom)

        logstr = self.__module__ + "." + self.__class__.__name__
        self.action.args.ccddata.header['HISTORY'] = logstr
        self.logger.info(logstr)

        return self.action.args
    # END: class SubtractOverscan()


class TrimOverscan(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):

        # parameters
        # image sections for each amp
        bsec, dsec, tsec, direc = self.action.args.map_ccd
        namps = len(bsec)
        # header keyword to update
        key = 'OSCANTRM'
        keycom = 'Overscan trimmed?'
        # get output image dimensions
        max_sec = max(tsec)
        # create new blank image
        new = np.zeros((max_sec[1]+1, max_sec[3]+1), dtype=np.float32)
        # loop over amps
        for ia in range(namps):
            # input range indices
            yi0 = dsec[ia][0]
            yi1 = dsec[ia][1] + 1
            xi0 = dsec[ia][2]
            xi1 = dsec[ia][3] + 1
            # output range indices
            yo0 = tsec[ia][0]
            yo1 = tsec[ia][1] + 1
            xo0 = tsec[ia][2]
            xo1 = tsec[ia][3] + 1
            # transfer to new image
            new[yo0:yo1, xo0:xo1] = self.action.args.ccddata.data[yi0:yi1,
                                                                  xi0:xi1]
            # update amp section
            sec = "[%d:" % (xo0+1)
            sec += "%d," % xo1
            sec += "%d:" % (yo0+1)
            sec += "%d]" % yo1
            self.logger.info("ADDING ATSEC%d" % (ia + 1))
            self.action.args.ccddata.header['ATSEC%d' % (ia+1)] = sec
            # remove obsolete sections
            self.action.args.ccddata.header.pop('ASEC%d' % (ia + 1))
            self.action.args.ccddata.header.pop('BSEC%d' % (ia + 1))
            self.action.args.ccddata.header.pop('DSEC%d' % (ia + 1))
            self.action.args.ccddata.header.pop('CSEC%d' % (ia + 1))
        # update with new image
        self.action.args.ccddata.data = new
        self.action.args.ccddata.header['NAXIS1'] = max_sec[3] + 1
        self.action.args.ccddata.header['NAXIS2'] = max_sec[1] + 1
        self.action.args.ccddata.header[key] = (True, keycom)

        logstr = self.__module__ + "." + self.__class__.__name__
        self.action.args.ccddata.header['HISTORY'] = logstr

        if self.config.instrument.saveintims:
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name, suffix="trim")
        return self.action.args
    # END: class TrimOverscan()


class CorrectGain(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        # Header keyword to update
        key = 'GAINCOR'
        keycom = 'Gain corrected?'
        # print(self.action.args.ccddata.header)
        namps = self.action.args.namps
        for ia in range(namps):
            # get amp section
            section = self.action.args.ccddata.header['ATSEC%d' % (ia + 1)]
            sec, rfor = parse_imsec(section)
            # get gain for this amp
            gain = self.context.data_set.get_info_column(
                self.action.args.name, 'GAIN%d' % (ia + 1))
            self.logger.info(
                "Applying gain correction of %.3f in section %s" %
                (gain, self.action.args.ccddata.header['ATSEC%d' % (ia + 1)]))
            self.action.args.ccddata.data[sec[0]:(sec[1]+1),
                                          sec[2]:(sec[3]+1)] *= gain

        self.action.args.ccddata.header[key] = (True, keycom)
        self.action.args.ccddata.header['BUNIT'] = ('electron', 'Pixel units')
        self.action.args.ccddata.unit = 'electron'

        logstr = self.__module__ + "." + self.__class__.__name__
        self.action.args.ccddata.header['HISTORY'] = logstr

        if self.config.instrument.saveintims:
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name, suffix="gain")
        return self.action.args
    # END: class CorrectGain()


class CorrectDefects(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Correcting detector defects")

        # Header keyword to update
        key = 'BPCLEAN'
        keycom = 'cleaned bad pixels?'

        # Create flags for bad columns fixed
        flags = np.zeros(self.action.args.ccddata.data.shape, dtype=np.uint8)

        # Does the defect file exist?
        path = "data/defect_%s_%dx%d.dat" % (self.action.args.ampmode.strip(),
                                             self.action.args.xbinsize,
                                             self.action.args.ybinsize)
        pkg = __name__.split('.')[0]
        defpath = pkg_resources.resource_filename(pkg, path)
        nbpix = 0   # count of defective pixels cleaned
        if os.path.exists(defpath):
            self.logger.info("Reading defect list in: %s" % defpath)
            deftab = pd.read_csv(defpath, sep='\s+')
            bcdel = 5   # range of pixels for calculating good value
            for indx, row in deftab.iterrows():
                # Get coords and adjust for python zero bias
                x0 = row['X0'] - 1
                x1 = row['X1']
                y0 = row['Y0'] - 1
                y1 = row['Y1']
                # Loop over y range
                for by in range(y0, y1):
                    # sample on low side of bad area
                    vals = list(self.action.args.ccddata.data[by,
                                x0-bcdel:x0])
                    # sample on high side
                    vals.extend(self.action.args.ccddata.data[by,
                                x1+1:x1+bcdel+1])
                    # get replacement value
                    gval = np.nanmedian(np.asarray(vals))
                    # Replace baddies with gval
                    for bx in range(x0, x1):
                        self.action.args.ccddata.data[by, bx] = gval
                        flags[by, bx] += 2
                        nbpix += 1
            self.action.args.ccddata.header[key] = (True, keycom)
            self.action.args.ccddata.header['BPFILE'] = (path, 'defect list')
        else:
            self.logger.error("Defect list not found for %s" % defpath)
            self.action.args.ccddata.header[key] = (False, keycom)

        self.logger.info("Cleaned %d bad pixels" % nbpix)
        self.action.args.ccddata.header['NBPCLEAN'] = \
            (nbpix, 'number of bad pixels cleaned')

        logstr = self.__module__ + "." + self.__class__.__name__
        self.action.args.ccddata.header['HISTORY'] = logstr

        # add flags array
        self.action.args.ccddata.mask = flags
        self.action.args.ccddata.flags = flags

        if self.config.instrument.saveintims:
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name, suffix="def")

        return self.action.args
    # END: class CorrectDefects()


class RemoveCosmicRays(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        # TODO: implement options from kcwi_stage1.pro
        self.logger.info("Finding and masking cosmic rays")

        # Header keyword to update
        key = 'CRCLEAN'
        keycom = 'cosmic rays cleaned?'

        header = self.action.args.ccddata.header
        namps = header['NVIDINP']
        read_noise = 0.
        for ia in range(namps):
            read_noise += header['BIASRN%d' % (ia + 1)]
        read_noise /= float(namps)

        # Set sigclip according to image parameters
        sigclip = self.context.config.instrument.CRR_SIGCLIP
        if 'FLATLAMP' in self.action.args.ccddata.header['IMTYPE']:
            if self.action.args.nasmask:
                sigclip = 10.
            else:
                sigclip = 7.
        if 'OBJECT' in self.action.args.ccddata.header['IMTYPE']:
            if self.action.args.ccddata.header['TTIME'] < 300.:
                sigclip = 10.

        if header['TTIME'] >= self.context.config.instrument.CRR_MINEXPTIME:
            mask, clean = _lacosmicx.lacosmicx(
                self.action.args.ccddata.data, gain=1.0, readnoise=read_noise,
                psffwhm=self.context.config.instrument.CRR_PSFFWHM,
                sigclip=sigclip,
                sigfrac=self.context.config.instrument.CRR_SIGFRAC,
                objlim=self.context.config.instrument.CRR_OBJLIM,
                fsmode=self.context.config.instrument.CRR_FSMODE,
                psfmodel=self.context.config.instrument.CRR_PSFMODEL,
                verbose=self.context.config.instrument.CRR_VERBOSE,
                sepmed=self.context.config.instrument.CRR_SEPMED,
                cleantype=self.context.config.instrument.CRR_CLEANTYPE)

            header['history'] = "LA CosmicX: cleaned cosmic rays"
            header[
                'history'] = "LA CosmicX params: " \
                             "sigclip=%5.2f sigfrac=%5.2f objlim=%5.2f" % (
                self.context.config.instrument.CRR_SIGCLIP,
                self.context.config.instrument.CRR_SIGFRAC,
                self.context.config.instrument.CRR_OBJLIM)
            header[
                'history'] = "LA CosmicX params: " \
                             "fsmode=%s psfmodel=%s psffwhm=%5.2f" % (
                self.context.config.instrument.CRR_FSMODE,
                self.context.config.instrument.CRR_PSFMODEL,
                self.context.config.instrument.CRR_PSFFWHM)
            header['history'] = "LA CosmicX params: sepmed=%s minexptime=%f" % (
                self.context.config.instrument.CRR_SEPMED,
                self.context.config.instrument.CRR_MINEXPTIME)
            # header['history'] = "LA CosmicX run on %s" % time.strftime("%c")

            mask = np.cast["bool"](mask)
            n_crs = mask.sum()
            self.action.args.ccddata.header[key] = (True, keycom)
            self.action.args.ccddata.header['NCRCLEAN'] = (
                n_crs, "number of cosmic ray pixels")
            self.action.args.ccddata.mask += mask
            self.action.args.ccddata.data = clean
        else:
            header[
                'history'] = "LA CosmicX: exptime < minexptime=%.1f" % \
                             self.context.config.instrument.CRR_MINEXPTIME

        logstr = self.__module__ + "." + self.__class__.__name__
        self.action.args.ccddata.header['HISTORY'] = logstr

        if self.config.instrument.saveintims:
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name, suffix="crr")

        return self.action.args
    # END: class RemoveCosmicRays()


class CreateUncertaintyImage(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        """Assumes units of image are electron"""

        # Header keyword to update
        key = 'UNCVAR'
        keycom = 'variance created?'

        self.logger.info("Create uncertainty image")
        # start with Poisson noise
        self.action.args.ccddata.uncertainty = VarianceUncertainty(
            self.action.args.ccddata.data, unit='electron^2', copy=True)
        # add readnoise, if known
        if 'BIASRN1' in self.action.args.ccddata.header:
            namps = self.action.args.ccddata.header['NVIDINP']
            for ia in range(namps):
                # get amp parameters
                biasrn = self.action.args.ccddata.header['BIASRN%d' % (ia + 1)]
                section = self.action.args.ccddata.header['ATSEC%d' % (ia + 1)]
                sec, rfor = parse_imsec(section)
                self.action.args.ccddata.uncertainty.array[
                    sec[0]:(sec[1]+1), sec[2]:(sec[3]+1)] += biasrn
        else:
            self.logger.warn("Readnoise undefined, uncertainty Poisson only")
        # document variance image creation
        self.action.args.ccddata.header[key] = (True, keycom)

        logstr = self.__module__ + "." + self.__class__.__name__
        self.action.args.ccddata.header['HISTORY'] = logstr

        return self.action.args
    # END: class CreateUncertaintyImage()


class RectifyImage(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):

        # Header keyword to update
        key = 'IMGRECT'
        keycom = 'Image rectified?'

        # get amp mode
        ampmode = self.action.args.ccddata.header['AMPMODE'].strip().upper()

        if '__B' in ampmode or '__G' in ampmode:
            newimg = np.rot90(self.action.args.ccddata.data, 2)
            self.action.args.ccddata.data = newimg
            if self.action.args.ccddata.uncertainty:
                newunc = np.rot90(self.action.args.ccddata.uncertainty.array, 2)
                self.action.args.ccddata.uncertainty.array = newunc
        elif '__D' in ampmode or '__F' in ampmode:
            newimg = np.fliplr(self.action.args.ccddata.data)
            self.action.args.ccddata.data = newimg
            if self.action.args.ccddata.uncertainty:
                newunc = np.fliplr(self.action.args.ccddata.uncertainty.array)
                self.action.args.ccddata.uncertainty.array = newunc
        elif '__A' in ampmode or '__H' in ampmode or 'TUP' in ampmode:
            newimg = np.flipud(self.action.args.ccddata.data)
            self.action.args.ccddata.data = newimg
            if self.action.args.ccddata.uncertainty:
                newunc = np.flipud(self.action.args.ccddata.uncertainty.array)
                self.action.args.ccddata.uncertainty.array = newunc

        self.action.args.ccddata.header[key] = (True, keycom)

        logstr = self.__module__ + "." + self.__class__.__name__
        self.action.args.ccddata.header['HISTORY'] = logstr

        # write out int image
        kcwi_fits_writer(self.action.args.ccddata,
                         table=self.action.args.table,
                         output_file=self.action.args.name, suffix="int")
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix="int")
        self.context.proctab.write_proctab()
        return self.action.args
    # END: class RectifyImage()


class SubtractBias(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):

        # Header keyword to update
        key = 'BIASSUB'
        keycom = 'master bias subtracted?'

        self.logger.info("Subtracting master bias")
        tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                             target_type='MBIAS',
                                             nearest=True)
        self.logger.info("%d master bias frames found" % len(tab))

        if len(tab) > 0:
            mbname = tab['OFNAME'][0].split('.')[0] + "_master_bias.fits"
            print("*************** READING IMAGE: %s" % mbname)
            mbias = kcwi_fits_reader(
                os.path.join(os.path.dirname(self.action.args.name), 'redux',
                             mbname))[0]

            # do the subtraction
            self.action.args.ccddata.data -= mbias.data

            # transfer bias read noise
            namps = self.action.args.ccddata.header['NVIDINP']
            for ia in range(namps):
                self.action.args.ccddata.header['BIASRN%d' % (ia + 1)] = \
                    mbias.header['BIASRN%d' % (ia + 1)]

            self.action.args.ccddata.header[key] = (True, keycom)
            self.action.args.ccddata.header['MBFILE'] = (mbname,
                                                         "Master bias filename")
        else:

            self.action.args.ccddata.header[key] = (False, keycom)

        logstr = self.__module__ + "." + self.__class__.__name__
        self.action.args.ccddata.header['HISTORY'] = logstr

        return self.action.args
    # END: class SubtractBias()


class SubtractDark(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):

        # Header keyword to update
        key = 'DARKSUB'
        keycom = 'master dark subtracted?'

        self.logger.info("Subtracting master dark")
        tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                             target_type='MDARK',
                                             nearest=True)
        self.logger.info("%d master dark frames found" % len(tab))

        if len(tab) > 0:
            mdname = tab['OFNAME'][0].split('.')[0] + "_master_dark.fits"
            print("*************** READING IMAGE: %s" % mdname)
            mdark = kcwi_fits_reader(
                os.path.join(os.path.dirname(self.action.args.name), 'redux',
                             mdname))[0]
            # scale by exposure time
            fac = 1.0
            if 'TTIME' in mdark.header and \
               'TTIME' in self.action.args.ccddata.header:
                fac = float(self.action.args.ccddata.header['TTIME']) / \
                      float(mdark.header['TTIME'])
                self.logger.info("dark scaled by %.3f" % fac)
            else:
                self.logger.warn("unable to scale dark by exposure time")

            # do the subtraction
            self.action.args.ccddata.data -= mdark.data * fac

            self.action.args.ccddata.header[key] = (True, keycom)
            self.action.args.ccddata.header['MDFILE'] = (mdname,
                                                         "Master dark filename")
            self.action.args.ccddata.header['DARKSCL'] = (fac,
                                                          "dark scale factor")
        else:
            self.logger.info("No master dark frame available, skipping")
            self.action.args.ccddata.header[key] = (False, keycom)

        logstr = self.__module__ + "." + self.__class__.__name__
        self.action.args.ccddata.header['HISTORY'] = logstr

        # Below now done in subtract_scattered_light()
        # write out int image
        # kcwi_fits_writer(self.action.args.ccddata,
        #                 table=self.action.args.table,
        #                 output_file=self.action.args.name, suffix="intd")
        # self.context.proctab.update_proctab(frame=self.action.args.ccddata,
        #                                    suffix="intd")
        # self.context.proctab.write_proctab()

        return self.action.args
    # END: class SubtractDark()


class ProcessBias(BaseImg):

    def __init__(self, action, context):
        BaseImg.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def hello(self):
        """
        Checks if we can build a stacked bias frame
        Expected arguments:
            want_type: ie. BIAS
            min_files: ie 10
            new_type: ie MASTER_BIAS
            new_file_name: master_bias.fits

        """
        try:
            args = self.action.args
            df = self.context.data_set.data_table
            files = df[(df.IMTYPE == args.want_type) &
                       (df.GROUPID == args.groupid)]
            nfiles = len(files)

            self.logger.info(f"pre condition got {nfiles},"
                             f" expecting {args.min_files}")
            if nfiles < 1 or nfiles < args.min_files:
                return False
            return True
        except Exception as e:
            self.logger.error(f"Exception in base_ccd_primitive: {e}")
            return False

    def _pre_condition(self):
        """
        Checks if we can build a stacked frame based on the processing table
        :return:
        """
        # Add to proctab
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix='RAW')
        self.context.proctab.write_proctab()
        # Get bias count
        self.logger.info("Checking precondition for process_bias")
        self.combine_list = self.context.proctab.n_proctab(
            frame=self.action.args.ccddata, target_type='BIAS',
            target_group=self.action.args.groupid)
        self.logger.info(f"pre condition got {len(self.combine_list)},"
                         f" expecting {self.action.args.min_files}")
        # Did we meet our pre-condition?
        if len(self.combine_list) >= self.action.args.min_files:
            return True
        else:
            return False

    def _perform(self):
        """
        Returns an Argument() with the parameters that depends on this operation
        """
        args = self.action.args
        method = 'average'
        suffix = 'master_bias'

        combine_list = list(self.combine_list['OFNAME'])
        # get master bias output name
        mbname = combine_list[0].split('.fits')[0] + '_master_bias.fits'
        stack = []
        stackf = []
        for bias in combine_list:
            stackf.append(bias)
            # using [0] drops the table
            stack.append(kcwi_fits_reader(bias)[0])

        stacked = ccdproc.combine(stack, method=method, sigma_clip=True,
                                  sigma_clip_low_thresh=None,
                                  sigma_clip_high_thresh=2.0)
        stacked.header.IMTYPE = args.new_type
        stacked.header['NSTACK'] = (len(combine_list),
                                    'number of images stacked')
        stacked.header['STCKMETH'] = (method, 'method used for stacking')
        for ii, fname in enumerate(stackf):
            stacked.header['STACKF%d' % (ii + 1)] = (fname, "stack input file")

        # for readnoise stats use 2nd and 3rd bias
        diff = stack[1].data.astype(np.float32) - \
            stack[2].data.astype(np.float32)
        namps = stack[1].header['NVIDINP']
        for ia in range(namps):
            # get gain
            gain = stacked.header['GAIN%d' % (ia + 1)]
            # get amp section
            sec, rfor = parse_imsec(stacked.header['DSEC%d' % (ia + 1)])
            noise = diff[sec[0]:(sec[1]+1), sec[2]:(sec[3]+1)]
            noise = np.reshape(noise, noise.shape[0]*noise.shape[1]) * \
                gain / 1.414
            # get stats on noise
            c, upp, low = sigmaclip(noise, low=3.5, high=3.5)
            bias_rn = c.std()
            self.logger.info("Amp%d read noise from bias in e-: %.3f" %
                             ((ia + 1), bias_rn))
            stacked.header['BIASRN%d' % (ia + 1)] = \
                (float("%.3f" % bias_rn), "RN in e- from bias")

        logstr = self.__module__ + "." + self.__class__.__name__
        stacked.header['HISTORY'] = logstr

        kcwi_fits_writer(stacked, output_file=mbname)
        self.context.proctab.update_proctab(frame=stacked, suffix=suffix,
                                            newtype=args.new_type)
        self.context.proctab.write_proctab()
        return Arguments(name=mbname)
    # END: class ProcessBias()


class StackDarks(BaseImg):

    def __init__(self, action, context):
        BaseImg.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if we can build a stacked frame based on the processing table
        :return:
        """
        # get current group id
        self.logger.info("Checking precondition for stack_darks")
        self.combine_list = self.context.proctab.n_proctab(
            frame=self.action.args.ccddata, target_type='DARK',
            target_group=self.action.args.groupid)
        self.logger.info(f"pre condition got {len(self.combine_list)},"
                         f" expecting {self.action.args.min_files}")
        # create master dark
        if len(self.combine_list) >= self.action.args.min_files:
            return True
        else:
            return False

    def _perform(self):
        """
        Returns an Argument() with the parameters that depends on this operation
        """
        args = self.action.args
        method = 'average'
        suffix = 'master_dark'

        combine_list = list(self.combine_list['OFNAME'])
        # get master dark output name
        mdname = combine_list[0].split('.fits')[0] + '_master_dark.fits'
        stack = []
        stackf = []
        for dark in combine_list:
            # get dark intensity (int) image file name in redux directory
            stackf.append(dark.split('.fits')[0] + '_int.fits')
            darkfn = os.path.join(args.in_directory, stackf[-1])
            # using [0] gets just the image data
            stack.append(kcwi_fits_reader(darkfn)[0])

        stacked = ccdproc.combine(stack, method=method, sigma_clip=True,
                                  sigma_clip_low_thresh=None,
                                  sigma_clip_high_thresh=2.0)
        stacked.unit = stack[0].unit
        stacked.header.IMTYPE = args.new_type
        stacked.header['NSTACK'] = (len(combine_list),
                                    'number of images stacked')
        stacked.header['STCKMETH'] = (method, 'method used for stacking')
        for ii, fname in enumerate(stackf):
            stacked.header['STACKF%d' % (ii + 1)] = (fname, "stack input file")

        logstr = self.__module__ + "." + self.__class__.__name__
        stacked.header['HISTORY'] = logstr

        kcwi_fits_writer(stacked, output_file=mdname)
        self.context.proctab.update_proctab(frame=stacked, suffix=suffix,
                                            newtype=args.new_type)
        self.context.proctab.write_proctab()
        return Arguments(name=mdname)
    # END: class StackDarks()


class SubtractScatteredLight(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):

        # Header keyword to update
        key = 'SCATSUB'
        keycom = "was scattered light subtracted?"

        # Skip if nod-and-shuffle
        if self.action.args.nasmask:
            self.logger.info("NAS Mask: skipping scattered light subtraction")
            self.action.args.ccddata.header[key] = (False, keycom)
        elif self.context.config.instrument.skipscat:
            self.logger.info("Skipping scattered light subtraction by request")
            self.action.args.ccddata.header[key] = (False, keycom)
        else:
            # Get size of image
            siz = self.action.args.ccddata.data.shape
            # Get x range for scattered light
            x0 = int(siz[1] / 2 - 180 / self.action.args.xbinsize)
            x1 = int(siz[1] / 2 + 180 / self.action.args.xbinsize)
            # Get y limits
            y0 = 0
            # y1 = int(siz[0] / 2 - 1)
            # y2 = y1 + 1
            y3 = siz[0]
            # print("x limits: %d, %d, y limits: %d, %d" % (x0, x1, y0, y3))
            # Y data values
            yvals = np.nanmedian(self.action.args.ccddata.data[y0:y3, x0:x1],
                                 axis=1)
            # X data values
            xvals = np.arange(len(yvals), dtype=np.float)
            # Break points
            nbkpt = int(siz[1] / 40.)
            bkpt = xvals[nbkpt:-nbkpt:nbkpt]
            # B-spline fit
            bspl = sp.interpolate.LSQUnivariateSpline(xvals, yvals, bkpt)
            if self.context.config.plot_level >= 1:
                # plot
                p = figure(title=self.action.args.plotlabel +
                           ", Scattered Light",
                           x_axis_label='y pixel', y_axis_label='e-',
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.circle(xvals, yvals, color='red', legend="Scat")
                xx = np.linspace(0, max(xvals), len(yvals) * 5)
                p.line(xx, bspl(xx), color='blue', line_width=3, legend="fit")
                bokeh_plot(p)
                if self.context.config.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.context.config.instrument.plot_pause)
            # Scattered light vector
            scat = bspl(xvals)
            # Subtract scattered light
            self.logger.info("Starting scattered light subtraction")
            for ix in range(0, siz[1]):
                self.action.args.ccddata.data[y0:y3, ix] = \
                    self.action.args.ccddata.data[y0:y3, ix] - scat
            self.action.args.ccddata.header[key] = (True, keycom)

        logstr = self.__module__ + "." + self.__class__.__name__
        self.action.args.ccddata.header['HISTORY'] = logstr

        # write out int image
        kcwi_fits_writer(self.action.args.ccddata,
                         table=self.action.args.table,
                         output_file=self.action.args.name, suffix="intd")
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix="intd")
        self.context.proctab.write_proctab()

        return self.action.args


class StackFlats(BaseImg):

    def __init__(self, action, context):
        BaseImg.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if we can build a stacked frame based on the processing table
        :return:
        """
        # get current group id
        self.logger.info("Checking precondition for stack_flats")
        self.combine_list = self.context.proctab.n_proctab(
            frame=self.action.args.ccddata, target_type='FLATLAMP',
            target_group=self.action.args.groupid)
        self.logger.info(f"pre condition got {len(self.combine_list)},"
                         f" expecting {self.action.args.min_files}")
        # create master flat
        if len(self.combine_list) >= self.action.args.min_files:
            return True
        else:
            return False

    def _perform(self):
        """
        Returns an Argument() with the parameters that depends on this operation
        """
        args = self.action.args
        method = 'average'
        suffix = 'flat_stack'

        combine_list = list(self.combine_list['OFNAME'])
        # get master dark output name
        mdname = combine_list[0].split('.fits')[0] + '_flat_stack.fits'
        stack = []
        stackf = []
        for flat in combine_list:
            # get dark intensity (int) image file name in redux directory
            stackf.append(flat.split('.fits')[0] + '_intd.fits')
            flatfn = os.path.join(args.in_directory, stackf[-1])
            # using [0] gets just the image data
            stack.append(kcwi_fits_reader(flatfn)[0])

        stacked = ccdproc.combine(stack, method=method, sigma_clip=True,
                                  sigma_clip_low_thresh=None,
                                  sigma_clip_high_thresh=2.0)
        stacked.header.IMTYPE = args.new_type
        stacked.header['NSTACK'] = (len(combine_list),
                                    'number of images stacked')
        stacked.header['STCKMETH'] = (method, 'method used for stacking')
        for ii, fname in enumerate(stackf):
            stacked.header['STACKF%d' % (ii + 1)] = (fname, "stack input file")

        logstr = self.__module__ + "." + self.__class__.__name__
        stacked.header['HISTORY'] = logstr

        kcwi_fits_writer(stacked, output_file=mdname)
        self.context.proctab.update_proctab(frame=stacked, suffix=suffix,
                                            newtype=args.new_type)
        self.context.proctab.write_proctab()
        return Arguments(name=mdname)
    # END: class StackFlats()


class ProcessDark(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        return self.action.args


class ProcessContbars(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        return self.action.args


class ProcessArc(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        return self.action.args


class ProcessFlat(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        return self.action.args


class ProcessObject(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        return self.action.args


class FindBars(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        basicConfig(level=logging.ERROR)

    def _perform(self):
        self.logger.info("Finding continuum bars")
        # Do we plot?
        if self.context.config.plot_level >= 1:
            do_plot = True
        else:
            do_plot = False
        # initialize
        refbar = self.context.config.instrument.REFBAR
        midcntr = []
        # get image dimensions
        nx = self.action.args.ccddata.data.shape[1]
        ny = self.action.args.ccddata.data.shape[0]
        # get binning
        ybin = self.action.args.ybinsize
        win = int(10 / ybin)
        # select from center rows of image
        midy = int(ny / 2)
        midvec = np.median(
            self.action.args.ccddata.data[(midy-win):(midy+win+1), :], axis=0)
        # set threshold for peak finding
        midavg = np.average(midvec)
        self.logger.info("peak threshold = %f" % midavg)
        # find peaks above threshold
        midpeaks, _ = find_peaks(midvec, height=midavg)
        # do we have the requisite number?
        if len(midpeaks) != self.context.config.instrument.NBARS:
            self.logger.error("Did not find %d peaks: n peaks = %d" %
                              (self.config.instrument.NBARS, len(midpeaks)))
        else:
            self.logger.info("found %d bars" % len(midpeaks))
            if do_plot:
                # plot the peak positions
                x = np.arange(len(midvec))
                # pl.plot(midvec, '-')
                p = figure(
                    title=self.action.args.plotlabel +
                    ", Thresh = %.2f" % midavg,
                    x_axis_label='CCD X (px)', y_axis_label='e-',
                    plot_width=self.context.config.instrument.plot_width,
                    plot_height=self.context.config.instrument.plot_height
                )
                p.line(x, midvec, color='blue')
                p.scatter(midpeaks, midvec[midpeaks], marker='x', color='red')
                p.line([0, nx], [midavg, midavg], color='grey',
                       line_dash='dashed')
                bokeh_plot(p)
                time.sleep(self.context.config.instrument.plot_pause)
                # calculate the bar centroids
                plotting_vector_x = []
                plotting_vector_y = []
                for peak in midpeaks:
                    xs = list(range(peak-win, peak+win+1))
                    ys = midvec[xs] - np.nanmin(midvec[xs])
                    xc = np.sum(xs*ys) / np.sum(ys)
                    midcntr.append(xc)
                    plotting_vector_x.append(xc)
                    plotting_vector_y.append(midavg)
                    plotting_vector_x.append(xc)
                    plotting_vector_y.append(midvec[peak])
                    plotting_vector_x.append(xc)
                    plotting_vector_y.append(midavg)
                p.line(plotting_vector_x, plotting_vector_y, color='grey')
                # p.line([xc, xc], [midavg, midvec[peak]], color='grey')
                p.scatter(midcntr, midvec[midpeaks], marker='x', color='green')
                bokeh_plot(p)
                if self.context.config.plot_level >= 2:
                    input("next: ")
                else:
                    time.sleep(self.context.config.instrument.plot_pause)
            self.logger.info("Found middle centroids for continuum bars")
        # store peaks
        self.action.args.midcntr = midcntr
        # store the row where we got them
        self.action.args.midrow = midy
        self.action.args.win = win
        # calculate reference delta x based on refbar
        self.action.args.refdelx = 0.
        for ib in range(refbar-1, refbar+3):
            self.action.args.refdelx += (midcntr[ib] - midcntr[ib-1])
        self.action.args.refdelx /= 4.
        # store image info
        self.action.args.cbarsno = self.action.args.ccddata.header['FRAMENO']
        self.action.args.cbarsfl = self.action.args.ccddata.header['OFNAME']
        return self.action.args
    # END: class FindBars()


class TraceBars(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Tracing continuum bars")
        if self.context.config.plot_level >= 1:
            do_plot = True
            pl.ion()
        else:
            do_plot = False
        if len(self.action.args.midcntr) < 1:
            self.logger.error("No bars found")
        else:
            # initialize
            samp = int(80 / self.action.args.ybinsize)
            win = self.action.args.win
            xi = []     # x input
            xo = []     # x output
            yi = []     # y input (and output)
            barid = []  # bar id number
            slid = []   # slice id number
            # loop over bars
            for barn, barx in enumerate(self.action.args.midcntr):
                # nearest pixel to bar center
                barxi = int(barx + 0.5)
                self.logger.info("bar number %d is at %.3f" % (barn, barx))
                # middle row data
                xi.append(barx)
                xo.append(barx)
                yi.append(self.action.args.midrow)
                barid.append(barn)
                slid.append(int(barn/5))
                # trace up
                samy = self.action.args.midrow + samp
                done = False
                while samy < (self.action.args.ccddata.data.shape[0] - win) \
                        and not done:
                    ys = np.median(
                        self.action.args.ccddata.data[(samy - win):
                                                      (samy + win + 1),
                                                      (barxi - win):
                                                      (barxi + win + 1)],
                        axis=0)
                    ys = ys - np.nanmin(ys)
                    xs = list(range(barxi - win, barxi + win + 1))
                    xc = np.sum(xs * ys) / np.sum(ys)
                    if np.nanmax(ys) > 255:
                        xi.append(xc)
                        xo.append(barx)
                        yi.append(samy)
                        barid.append(barn)
                        slid.append(int(barn/5))
                    else:
                        done = True
                    samy += samp
                # trace down
                samy = self.action.args.midrow - samp
                done = False
                while samy >= win and not done:
                    ys = np.median(
                        self.action.args.ccddata.data[(samy - win):
                                                      (samy + win + 1),
                                                      (barxi - win):
                                                      (barxi + win + 1)],
                        axis=0)
                    ys = ys - np.nanmin(ys)
                    xs = list(range(barxi - win, barxi + win + 1))
                    xc = np.sum(xs * ys) / np.sum(ys)
                    if np.nanmax(ys) > 255:
                        xi.append(xc)
                        xo.append(barx)
                        yi.append(samy)
                        barid.append(barn)
                        slid.append(int(barn / 5))
                    else:
                        done = True
                    # disable for now
                    samy -= samp
            # end loop over bars
            # create source and destination coords
            yo = yi
            dst = np.column_stack((xi, yi))
            src = np.column_stack((xo, yo))
            if do_plot:
                # plot them
                # pl.ioff()
                p = figure(title=self.action.args.plotlabel,
                           x_axis_label="CCD X (px)", y_axis_label="CCD Y (px)",
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.scatter(xi, yi, marker='x', size=2, color='blue')
                p.scatter(self.action.args.midcntr,
                          [self.action.args.midrow]*120, color='red')
                bokeh_plot(p)
                if self.context.config.plot_level >= 2:
                    input("next: ")
                else:
                    time.sleep(self.context.config.instrument.plot_pause)
            trace = {
                'src': src,
                'dst': dst,
                'barid': barid,
                'slid': slid,
                'MIDROW': self.action.args.midrow,
                'WINDOW': self.action.args.win,
                'REFDELX': self.action.args.refdelx,
                'CBARSNO': self.action.args.cbarsno,
                'CBARSFL': self.action.args.cbarsfl}

            # in this line we pass the trace information to an argument
            # instead of writing it to a table
            self.context.trace = trace
            ofname = self.action.args.cbarsfl.split('.')[0] + "_trace.fits"
            write_table(table=[src, dst, barid, slid],
                        names=('src', 'dst', 'barid', 'slid'),
                        output_dir=os.path.dirname(self.action.args.name),
                        output_name=ofname,
                        comment=['Source and destination fiducial points',
                                 'Derived from KCWI continuum bars images',
                                 'For defining spatial transformation'],
                        keywords={'MIDROW': (self.action.args.midrow,
                                             "Middle Row of image"),
                                  'WINDOW': (self.action.args.win,
                                             "Window for bar"),
                                  'REFDELX': (self.action.args.refdelx,
                                              "Reference bar sep in px"),
                                  'CBARSNO': (self.action.args.cbarsno,
                                              "Cont. bars image number"),
                                  'CBARSFL': (self.action.args.cbarsfl,
                                              "Cont. bars image")})

            if self.context.config.instrument.saveintims:
                # fit transform
                self.logger.info("Fitting spatial control points")
                tform = tf.estimate_transform('polynomial', src, dst, order=3)
                self.logger.info("Transforming bars image")
                warped = tf.warp(self.action.args.ccddata.data, tform)
                # write out warped image
                self.action.args.ccddata.data = warped
                kcwi_fits_writer(self.action.args.ccddata,
                                 self.action.args.table,
                                 output_file=self.action.args.name,
                                 suffix='warped')
                self.logger.info("Transformed bars produced")
            return self.action.args
    # END: class TraceBars()


class ExtractArcs(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Extracting arc spectra")
        tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                             target_type='CONTBARS',
                                             nearest=True)
        self.logger.info("%d continuum bars frames found" % len(tab))
        # ofname = tab['OFNAME'][0]
        ofname = tab['OFNAME'][0].split('.')[0] + "_trace.fits"
        print("*************** READING TABLE: %s" % ofname)
        # trace = read_table(tab=tab, indir='redux', suffix='trace')
        # Find  and read control points from continuum bars
        if hasattr(self.context, 'trace'):
            trace = self.context.trace
        else:
            trace = read_table(input_dir=os.path.dirname(self.action.args.name),
                               file_name=ofname)
            self.context.trace = {}
            for key in trace.meta.keys():
                self.context.trace[key] = trace.meta[key]
        midrow = self.context.trace['MIDROW']
        win = self.context.trace['WINDOW']
        self.action.args.refdelx = self.context.trace['REFDELX']
        self.action.args.cbarsno = self.context.trace['CBARSNO']
        self.action.args.cbarsfl = self.context.trace['CBARSFL']
        self.action.args.arcno = self.action.args.ccddata.header['FRAMENO']
        self.action.args.arcfl = self.action.args.ccddata.header['OFNAME']

        src = trace['src']  # source control points
        dst = trace['dst']  # destination control points

        self.logger.info("Fitting spatial control points")
        tform = tf.estimate_transform('polynomial', src, dst, order=3)

        self.logger.info("Transforming arc image")
        warped = tf.warp(self.action.args.ccddata.data, tform)
        # Write warped arcs if requested
        # if self.frame.saveintims():
        #    # write out warped image
        #    self.frame.data = warped
        #    self.write_image(suffix='warped')
        #    self.log.info("Transformed arcs produced")
        # extract arcs
        self.logger.info("Extracting arcs")
        arcs = []
        for xyi, xy in enumerate(src):
            if xy[1] == midrow:
                xi = int(xy[0]+0.5)
                arc = np.median(
                    warped[:, (xi - win):(xi + win + 1)], axis=1)
                arc = arc - np.nanmin(arc[100:-100])    # avoid ends
                arcs.append(arc)
        # Did we get the correct number of arcs?
        if len(arcs) == self.context.config.instrument.NBARS:
            self.logger.info("Extracted %d arcs" % len(arcs))
            self.context.arcs = arcs
        else:
            self.logger.error("Did not extract %d arcs, extracted %d" %
                              (self.context.config.instrument.NBARS, len(arcs)))
        return self.action.args
# END: class ExtractArcs()


class ArcOffsets(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Finding inter-bar offsets")
        arcs = self.context.arcs
        if arcs is not None:
            # Do we plot?
            if self.context.config.plot_level >= 2:
                do_plot = True
            else:
                do_plot = False
            # Compare with reference arc
            refarc = arcs[self.context.config.instrument.REFBAR][:]
            # number of cross-correlation samples (avoiding ends)
            nsamp = len(refarc[10:-10])
            # possible offsets
            offar = np.arange(1-nsamp, nsamp)
            # Collect offsets
            offsets = []
            for na, arc in enumerate(arcs):
                # Cross-correlate, avoiding junk on the ends
                xcorr = np.correlate(refarc[10:-10], arc[10:-10], mode='full')
                # Calculate offset
                offset = offar[xcorr.argmax()]
                offsets.append(offset)
                self.logger.info("Arc %d Slice %d XCorr shift = %d" %
                                 (na, int(na/5), offset))
                # display if requested
                if do_plot:
                    p = figure(title=self.action.args.plotlabel +
                               ", Arc %d Slice %d XCorr, Shift = %d" %
                               (na, int(na/5), offset),
                               x_axis_label="CCD y (px)", y_axis_label="e-",
                               plot_width=self.config.instrument.plot_width,
                               plot_height=self.config.instrument.plot_height)
                    x = range(len(refarc))
                    p.line(x, refarc, color='green', legend='ref bar (%d)' %
                           self.context.config.instrument.REFBAR)
                    p.line(x, np.roll(arc, offset), color='red',
                           legend='bar %d' % na)
                    bokeh_plot(p)
                    q = input("<cr> - Next, q to quit: ")
                    if 'Q' in q.upper():
                        do_plot = False
            self.context.baroffs = offsets
        else:
            self.logger.error("No extracted arcs found")
        return self.action.args
    # END: class ArcOffsets()


class CalcPrelimDisp(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        # get binning
        ybin = self.action.args.ybinsize
        # 0 - compute alpha
        prelim_alpha = self.action.args.grangle - 13.0 - self.action.args.adjang
        # 1 - compute preliminary angle of diffraction
        prelim_beta = self.action.args.camangle - prelim_alpha
        # 2 - compute preliminary dispersion
        prelim_disp = math.cos(prelim_beta/math.degrees(1.)) / \
            self.action.args.rho / self.context.config.instrument.FCAM * \
            (self.context.config.instrument.PIX*ybin) * 1.e4
        prelim_disp *= math.cos(
            self.context.config.instrument.GAMMA/math.degrees(1.))
        self.logger.info("Initial alpha, beta (deg): %.3f, %.3f" %
                         (prelim_alpha, prelim_beta))
        self.logger.info("Initial calculated dispersion (A/binned pix): %.3f" %
                         prelim_disp)
        self.context.prelim_disp = prelim_disp
        return self.action.args
    # END: class CalcPrelimDisp()


class ReadAtlas(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        # What lamp are we using?
        lamp = self.action.args.illum
        # rez factor
        # if 'fear' in lamp.lower():
        #    rezfact = 0.5
        # else:
        #    rezfact = 1.0
        # atpath = os.path.join("../data", "%s.fits" % lamp.lower())
        # Does the atlas file exist?
        path = "data/%s.fits" % lamp.lower()  # always use slash
        pkg = __name__.split('.')[0]
        atpath = pkg_resources.resource_filename(pkg, path)
        if os.path.exists(atpath):
            self.logger.info("Reading atlas spectrum in: %s" % atpath)
        else:
            self.logger.error("Atlas spectrum not found for %s" % atpath)
        # Read the atlas
        ff = pf.open(atpath)
        reflux = ff[0].data
        refdisp = ff[0].header['CDELT1']
        refwav = np.arange(0, len(reflux)) * refdisp + ff[0].header['CRVAL1']
        ff.close()
        # Convolve with appropriate Gaussian
        resolution = self.action.args.resolution
        atrespix = resolution / refdisp
        self.logger.info("Resolution = %.3f Ang, or %.2f Atlas px" %
                         (resolution, atrespix))
        reflux = gaussian_filter1d(reflux, atrespix/2.354)
        # Observed arc spectrum
        obsarc = self.context.arcs[self.context.config.instrument.REFBAR]
        # Preliminary wavelength solution
        xvals = np.arange(0, len(obsarc)) - int(len(obsarc)/2)
        obswav = xvals * self.context.prelim_disp + self.action.args.cwave
        # Get central third
        minow = int(len(obsarc)/3)
        maxow = int(2.*len(obsarc)/3)
        # Unless we are low dispersion, then get central 3 5ths
        if 'BL' in self.action.args.grating or 'RL' in self.action.args.grating:
            minow = int(len(obsarc)/5)
            maxow = int(4.*len(obsarc)/5)
        minwav = obswav[minow]
        maxwav = obswav[maxow]
        # Get corresponding ref range
        minrw = [i for i, v in enumerate(refwav) if v >= minwav][0]
        maxrw = [i for i, v in enumerate(refwav) if v <= maxwav][-1]
        # Subsample for cross-correlation
        cc_obsarc = obsarc[minow:maxow]
        cc_obswav = obswav[minow:maxow]
        cc_reflux = reflux[minrw:maxrw]
        cc_refwav = refwav[minrw:maxrw]
        # Resample onto reference wavelength scale
        obsint = interpolate.interp1d(cc_obswav, cc_obsarc, kind='cubic',
                                      bounds_error=False,
                                      fill_value='extrapolate'
                                      )
        cc_obsarc = obsint(cc_refwav)
        # Apply cosign bell taper to both
        cc_obsarc *= signal.windows.tukey(
            len(cc_obsarc), alpha=self.config.instrument.TAPERFRAC)
        cc_reflux *= signal.windows.tukey(
            len(cc_reflux), alpha=self.config.instrument.TAPERFRAC)
        nsamp = len(cc_refwav)
        offar = np.arange(1 - nsamp, nsamp)
        # Cross-correlate
        xcorr = np.correlate(cc_obsarc, cc_reflux, mode='full')
        # Get central region
        x0c = int(len(xcorr)/3)
        x1c = int(2*(len(xcorr)/3))
        xcorr_central = xcorr[x0c:x1c]
        offar_central = offar[x0c:x1c]
        # Calculate offset
        offset_pix = offar_central[xcorr_central.argmax()]
        offset_wav = offset_pix * refdisp
        self.logger.info("Initial arc-atlas offset (px, Ang): %d, %.1f" %
                         (offset_pix, offset_wav))
        if self.context.config.plot_level >= 1:
            # Plot
            p = figure(title=self.action.args.plotlabel +
                       ", (%s), Offset = %d px" % (lamp, offset_pix),
                       x_axis_label="Offset(px)", y_axis_label="X-corr",
                       plot_width=self.context.config.instrument.plot_width,
                       plot_height=self.context.config.instrument.plot_height)

            p.line(offar_central, xcorr_central, legend='Data')
            ylim_min = min(xcorr_central)
            ylim_max = max(xcorr_central)
            p.line([offset_pix, offset_pix], [ylim_min, ylim_max],
                   color='green', legend='Peak')
            bokeh_plot(p)
            if self.context.config.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.context.config.instrument.plot_pause)
            # Get central wavelength
            cwave = self.action.args.cwave
            # Set up offset tweaking
            q = 'test'
            while q:
                # Plot the two spectra
                p = figure(title=self.action.args.plotlabel +
                           ", (%s), Offset = %.1f Ang (%d px)" % (lamp,
                                                                  offset_wav,
                                                                  offset_pix),
                           x_axis_label="Wave(A)", y_axis_label="Rel. Flux",
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.line(obswav[minow:maxow] - offset_wav,
                       obsarc[minow:maxow]/np.nanmax(obsarc[minow:maxow]),
                       legend="ref bar (%d)" %
                              self.context.config.instrument.REFBAR)
                p.line(refwav[minrw:maxrw],
                       reflux[minrw:maxrw]/np.nanmax(reflux[minrw:maxrw]),
                       color="red", legend="Atlas")
                p.x_range = Range1d(np.nanmin(obswav[minow:maxow]),
                                    np.nanmax(obswav[minow:maxow]))
                ylim_min = min(
                    obsarc[minow:maxow]/np.nanmax(obsarc[minow:maxow]))
                ylim_max = max(
                    obsarc[minow:maxow]/np.nanmax(obsarc[minow:maxow]))
                p.line([cwave, cwave], [ylim_min, ylim_max], color="green",
                       legend="CWAVE")
                bokeh_plot(p)

                if self.context.config.plot_level >= 2:
                    q = input("Enter: <cr> - next, new offset (px): ")
                    if q:
                        offset_pix = int(q)
                        offset_wav = offset_pix * refdisp
                else:
                    time.sleep(self.context.config.instrument.plot_pause)
                    q = None
            self.logger.info("Final   arc-atlas offset (px, Ang): %d, %.1f" %
                             (offset_pix, offset_wav))
        # Store atlas spectrum
        self.action.args.reflux = reflux
        self.action.args.refwave = refwav
        # Store offsets
        self.action.args.offset_pix = offset_pix
        self.action.args.offset_wave = offset_wav
        # Store reference dispersion
        self.action.args.refdisp = refdisp
        # Store central limits
        self.action.args.minrow = minow
        self.action.args.maxrow = maxow
        # Store x values
        self.action.args.xvals = xvals
        self.action.args.x0 = int(len(obsarc)/2)
        return self.action.args
    # END: class ReadAtlas()


def bar_fit_helper(argument):

    b = argument['b']
    bs = argument['bs']
    # wavelength coefficients
    coeff = [0., 0., 0., 0., 0.]
    # container for maxima, shifts
    maxima = []
    shifts = []
    # get sub spectrum for this bar
    subspec = bs[argument['minrow']:argument['maxrow']]
    # now loop over dispersions
    for di, disp in enumerate(argument['disps']):
        # populate the coefficients
        coeff[4] = argument['p0'][b]
        coeff[3] = disp
        cosbeta = disp / (argument['PIX'] * argument['ybin']) * \
            argument['rho'] * argument['FCAM'] * 1.e-4
        if cosbeta > 1.:
            cosbeta = 1.
        beta = math.acos(cosbeta)
        coeff[2] = -(argument['PIX'] * argument['ybin'] /
                     argument['FCAM']) ** 2 * math.sin(beta) / 2. / \
            argument['rho'] * 1.e4
        coeff[1] = -(argument['PIX'] * argument['ybin'] /
                     argument['FCAM']) ** 3 * math.cos(beta) / 6. / \
            argument['rho'] * 1.e4
        coeff[0] = (argument['PIX'] * argument['ybin'] /
                    argument['FCAM']) ** 4 * math.sin(beta) / 24. / \
            argument['rho'] * 1.e4
        # what are the min and max wavelengths to consider?
        wl0 = np.polyval(coeff, argument['xvals'][argument['minrow']])
        wl1 = np.polyval(coeff, argument['xvals'][argument['maxrow']])
        minwvl = np.nanmin([wl0, wl1])
        maxwvl = np.nanmax([wl0, wl1])
        # where will we need to interpolate to cross-correlate?
        minrw = [i for i, v in enumerate(argument['refwave'])
                 if v >= minwvl][0]
        maxrw = [i for i, v in enumerate(argument['refwave'])
                 if v <= maxwvl][-1]
        subrefwvl = argument['refwave'][minrw:maxrw]
        subrefspec = argument['reflux'][minrw:maxrw]
        # get bell cosine taper to avoid nasty edge effects
        tkwgt = signal.windows.tukey(len(subrefspec),
                                     alpha=argument['taperfrac'])
        # apply taper to atlas spectrum
        subrefspec *= tkwgt
        # adjust wavelengths
        waves = np.polyval(coeff, argument['subxvals'])
        # interpolate the bar spectrum
        obsint = interpolate.interp1d(waves, subspec, kind='cubic',
                                      bounds_error=False,
                                      fill_value='extrapolate')
        intspec = obsint(subrefwvl)
        # apply taper to bar spectrum
        intspec *= tkwgt
        # get a label
        # cross correlate the interpolated spectrum with the atlas spec
        nsamp = len(subrefwvl)
        offar = np.arange(1 - nsamp, nsamp)
        # Cross-correlate
        xcorr = np.correlate(intspec, subrefspec, mode='full')
        # Get central region
        x0c = int(len(xcorr) / 3)
        x1c = int(2 * (len(xcorr) / 3))
        xcorr_central = xcorr[x0c:x1c]
        offar_central = offar[x0c:x1c]
        # Calculate offset
        maxima.append(xcorr_central[xcorr_central.argmax()])
        shifts.append(offar_central[xcorr_central.argmax()])
    # Get interpolations
    int_max = interpolate.interp1d(argument['disps'], maxima, kind='cubic',
                                   bounds_error=False,
                                   fill_value='extrapolate')
    int_shift = interpolate.interp1d(argument['disps'], shifts, kind='cubic',
                                     bounds_error=False,
                                     fill_value='extrapolate')
    xdisps = np.linspace(min(argument['disps']), max(argument['disps']),
                         num=argument['nn'] * 100)
    # get peak values
    maxima_res = int_max(xdisps)
    shifts_res = int_shift(xdisps) * argument['refdisp']
    bardisp = xdisps[maxima_res.argmax()]
    barshift = shifts_res[maxima_res.argmax()]
    # update coeffs
    coeff[4] = argument['p0'][b] - barshift
    coeff[3] = bardisp
    cosbeta = coeff[3] / (argument['PIX'] * argument['ybin']) * \
        argument['rho'] * argument['FCAM'] * 1.e-4
    if cosbeta > 1.:
        cosbeta = 1.
    beta = math.acos(cosbeta)
    coeff[2] = -(argument['PIX'] * argument['ybin'] / argument['FCAM']) ** 2 * \
        math.sin(beta) / 2. / argument['rho'] * 1.e4
    coeff[1] = -(argument['PIX'] * argument['ybin'] / argument['FCAM']) ** 3 * \
        math.cos(beta) / 6. / argument['rho'] * 1.e4
    coeff[0] = (argument['PIX'] * argument['ybin'] / argument['FCAM']) ** 4 * \
        math.sin(beta) / 24. / argument['rho'] * 1.e4
    scoeff = pascal_shift(coeff, argument['x0'])
    print("Central Fit: Bar#, Cdisp, Coefs: "
          "%3d  %.4f  %.2f  %.4f  %13.5e %13.5e" %
          (b, bardisp, scoeff[4], scoeff[3], scoeff[2], scoeff[1]))
    # Return results
    return b, scoeff, coeff[4], coeff[3], maxima
    # END: def bar_fit_helper()


class FitCenter(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        self.action.args.twkcoeff = []

    def _perform(self):
        """ Fit central region

        At this point we have the offsets between bars and the approximate
        offset from the reference bar to the atlas spectrum and the approximate
        dispersion.
        """
        self.logger.info("Finding wavelength solution for central region")
        # Are we interactive?
        do_inter = (self.context.config.plot_level >= 2)

        # y binning
        ybin = self.action.args.ybinsize
        # let's populate the 0 points vector
        p0 = self.action.args.cwave + np.array(self.context.baroffs) * \
            self.context.prelim_disp - self.action.args.offset_wave
        # next we are going to brute-force scan around the preliminary
        # dispersion for a better solution. We will wander 5% away from it.
        max_ddisp = 0.05    # fraction
        # we will try nn values
        nn = (int(max_ddisp*abs(self.context.prelim_disp) /
                  self.action.args.refdisp * (self.action.args.maxrow -
                                              self.action.args.minrow)/3.0))
        if nn < 10:
            nn = 10
        if nn > 25:
            nn = 25
        self.logger.info("N disp. samples: %d" % nn)
        # dispersions to try
        disps = self.context.prelim_disp * (1.0 + max_ddisp *
                                            (np.arange(0, nn+1) - nn/2.) *
                                            2.0 / nn)
        # values for central fit
        subxvals = self.action.args.xvals[
                   self.action.args.minrow:self.action.args.maxrow]

        # loop over bars and assemble input arguments
        my_arguments = []
        for b, bs in enumerate(self.context.arcs):
            arguments = {
                'b': b, 'bs': bs, 'minrow': self.action.args.minrow,
                'maxrow': self.action.args.maxrow, 'disps': disps, 'p0': p0,
                'PIX': self.context.config.instrument.PIX, 'ybin': ybin,
                'rho': self.action.args.rho,
                'FCAM': self.context.config.instrument.FCAM,
                'xvals': self.action.args.xvals,
                'refwave': self.action.args.refwave,
                'reflux': self.action.args.reflux,
                'taperfrac': self.context.config.instrument.TAPERFRAC,
                'refdisp': self.action.args.refdisp, 'subxvals': subxvals,
                'nn': nn, 'x0': self.action.args.x0
            }
            my_arguments.append(arguments)

        twkcoeff = {}
        centwave = []
        centdisp = []

        p = Pool()
        results = p.map(bar_fit_helper, list(my_arguments))
        p.close()

        for result in results:
            b = result[0]
            scoeff = result[1]
            _centwave = result[2]
            _centdisp = result[3]
            twkcoeff[b] = scoeff
            centwave.append(_centwave)
            centdisp.append(_centdisp)
            maxima = result[4]
            if do_inter:
                # plot maxima
                # image label
                imlab = "Img # %d (%s) Sl: %s Fl: %s Gr: %s" % \
                        (self.action.args.ccddata.header['FRAMENO'],
                         self.action.args.illum,
                         self.action.args.ifuname, self.action.args.filter,
                         self.action.args.grating)
                p = figure(title=imlab+": Bar %d, Slice %d" % (b, int(b / 5)),
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height,
                           x_axis_label="Central dispersion (Ang/px)",
                           y_axis_label="X-Corr Peak Value")
                p.scatter(disps, maxima, color='red', legend="Data")
                p.line(disps, maxima, color='blue')
                ylim = [min(maxima), max(maxima)]
                p.line([_centdisp, _centdisp], ylim, color='green',
                       legend="Fit Disp")
                p.line([self.context.prelim_disp, self.context.prelim_disp],
                       ylim, color='red', legend="Calc Disp")
                bokeh_plot(p)
                q = input("<cr> - Next, q to quit: ")
                if 'Q' in q.upper():
                    do_inter = False

        self.action.args.twkcoeff = twkcoeff
        # Plot results
        if self.context.config.plot_level >= 1:
            # Plot central wavelength
            p = figure(title=self.action.args.plotlabel, x_axis_label="Bar #",
                       y_axis_label="Central Wavelength (A)",
                       plot_width=self.config.instrument.plot_width,
                       plot_height=self.config.instrument.plot_height)
            x = range(len(centwave))
            p.scatter(x, centwave, marker='x', legend='bar wave')
            p.line([0, 120], [self.action.args.cwave, self.action.args.cwave],
                   color='red', legend='CWAVE')
            ylim = [min(centwave), max(centwave)]
            for ix in range(1, 24):
                sx = ix*5 - 0.5
                p.line([sx, sx], ylim, color='black')
            p.x_range = Range1d(-1, 120)
            p.legend.location = "top_center"
            bokeh_plot(p)
            if self.context.config.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.context.config.instrument.plot_pause)
            # Plot central dispersion
            p = figure(title=self.action.args.plotlabel, x_axis_label="Bar #",
                       y_axis_label="Central Dispersion (A)",
                       plot_width=self.config.instrument.plot_width,
                       plot_height=self.config.instrument.plot_height)
            x = range(len(centdisp))
            p.scatter(x, centdisp, marker='x', legend='bar disp')
            p.line([0, 120], [self.context.prelim_disp,
                              self.context.prelim_disp], color='red',
                   legend='Calc Disp')
            ylim = [min(centdisp), max(centdisp)]
            for ix in range(1, 24):
                sx = ix * 5 - 0.5
                p.line([sx, sx], ylim,  color='black')
            p.x_range = Range1d(-1, 120)
            p.legend.location = "bottom_center"
            bokeh_plot(p)
            if self.context.config.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.context.config.instrument.plot_pause)

        return self.action.args
    # END: class FitCenter()


class GetAtlasLines(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        self.action.args.atminrow = None
        self.action.args.atmaxrow = None
        self.action.args.atminwave = None
        self.action.args.atmaxwave = None
        self.action.args.at_wave = None

    def _perform(self):
        self.logger.info("Finding isolated atlas lines")
        """Get relevant atlas line positions and wavelengths"""
        # get atlas wavelength range
        #
        # get pixel values (no longer centered in the middle)
        specsz = len(self.context.arcs[self.context.config.instrument.REFBAR])
        xvals = np.arange(0, specsz)
        # min, max rows
        minrow = 50
        maxrow = specsz - 50
        # wavelength range
        mnwvs = []
        mxwvs = []
        # Get wavelengths
        for b in range(self.context.config.instrument.NBARS):
            waves = np.polyval(self.action.args.twkcoeff[b], xvals)
            mnwvs.append(np.min(waves))
            mxwvs.append(np.max(waves))
        minwav = min(mnwvs) + 10.
        maxwav = max(mxwvs) - 10.
        # Get corresponding atlas range
        minrw = [i for i, v in enumerate(self.action.args.refwave)
                 if v >= minwav][0]
        maxrw = [i for i, v in enumerate(self.action.args.refwave)
                 if v <= maxwav][-1]
        self.logger.info("Min, Max wave (A): %.2f, %.2f" % (minwav, maxwav))
        # store atlas ranges
        self.action.args.atminrow = minrw
        self.action.args.atmaxrow = maxrw
        self.action.args.atminwave = minwav
        self.action.args.atmaxwave = maxwav
        # get atlas sub spectrum
        atspec = self.action.args.reflux[minrw:maxrw]
        atwave = self.action.args.refwave[minrw:maxrw]
        # get reference bar spectrum
        subxvals = xvals[minrow:maxrow]
        subyvals = self.context.arcs[self.context.config.instrument.REFBAR][
                   minrow:maxrow].copy()
        subwvals = np.polyval(
            self.action.args.twkcoeff[self.context.config.instrument.REFBAR],
            subxvals)
        # smooth subyvals
        win = boxcar(3)
        subyvals = sp.signal.convolve(subyvals, win, mode='same') / sum(win)
        # find good peaks in arc spectrum
        smooth_width = 4  # in pixels
        ampl_thresh = 0.
        # slope_thresh = 0.7 * smooth_width / 1000.   # more severe for arc
        # for fitting peaks
        peak_width = int(self.action.args.resolution /
                         self.action.args.refdisp)  # in pixels
        if peak_width < 4:
            peak_width = 4
        # slope_thresh = peak_width / 12000.
        slope_thresh = 0.016 / peak_width
        self.logger.info("Using a peak_width of %d px, a slope_thresh of %.5f "
                         "a smooth_width of %d and an ampl_thresh of %.3f" %
                         (peak_width, slope_thresh, smooth_width, ampl_thresh))
        init_cent, avwsg, init_hgt = findpeaks(atwave, atspec, smooth_width,
                                               slope_thresh, ampl_thresh,
                                               peak_width)
        avwfwhm = avwsg * 2.354
        self.logger.info("Found %d peaks with <sig> = %.3f (A),"
                         " <FWHM> = %.3f (A)" % (len(init_cent), avwsg,
                                                 avwfwhm))
        if 'BH' in self.action.args.grating or 'BM' in self.action.args.grating:
            fwid = avwfwhm
        else:
            fwid = avwfwhm
        # clean near neighbors
        diffs = np.diff(init_cent)
        spec_cent = []
        spec_hgt = []
        rej_neigh_w = []
        rej_neigh_y = []
        neigh_fact = 1.25
        for i, w in enumerate(init_cent):
            if i == 0:
                if diffs[i] < avwfwhm * neigh_fact:
                    rej_neigh_w.append(w)
                    rej_neigh_y.append(init_hgt[i])
                    continue
            elif i == len(diffs):
                if diffs[i - 1] < avwfwhm * neigh_fact:
                    rej_neigh_w.append(w)
                    rej_neigh_y.append(init_hgt[i])
                    continue
            else:
                if diffs[i - 1] < avwfwhm * neigh_fact or \
                        diffs[i] < avwfwhm * neigh_fact:
                    rej_neigh_w.append(w)
                    rej_neigh_y.append(init_hgt[i])
                    continue
            spec_cent.append(w)
            spec_hgt.append(init_hgt[i])
        self.logger.info("Found %d isolated peaks" % len(spec_cent))
        #
        # generate an atlas line list
        refws = []
        refas = []
        rej_fit_w = []
        rej_fit_y = []
        rej_par_w = []
        rej_par_a = []
        nrej = 0
        for i, pk in enumerate(spec_cent):
            # Fit Atlas Peak
            line_x = [ii for ii, v in enumerate(atwave) if v >= pk][0]
            minow, maxow, count = get_line_window(atspec, line_x)
            if count < 5 or not minow or not maxow:
                rej_fit_w.append(pk)
                rej_fit_y.append(spec_hgt[i])
                nrej += 1
                self.logger.info("Atlas window rejected for line %.3f" % pk)
                continue
            yvec = atspec[minow:maxow + 1]
            xvec = atwave[minow:maxow + 1]
            try:
                fit, _ = curve_fit(gaus, xvec, yvec, p0=[100., pk, 1.])
            except RuntimeError:
                rej_fit_w.append(pk)
                rej_fit_y.append(spec_hgt[i])
                nrej += 1
                self.logger.info("Atlas Gaussian fit rejected for line %.3f" %
                                 pk)
                continue
            int_line = interpolate.interp1d(xvec, yvec, kind='cubic',
                                            bounds_error=False,
                                            fill_value='extrapolate')
            x_dense = np.linspace(min(xvec), max(xvec), num=1000)
            # get peak value
            y_dense = int_line(x_dense)
            pki = y_dense.argmax()
            pkw = x_dense[pki]
            # pka = yvec.argmax()
            xoff = abs(pkw - fit[1]) / self.action.args.refdisp  # in pixels
            woff = abs(pkw - pk)  # in Angstroms
            wrat = abs(fit[2]) / fwid  # can be neg or pos
            if woff > 1. or xoff > 1. or wrat > 1.1:
                rej_par_w.append(pkw)
                rej_par_a.append(y_dense[pki])
                nrej += 1
                self.logger.info("Atlas line parameters rejected for line %.3f"
                                 % pk)
                self.logger.info("woff = %.3f, xoff = %.2f, wrat = %.3f" %
                                 (woff, xoff, wrat))
                continue
            refws.append(pkw)
            refas.append(y_dense[pki])
        # store wavelengths
        self.action.args.at_wave = refws
        # plot results
        # image label
        imlab = "Img # %d (%s) Sl: %s Fl: %s Gr: %s" % \
                (self.action.args.ccddata.header['FRAMENO'],
                 self.action.args.illum,
                 self.action.args.ifuname, self.action.args.filter,
                 self.action.args.grating)
        norm_fac = np.nanmax(atspec)
        p = figure(title=imlab + ": Ngood = %d, Nrej = %d" % (len(refws), nrej),
                   x_axis_label="Wavelength (A)",
                   y_axis_label="Normalized Flux",
                   plot_width=self.config.instrument.plot_width,
                   plot_height=self.config.instrument.plot_height)
        p.line(subwvals, subyvals / np.nanmax(subyvals), legend='RefArc',
               color='lightgray')
        p.line(atwave, atspec / norm_fac, legend='Atlas', color='blue')
        # Initial line list
        # p.diamond(init_cent, init_hgt / norm_fac, legend='Fpks',
        #          color='black', size=10, fill_alpha=0.1)
        # Rejected: nearby neighbor
        p.diamond(rej_neigh_w, rej_neigh_y / norm_fac, legend='NeighRej',
                  color='red', size=8)
        # Rejected: fit failure
        p.diamond(rej_fit_w, rej_fit_y / norm_fac, legend='FitRej',
                  color='cyan', size=8)
        # Rejected: line parameter outside range
        p.diamond(rej_par_w, rej_par_a / norm_fac, legend='ParRej',
                  color='orange', size=8)
        p.diamond(refws, refas / norm_fac, legend='Kept', color='green',
                  size=10)
        p.x_range = Range1d(min(subwvals), max(subwvals))
        if self.context.config.plot_level >= 1:
            bokeh_plot(p)
            if self.context.config.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                pl.pause(self.context.config.instrument.plot_pause)
        export_png(p, "atlas_lines_%s_%s_%s_%05d.png" %
                   (self.action.args.illum, self.action.args.grating,
                    self.action.args.ifuname,
                    self.action.args.ccddata.header['FRAMENO']))
        self.logger.info("Final atlas list has %d lines" % len(refws))

        return self.action.args
    # END: class GetAtlasLines()


class SolveArcs(BasePrimitive):

    def __init__(self, action, context):
        self.logger.info("Solving individual arc spectra (not yet implemented)")
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        self.fincoeff = []

    def _perform(self):
        """Solve the bar arc wavelengths"""
        if self.context.config.plot_level >= 2:
            master_inter = True
        else:
            master_inter = False
        if self.context.config.plot_level >= 3:
            do_inter = True
        else:
            do_inter = False
        verbose = False
        # Bar statistics
        bar_sig = []
        bar_nls = []
        # set thresh for finding lines
        hgt = 50.
        self.logger.info("line thresh = %.2f" % hgt)
        # get relevant part of atlas spectrum
        atwave = self.action.args.refwave[self.action.args.atminrow:
                                          self.action.args.atmaxrow]
        atspec = self.action.args.reflux[self.action.args.atminrow:
                                         self.action.args.atmaxrow]
        # get x values starting at zero pixels
        self.xsvals = np.arange(0, len(
            self.context.arcs[self.context.config.instrument.REFBAR]))
        # loop over arcs and generate a wavelength solution for each
        for ib, b in enumerate(self.context.arcs):
            # print("")
            # self.logger.info("BAR %d" % ib)
            # print("")
            # Starting pascal shifted coeffs from fit_center()
            coeff = self.action.args.twkcoeff[ib]
            # get bar wavelengths
            bw = np.polyval(coeff, self.xsvals)
            # smooth spectrum according to slicer
            if 'Small' in self.action.args.ifuname():
                bspec = b
            else:
                if 'Large' in self.action.args.ifuname():
                    win = boxcar(5)
                else:
                    win = boxcar(3)
                bspec = sp.signal.convolve(b, win, mode='same') / sum(win)
            # spmode = mode(np.round(bspec))
            # spmed = np.nanmedian(bspec)
            # self.logger.info("Arc spec median = %.3f, mode = %d" %
            #              (spmed, spmode.mode[0]))
            # store values to fit
            at_wave_dat = []
            arc_pix_dat = []
            rej_wave = []
            nrej = 0
            # loop over lines
            for iw, aw in enumerate(self.action.args.at_wave):
                # get window for this line
                try:
                    line_x = [i for i, v in enumerate(bw) if v >= aw][0]
                    minow, maxow, count = get_line_window(bspec, line_x,
                                                          thresh=hgt)
                    if count < 5 or not minow or not maxow:
                        rej_wave.append(aw)
                        nrej += 1
                        if verbose:
                            self.logger.info("Arc window rejected for line %.3f"
                                             % aw)
                        continue
                    # check if window no longer contains initial value
                    if minow > line_x > maxow:
                        rej_wave.append(aw)
                        nrej += 1
                        if verbose:
                            self.logger.info(
                                "Arc window wandered off for line %.3f" % aw)
                        continue
                    yvec = bspec[minow:maxow + 1]
                    xvec = self.xsvals[minow:maxow + 1]
                    wvec = bw[minow:maxow + 1]
                    max_value = yvec[yvec.argmax()]
                    # Gaussian fit
                    try:
                        fit, _ = curve_fit(gaus, xvec, yvec,
                                           p0=[100., line_x, 1.])
                    except RuntimeError:
                        nrej += 1
                        if verbose:
                            self.logger.info(
                                "Arc Gaussian fit rejected for line %.3f" % aw)
                        continue
                    sp_pk_x = fit[1]
                    # Get interpolation
                    int_line = interpolate.interp1d(xvec, yvec, kind='cubic',
                                                    bounds_error=False,
                                                    fill_value='extrapolate')
                    xplot = np.linspace(min(xvec), max(xvec), num=1000)
                    # get peak value
                    plt_line = int_line(xplot)
                    max_index = plt_line.argmax()
                    peak = xplot[max_index]
                    # Calculate centroid
                    cent = np.sum(xvec * yvec) / np.sum(yvec)
                    if abs(cent - peak) > 0.7:
                        rej_wave.append(aw)
                        nrej += 1
                        if verbose:
                            self.logger.info("Arc peak - cent offset rejected "
                                             "for line %.3f" % aw)
                        continue
                    # store data
                    arc_pix_dat.append(peak)
                    at_wave_dat.append(aw)
                    # plot, if requested
                    if do_inter:
                        ptitle = "Bar: %d - %3d/%3d: x0, x1, Cent, Wave = " \
                                 "%d, %d, %8.1f, %9.2f" % \
                                 (ib, (iw + 1), len(self.action.args.at_wave),
                                  minow, maxow, cent, aw)
                        atx0 = [i for i, v in enumerate(atwave)
                                if v >= min(wvec)][0]
                        atx1 = [i for i, v in enumerate(atwave)
                                if v >= max(wvec)][0]
                        atnorm = np.nanmax(yvec) / np.nanmax(atspec[atx0:atx1])
                        p = figure(
                            title=ptitle, x_axis_label="Wavelength (A)",
                            y_axis_label="Relative Flux",
                            plot_width=self.config.instrument.plot_width,
                            plot_height=self.config.instrument.plot_height)
                        p.line(wvec, yvec, legend='Arc', color='black')
                        p.circle(wvec, yvec, color='red')
                        ylim = [0, np.nanmax(yvec)]
                        p.circle(atwave[atx0:atx1], atspec[atx0:atx1] * atnorm,
                                 color='green', legend='Atlas')
                        p.line([aw, aw], ylim, color='red', legend='W in')
                        input("next - <cr>: ")
                        pl.clf()
                        pl.plot(xvec, yvec, 'r.', label='Data')
                        pl.plot(xplot, plt_line, label='Interp')
                        ylim = [0, pl.gca().get_ylim()[1]]
                        xlim = pl.gca().get_xlim()
                        pl.plot(xlim, [max_value * 0.5, max_value * 0.5], 'k--')
                        pl.plot([cent, cent], ylim, 'g--', label='Cntr')
                        pl.plot([line_x, line_x], ylim, 'r-.', label='X in')
                        pl.plot([peak, peak], ylim, 'c-.', label='Peak')
                        pl.plot([sp_pk_x, sp_pk_x], ylim, 'm-.', label='Gpeak')
                        pl.xlabel("CCD Y (px)")
                        pl.ylabel("Flux (DN)")
                        pl.ylim(ylim)
                        pl.title(ptitle)
                        pl.legend()

                        q = input(ptitle + "; <cr> - Next, q to quit: ")
                        if 'Q' in q.upper():
                            do_inter = False
                            pl.ioff()
                except IndexError:
                    if verbose:
                        self.logger.info(
                            "Atlas line not in observation: %.2f" % aw)
                    rej_wave.append(aw)
                    nrej += 1
                    continue
                except ValueError:
                    if verbose:
                        self.logger.info(
                            "Interpolation error for line at %.2f" % aw)
                    rej_wave.append(aw)
                    nrej += 1
            self.logger.info("Fitting wavelength solution starting with %d "
                             "lines after rejecting %d lines" %
                             (len(arc_pix_dat), nrej))
            # Fit wavelengths
            # Initial fit
            wfit = np.polyfit(arc_pix_dat, at_wave_dat, 4)
            pwfit = np.poly1d(wfit)
            arc_wave_fit = pwfit(arc_pix_dat)
            resid = arc_wave_fit - at_wave_dat
            resid_c, low, upp = sigmaclip(resid, low=3., high=3.)
            wsig = resid_c.std()
            rej_rsd = []
            rej_rsd_wave = []
            # Iteratively remove outliers
            for it in range(4):
                # self.logger.info("Iteration %d" % it)
                arc_dat = []
                at_dat = []
                # Trim outliers
                for il, rsd in enumerate(resid):
                    if low < rsd < upp:
                        arc_dat.append(arc_pix_dat[il])
                        at_dat.append(at_wave_dat[il])
                    else:
                        # self.logger.info("REJ: %d, %.2f, %.3f" %
                        #              (il, arc_pix_dat[il], at_wave_dat[il]))
                        rej_rsd_wave.append(at_wave_dat[il])
                        rej_rsd.append(rsd)
                # refit
                arc_pix_dat = arc_dat.copy()
                at_wave_dat = at_dat.copy()
                wfit = np.polyfit(arc_pix_dat, at_wave_dat, 4)
                pwfit = np.poly1d(wfit)
                arc_wave_fit = pwfit(arc_pix_dat)
                resid = arc_wave_fit - at_wave_dat
                resid_c, low, upp = sigmaclip(resid, low=3., high=3.)
                wsig = np.nanstd(resid)
            # store results
            # print("")
            self.logger.info("Bar %03d, Slice = %02d, RMS = %.3f, N = %d" %
                             (ib, int(ib / 5), wsig, len(arc_pix_dat)))
            # print("")
            self.fincoeff.append(wfit)
            bar_sig.append(wsig)
            bar_nls.append(len(arc_pix_dat))
            # plot bar fit residuals
            if master_inter:
                pl.ion()
                pl.clf()
                pl.plot(at_wave_dat, resid, 'd', label='Rsd')
                ylim = pl.gca().get_ylim()
                if rej_rsd_wave:
                    pl.plot(rej_rsd_wave, rej_rsd, 'rd', label='Rej')
                pl.xlabel("Wavelength (A)")
                pl.ylabel("Fit - Inp (A)")
                pl.title(self.action.args.plotlabel() +
                         " Bar = %03d, Slice = %02d, RMS = %.3f, N = %d" %
                         (ib, int(ib / 5), wsig, len(arc_pix_dat)))
                xlim = [self.action.args.atminwave, self.action.args.atmaxwave]
                pl.plot(xlim, [0., 0.], '-')
                pl.plot(xlim, [wsig, wsig], '-.', color='gray')
                pl.plot(xlim, [-wsig, -wsig], '-.', color='gray')
                pl.plot([self.action.args.cwave(), self.action.args.cwave()],
                        ylim, '-.', label='CWAV')
                pl.xlim(xlim)
                pl.ylim(ylim)
                pl.legend()
                input("Next? <cr>: ")

                # overplot atlas and bar using fit wavelengths
                pl.clf()
                bwav = pwfit(self.xsvals)
                pl.plot(bwav, b, label='Arc')
                ylim = pl.gca().get_ylim()
                atnorm = np.nanmax(b) / np.nanmax(atspec)
                pl.plot(atwave, atspec * atnorm, label='Atlas')
                pl.plot([self.action.args.cwave(), self.action.args.cwave()],
                        ylim, 'm-.', label='CWAV')
                pl.xlim(xlim)
                pl.ylim(ylim)
                pl.xlabel("Wavelength (A)")
                pl.ylabel("Flux")
                pl.title(self.action.args.plotlabel() +
                         " Bar = %03d, Slice = %02d, RMS = %.3f, N = %d" %
                         (ib, int(ib / 5), wsig, len(arc_pix_dat)))
                leg_first = True
                for w in self.action.args.at_wave:
                    if leg_first:
                        pl.plot([w, w], ylim, 'k-.', label='Orig')
                        leg_first = False
                    else:
                        pl.plot([w, w], ylim, 'k-.')
                leg_first = True
                for w in arc_wave_fit:
                    if leg_first:
                        pl.plot([w, w], ylim, 'c-.', label='Kept')
                        leg_first = False
                    else:
                        pl.plot([w, w], ylim, 'c-.')
                leg_first = True
                for w in rej_rsd_wave:
                    if leg_first:
                        pl.plot([w, w], ylim, 'r-.', label='RejRsd')
                        leg_first = False
                    else:
                        pl.plot([w, w], ylim, 'r-.')
                leg_first = True
                for w in rej_wave:
                    if leg_first:
                        pl.plot([w, w], ylim, 'y-.', label='RejFit')
                        leg_first = False
                    else:
                        pl.plot([w, w], ylim, 'y-.')
                pl.legend()
                q = input("Next? <cr>, q - quit: ")
                if 'Q' in q.upper():
                    master_inter = False
                    pl.ioff()
        # Plot final results
        # Plot fit sigmas
        pl.clf()
        pl.plot(bar_sig, 'd', label='RMS')
        xlim = [-1, 120]
        ylim = pl.gca().get_ylim()
        self.av_bar_sig = float(np.nanmean(bar_sig))
        self.st_bar_sig = float(np.nanstd(bar_sig))
        self.logger.info("<STD>     = %.3f +- %.3f (A)" % (self.av_bar_sig,
                                                           self.st_bar_sig))
        pl.plot(xlim, [self.av_bar_sig, self.av_bar_sig], 'k--')
        pl.plot(xlim, [(self.av_bar_sig - self.st_bar_sig),
                       (self.av_bar_sig - self.st_bar_sig)], 'k:')
        pl.plot(xlim, [(self.av_bar_sig + self.st_bar_sig),
                       (self.av_bar_sig + self.st_bar_sig)], 'k:')
        for ix in range(1, 24):
            sx = ix * 5 - 0.5
            pl.plot([sx, sx], ylim, '-.', color='black')
        pl.xlabel("Bar #")
        pl.ylabel("RMS (A)")
        pl.title(self.action.args.plotlabel() +
                 " <RMS>: %.3f +- %.3f" % (self.av_bar_sig, self.st_bar_sig))
        pl.xlim(xlim)
        pl.gca().margins(0)
        pl.savefig("arc_%05d_resid_%s_%s_%s.png" %
                   (self.action.args.ccddata.header['FRAMENO'],
                    self.action.args.illum(),
                    self.action.args.grating(), self.action.args.ifuname()))
        if self.context.config.plot_level >= 1:
            bokeh_plot(p)
            if self.context.config.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                pl.pause(self.context.config.instrument.plot_pause)
        # Plot number of lines fit
        pl.clf()
        pl.plot(bar_nls, 'd', label='N lines')
        ylim = pl.gca().get_ylim()
        self.av_bar_nls = float(np.nanmean(bar_nls))
        self.st_bar_nls = float(np.nanstd(bar_nls))
        self.logger.info("<N Lines> = %.1f +- %.1f" % (self.av_bar_nls,
                                                       self.st_bar_nls))
        pl.plot(xlim, [self.av_bar_nls, self.av_bar_nls], 'k--')
        pl.plot(xlim, [(self.av_bar_nls - self.st_bar_nls),
                       (self.av_bar_nls - self.st_bar_nls)], 'k:')
        pl.plot(xlim, [(self.av_bar_nls + self.st_bar_nls),
                       (self.av_bar_nls + self.st_bar_nls)], 'k:')
        for ix in range(1, 24):
            sx = ix * 5 - 0.5
            pl.plot([sx, sx], ylim, '-.', color='black')
        pl.xlabel("Bar #")
        pl.ylabel("N Lines")
        pl.title(self.action.args.plotlabel() + " <N Lines>: %.3f +- %.3f" %
                 (self.av_bar_nls, self.st_bar_nls))
        pl.xlim(xlim)
        pl.gca().margins(0)
        pl.savefig("arc_%05d_nlines_%s_%s_%s.png" %
                   (self.action.args.ccddata.header['FRAMENO'],
                    self.action.args.illum(),
                    self.action.args.grating(), self.action.args.ifuname()))
        if self.context.config.plot_level >= 1:
            bokeh_plot(p)
            if self.context.config.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                pl.pause(self.context.config.instrument.plot_pause)
        # Plot coefs
        ylabs = ['Ang/px^4', 'Ang/px^3', 'Ang/px^2', 'Ang/px', 'Ang']
        for ic in reversed(range(len(self.fincoeff[0]))):
            pl.clf()
            coef = []
            for c in self.fincoeff:
                coef.append(c[ic])
            pl.plot(coef, 'd')
            ylim = pl.gca().get_ylim()
            for ix in range(1, 24):
                sx = ix * 5 - 0.5
                pl.plot([sx, sx], ylim, '-.', color='black')
            pl.xlabel("Bar #")
            pl.ylabel("Coef %d (%s)" % (ic, ylabs[ic]))
            pl.title(self.action.args.plotlabel() + " Coef %d" % ic)
            pl.xlim(xlim)
            pl.gca().margins(0)
            pl.savefig("arc_%05d_coef%d_%s_%s_%s.png" %
                       (self.action.args.ccddata.header['FRAMENO'], ic,
                        self.action.args.illum(), self.action.args.grating(),
                        self.action.args.ifuname()))
            if self.context.config.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                pl.pause(self.context.config.instrument.plot_pause)

        return self.action.args
    # END: class SolveArcs()


class SolveGeom(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Solving overall geometry (not yet implemented)")
        return self.action.args
    # END: class SolveGeom()


class GenerateMaps(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Generating geometry maps (not yet implemented)")
        return self.action.args
    # END: class GenerateMaps()


class ApplyFlat(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Applying flat field (not yet implemented)")
        return self.action.args
    # END: class ApplyFlat()


class SubtractSky(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Subtracting sky background (not yet implemented)")
        return self.action.args
    # END: class SubtractSky()


class MakeCube(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Creating data cube (not yet implemented)")
        return self.action.args
    # END: class MakeCube()


class CorrectDar(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Correcting for DAR (not yet implemented)")
        return self.action.args
    # END: class CorrectDar()


class FluxCalibrate(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Calibrating object flux (not yet implemented)")
        return self.action.args
    # END: class FluxCalibrate()


class MakeInvsens(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger0

    def _perform(self):
        self.logger.info("Making inverse sensitivity curve "
                         "(not yet implemented)")
        return self.action.args
    # END: class MakeInvsens()
