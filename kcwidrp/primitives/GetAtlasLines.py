from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_plotting import save_plot
from kcwidrp.primitives.kcwi_file_primitives import plotlabel

from bokeh.plotting import figure
from bokeh.models import Range1d
import numpy as np
import scipy as sp
from scipy.interpolate import interpolate
from scipy.signal.windows import boxcar
from scipy.optimize import curve_fit
from scipy.stats import sigmaclip
import time
import os
import traceback


def gaus(x, a, mu, sigma):
    """
    Gaussian fitting function.

    Args:
        x (float): x value
        a (float): amplitude
        mu (float): x position of center of Gaussian
        sigma (float): Gaussian width

    Returns:
        Gaussian value at x position

    """
    return a * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))


def get_line_window(y, c, thresh=0., logger=None, strict=False, maxwin=100,
                    frac_max=0.5):
    """
    Find a window that includes the specified amount of the line.

    Get a useful window on a given line to allow finding an accurate peak for
    each line.  Start at input center and iterate above and below until the
    window encompasses the specified fraction of the line maximum flux.
    Rejects lines that are problematic.

    Args:
        y (list of floats): flux values
        c (float): the nominal center of the line in pixels
        thresh (float): threshhold for rejection (default: 0.)
        logger (log object): logging instance for output
        strict (bool): should we apply strict criterion? (default: False)
        maxwin (int): largest possible window in pixels (default: 100)
        frac_max (float): fraction of maximum flux to encompass (default: 0.5)

    :returns:
        - x0 (int): lower limit of window in pixels
        - x1 (int): upper limit of window in pixels
        - count (int): number of pixels in window

    """
    verbose = logger is not None
    nx = len(y)
    # check edges
    if c < 2 or c > nx - 2:
        if verbose:
            logger.info("input center too close to edge")
        return None, None, 0
    # get initial values
    x0 = c - 2
    x1 = c + 2
    mx = np.nanmax(y[x0:x1+1])
    count = 5
    # check low side
    if x0 - 1 < 0:
        if verbose:
            logger.info("max check: low edge hit")
        return None, None, 0
    while y[x0-1] > mx:
        x0 -= 1
        count += 1
        if x0 - 1 < 0:
            if verbose:
                logger.info("Max check: low edge hit")
            return None, None, 0

    # check high side
    if x1 + 1 >= nx:
        if verbose:
            logger.info("max check: high edge hit")
        return None, None, 0
    while y[x1+1] > mx:
        x1 += 1
        count += 1
        if x1 + 1 >= nx:
            if verbose:
                logger.info("Max check: high edge hit")
            return None, None, 0
    # how big is our window?
    if (x1 - x0) > maxwin:
        if verbose:
            logger.info("Window expanded beyond limit")
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
    # expand until we get to max_frac
    fmx = mx * frac_max
    #
    # Low index side
    prev = mx
    while y[x0] > fmx:
        if y[x0] > mx or x0 <= 0 or y[x0] > prev:
            if verbose:
                if y[x0] > mx:
                    logger.info("frac_max check: low index err - missed max")
                if x0 <= 0:
                    logger.info("frac_max check: low index err - at edge")
                if y[x0] > prev:
                    logger.info("frac_max check: low index err - wiggly")
            return None, None, 0
        prev = y[x0]
        x0 -= 1
        count += 1
    # High index side
    prev = mx
    while y[x1] > fmx:
        if y[x1] > mx or x1 >= nx or y[x1] > prev:
            if verbose:
                if y[x1] > mx:
                    logger.info("frac_max check: high index err - missed max")
                if x1 >= nx:
                    logger.info("frac_max check: high index err - at edge")
                if y[x1] > prev:
                    logger.info("frac_max check: high index err - wiggly")
            return None, None, 0
        prev = y[x1]
        if x1 < (nx-1):
            x1 += 1
            count += 1
        else:
            if verbose:
                logger.info("Edge encountered")
            return None, None, 0
    if strict:
        # where did we end up?
        if c < x0 or x1 < c:
            if verbose:
                logger.info("initial position outside final window")
            return None, None, 0

    return x0, x1, count
    # END: get_line_window()


def findpeaks(x, y, wid, sth, ath, pkg=None, verbose=False):
    """
    Find peaks in spectrum

    Uses a gradient to determine where peaks are and then performs some
    cleaning of blended lines using sigma rejection.

    Args:
        x (list of floats): x values of vector
        y (list of floats): y values of vector
        wid (int): boxcar smoothing width in pixels
        sth (float): slope threshhold
        ath (float): amplitude threshhold
        pkg (float): limit on peak wandering in pixels
        verbose (bool): control output messages

    :returns:
        - list of peaks
        - average sigma of cleaned peaks
        - list of y values at peaks

    """
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
                            # gaussian fit
                            res, _ = curve_fit(gaus, xx, yy,
                                               p0=[y[i], x[i], 1.])
                            # check offset of fit from initial peak
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
    sgmn = None
    if len(pks) > 0:
        cln_sgs, low, upp = sigmaclip(sgs, low=3., high=3.)
        for i in range(len(pks)):
            # clean only blends (overly wide lines)
            if sgs[i] < upp:
                cpks.append(pks[i])
                cvals.append(hgt[i])
        sgmn = cln_sgs.mean()
        # sgmd = float(np.nanmedian(cln_sgs))
    else:
        print("No peaks found!")
    return cpks, sgmn, cvals
    # END: findpeaks()


class GetAtlasLines(BasePrimitive):
    """
    Get relevant atlas line positions and wavelengths

    Generates a list of well observed atlas lines for generating the final
    wavelength solution.

    Uses the following configuration parameters:

    * FRACMAX: fraction of max to use for window for fitting (defaults to 0.5)
    * LINELIST: an input line list to use instead of determining one on the fly

    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        self.action.args.atminrow = None
        self.action.args.atmaxrow = None
        self.action.args.atminwave = None
        self.action.args.atmaxwave = None
        self.action.args.at_wave = None
        self.action.args.at_flux = None

    def _pre_condition(self):
        self.logger.info("Checking for master arc")
        if 'MARC' in self.action.args.ccddata.header['IMTYPE']:
            return True
        else:
            return False

    def _perform(self):
        """Get atlas line positions for wavelength fitting"""
        self.logger.info("Finding isolated atlas lines")
        verbose = (self.config.instrument.verbose > 1)
        frac_max = self.config.instrument.FRACMAX
        self.logger.info("Finding line windows using fraction of line max of "
                         "%.2f" % frac_max)

        # get atlas wavelength range
        # get pixel values (no longer centered in the middle)
        specsz = len(self.context.arcs[self.config.instrument.REFBAR])
        xvals = np.arange(0, specsz)
        # min, max rows, trimming the ends
        minrow = 50
        maxrow = specsz - 50
        # wavelength range
        mnwvs = []
        mxwvs = []
        refbar_disp = 1.
        # Get wavelengths for each bar
        for b in range(self.config.instrument.NBARS):
            waves = np.polyval(self.action.args.twkcoeff[b], xvals)
            mnwvs.append(np.min(waves))
            mxwvs.append(np.max(waves))
            if b == self.config.instrument.REFBAR:
                refbar_disp = self.action.args.twkcoeff[b][-2]
        self.logger.info("Ref bar (%d) dispersion = %.3f Ang/px" %
                         (self.config.instrument.REFBAR, refbar_disp))
        # Get extrema (trim ends a bit)
        minwav = min(mnwvs) + 10.
        maxwav = max(mxwvs) - 10.
        wave_range = maxwav - minwav
        # Do we have a dichroic?
        if self.action.args.dich:
            if self.action.args.camera == 0:  # Blue
                maxwav = min([maxwav, 5620.])
            elif self.action.args.camera == 1:  # Red
                minwav = max([minwav, 5580.])
            else:
                self.logger.error("Camera keyword not defined!")
        dichroic_fraction = (maxwav - minwav) / wave_range
        # Get corresponding atlas range
        minrw = [i for i, v in enumerate(self.action.args.refwave)
                 if v >= minwav][0]
        maxrw = [i for i, v in enumerate(self.action.args.refwave)
                 if v <= maxwav][-1]
        self.logger.info("Min, Max wave (A): %.2f, %.2f" % (minwav, maxwav))
        if self.action.args.dich:
            self.logger.info("Dichroic fraction: %.3f" % dichroic_fraction)
        # store atlas ranges
        self.action.args.atminrow = minrw
        self.action.args.atmaxrow = maxrw
        self.action.args.atminwave = minwav
        self.action.args.atmaxwave = maxwav
        self.action.args.dichroic_fraction = dichroic_fraction
        # output filename stub
        atfnam = "arc_%05d_%s_%s_%s_atlines" % \
                 (self.action.args.ccddata.header['FRAMENO'],
                  self.action.args.illum, self.action.args.grating,
                  self.action.args.ifuname)
        # check if line list was given on command line
        if self.config.instrument.LINELIST:
            with open(self.config.instrument.LINELIST) as llfn:
                atlines = llfn.readlines()
            refws = []
            refas = []
            for line in atlines:
                if '#' in line:
                    continue
                refws.append(float(line.split()[0]))
                refas.append(float(line.split()[1]))
            self.logger.info("Read %d lines from %s" %
                             (len(refws), self.config.instrument.LINELIST))
        else:
            # get atlas sub spectrum
            atspec = self.action.args.reflux[minrw:maxrw]
            atwave = self.action.args.refwave[minrw:maxrw]
            # get ref bar arc spectrum, pixel values, and prelim wavelengths
            subxvals = xvals[minrow:maxrow]
            subyvals = self.context.arcs[self.config.instrument.REFBAR][
                       minrow:maxrow].copy()
            subwvals = np.polyval(
                self.action.args.twkcoeff[self.config.instrument.REFBAR],
                subxvals)
            # smooth subyvals
            win = boxcar(3)
            subyvals = sp.signal.convolve(subyvals, win, mode='same') / sum(win)
            # find good peaks in arc spectrum
            smooth_width = 4  # in pixels
            # peak width
            peak_width = int(self.action.args.atsig/abs(refbar_disp))
            if peak_width < 4:
                peak_width = 4
            # slope threshold
            slope_thresh = 0.7 * smooth_width / 2. / 100.
            # slope_thresh = 0.7 * smooth_width / 1000.   # more severe for arc
            # slope_thresh = 0.016 / peak_width
            # get amplitude threshold
            ampl_thresh = 0.
            self.logger.info("Using a peak_width of %d px, a slope_thresh of "
                             "%.5f a smooth_width of %d and an ampl_thresh of "
                             "%.3f" % (peak_width, slope_thresh, smooth_width,
                                       ampl_thresh))
            arc_cent, avwsg, arc_hgt = findpeaks(subwvals, subyvals,
                                                 smooth_width, slope_thresh,
                                                 ampl_thresh,
                                                 peak_width)
            avwfwhm = avwsg * 2.354
            self.logger.info("Found %d lines with <sig> = %.3f (A),"
                             " <FWHM> = %.3f (A)" % (len(arc_cent), avwsg,
                                                     avwfwhm))
            # fitting window based on grating type
            if 'H' in self.action.args.grating or \
                    'M' in self.action.args.grating:
                fwid = avwfwhm
            else:
                fwid = avwsg
            # clean near neighbors
            spec_cent = arc_cent
            spec_hgt = arc_hgt
            #
            # generate an atlas line list
            refws = []      # atlas line wavelength
            refas = []      # atlas line amplitude
            rej_win_w = []  # win rejected atlas line wavelength
            rej_win_a = []  # win rejected atlas line amplitude
            rej_fit_w = []  # fit rejected atlas line wavelength
            rej_fit_a = []  # fit rejected atlas line amplitude
            rej_par_w = []  # par rejected atlas line wavelength
            rej_par_a = []  # par rejected atlas line amplitude
            rej_dup_w = []  # dup rejected atlas line wavelength
            rej_dup_a = []  # dup rejected atlas line wavelength
            rej_fnt_w = []  # fnt rejected atlas line wavelength
            rej_fnt_a = []  # fnt rejected atlas line wavelength
            nrej = 0
            # look at each arc spectrum line
            for i, pk in enumerate(spec_cent):
                if pk <= minwav or pk >= maxwav:
                    continue
                # get atlas pixel position corresponding to arc line
                try:
                    line_x = [ii for ii, v in enumerate(atwave) if v >= pk][0]
                    # get window around atlas line to fit
                    minow, maxow, count = get_line_window(
                        atspec, line_x, logger=(self.logger if verbose
                                                else None),
                        frac_max=frac_max)
                except IndexError:
                    count = 0
                    minow = None
                    maxow = None
                    self.logger.warning("line at edge: %d, %.2f, %.f2f" %
                                        (i, pk, max(atwave)))
                # is resulting window large enough for fitting?
                if count < 5 or not minow or not maxow:
                    # keep track of fit rejected lines
                    rej_win_w.append(pk)
                    rej_win_a.append(spec_hgt[i])
                    nrej += 1
                    self.logger.info("Atlas window rejected for line %.3f" % pk)
                    continue
                # get data to fit
                yvec = atspec[minow:maxow + 1]
                xvec = atwave[minow:maxow + 1]
                # attempt Gaussian fit
                try:
                    fit, _ = curve_fit(gaus, xvec, yvec,
                                       p0=[spec_hgt[i], pk, 1.],
                                       maxfev=5000)
                except RuntimeError as e:
                    # keep track of Gaussian fit rejected lines
                    rej_fit_w.append(pk)
                    rej_fit_a.append(spec_hgt[i])
                    nrej += 1
                    self.logger.info("Atlas Gaussian fit rejected for line "
                                     "%.3f" % pk)
                    if verbose:
                        tb_str = traceback.format_exception(etype=type(e),
                                                            value=e,
                                                            tb=e.__traceback__)
                        self.logger.info("".join(tb_str))
                    continue
                # get interpolation function of atlas line
                int_line = interpolate.interp1d(xvec, yvec, kind='cubic',
                                                bounds_error=False,
                                                fill_value='extrapolate')
                # use very dense pixel sampling
                x_dense = np.linspace(min(xvec), max(xvec), num=1000)
                # resample line with dense sampling
                y_dense = int_line(x_dense)
                # get peak amplitude and wavelength
                pki = y_dense.argmax()
                pkw = x_dense[pki]
                # calculate some diagnostic parameters for the line
                # how many atlas pixels have we moved?
                xoff = abs(pkw - fit[1]) / self.action.args.refdisp
                # what is the wavelength offset in Angstroms?
                woff = abs(pkw - pk)
                # what fraction of the canonical fit width is the line?
                wrat = abs(fit[2]) / fwid  # can be neg or pos
                # current criteria for these diagnostic parameters
                if woff > 5. or xoff > 1.5 or wrat > 1.1:
                    # keep track of par rejected atlas lines
                    rej_par_w.append(pk)
                    rej_par_a.append(spec_hgt[i])
                    nrej += 1
                    self.logger.info("Atlas line parameters rejected for line "
                                     "%.3f" % pk)
                    self.logger.info("woff = %.3f, xoff = %.2f, wrat = %.3f" %
                                     (woff, xoff, wrat))
                    continue
                # check for duplicated lines
                if pkw in refws:
                    rej_dup_w.append(pk)
                    rej_dup_a.append(spec_hgt[i])
                    nrej += 1
                    self.logger.info("Atlas line duplicated for line %.3f" % pk)
                    continue
                self.logger.info("Atlas line accepted: %.3f" % pkw)
                refws.append(pkw)
                refas.append(y_dense[pki])
            # eliminate faintest lines if we have a large number
            self.logger.info("number of remaining lines: %d" % len(refas))
            if len(refas) > 400:
                # sort on flux
                sf = np.argsort(refas)
                refws = np.asarray(refws)[sf]
                refas = np.asarray(refas)[sf]
                # remove faintest two-thirds
                hlim = int(len(refas) * 0.67)
                refws = refws[hlim:]
                refas = refas[hlim:]
                rej_fnt_w = refws[:hlim]
                rej_fnt_a = refas[:hlim]
                # sort back onto wavelength
                sw = np.argsort(refws)
                refws = refws[sw].tolist()
                refas = refas[sw].tolist()

            self.logger.info("Using %d generated lines" % len(refws))

            # plot final list of Atlas lines and show rejections
            norm_fac = np.nanmax(atspec)
            if self.config.instrument.plot_level >= 1:
                p = figure(title=plotlabel(self.action.args) +
                           "ATLAS LINES Ngood = %d, Nrej = %d" % (len(refws),
                                                                  nrej),
                           x_axis_label="Wavelength (A)",
                           y_axis_label="Normalized Flux",
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.line(subwvals, subyvals / np.nanmax(subyvals),
                       legend_label='RefArc', color='lightgray')
                p.line(atwave, atspec / norm_fac, legend_label='Atlas',
                       color='blue')
                # Rejected: window bad
                if len(rej_win_w) > 0:
                    p.diamond(rej_win_w, rej_win_a / norm_fac,
                              legend_label='WinRej', color='cyan', size=8)
                # Rejected: fit failure
                if len(rej_fit_w) > 0:
                    p.diamond(rej_fit_w, rej_fit_a / norm_fac,
                              legend_label='FitRej', color='red', size=8)
                # Rejected: line parameter outside range
                if len(rej_par_w) > 0:
                    p.diamond(rej_par_w, rej_par_a / norm_fac,
                              legend_label='ParRej', color='orange', size=8)
                # Rejected: duplicated line
                if len(rej_dup_w) > 0:
                    p.diamond(rej_dup_w, rej_dup_a / norm_fac,
                              legend_label='DupRej', color='brown', size=8)
                # Rejected: faint
                if len(rej_fnt_w) > 0:
                    p.diamond(rej_fnt_w, rej_fnt_a / norm_fac,
                              legend_label='FntRej', color='magenta', size=8)
                # Kept
                p.diamond(refws, refas / norm_fac, legend_label='Kept',
                          color='green', size=16)
                p.line([minwav, minwav], [-0.1, 1.1], legend_label='WavLim',
                       color='brown')
                p.line([maxwav, maxwav], [-0.1, 1.1], color='brown')
                p.x_range = Range1d(min([min(subwvals), minwav-10.]),
                                    max(subwvals))
                p.y_range = Range1d(-0.04, 1.04)
                bokeh_plot(p, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)
                save_plot(p, filename=atfnam+".png")

        # output directory
        output_dir = os.path.join(self.config.instrument.cwd,
                                  self.config.instrument.output_directory)
        # write out final atlas line list
        atlines = np.array([refws, refas])
        atlines = atlines.T
        with open(os.path.join(output_dir, atfnam + '.txt'), 'w') as atlfn:
            np.savetxt(atlfn, atlines, fmt=['%12.3f', '%12.3f'])
        # store wavelengths, fluxes
        self.action.args.at_wave = refws
        self.action.args.at_flux = refas
        self.logger.info("Final atlas list has %d lines" % len(refws))

        log_string = GetAtlasLines.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class GetAtlasLines()
