from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_plotting import save_plot

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
    if x1 + 1 >= nx:
        if verbose:
            print("max check: high edge hit")
        return None, None, 0
    while y[x1+1] > mx:
        x1 += 1
        count += 1
        if x1 + 1 >= nx:
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
                                               p0=[y[i], x[i], 1.])
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


class GetAtlasLines(BasePrimitive):
    """Get relevant atlas line positions and wavelengths"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        self.action.args.atminrow = None
        self.action.args.atmaxrow = None
        self.action.args.atminwave = None
        self.action.args.atmaxwave = None
        self.action.args.at_wave = None
        self.action.args.at_flux = None

    def _perform(self):
        self.logger.info("Finding isolated atlas lines")
        # get atlas wavelength range
        #
        # get pixel values (no longer centered in the middle)
        specsz = len(self.context.arcs[self.config.instrument.REFBAR])
        xvals = np.arange(0, specsz)
        # min, max rows
        minrow = 50
        maxrow = specsz - 50
        # wavelength range
        mnwvs = []
        mxwvs = []
        # Get wavelengths
        for b in range(self.config.instrument.NBARS):
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
        # get amplitude threshold
        c, low, upp = sigmaclip(atspec, low=3.5, high=3.5)
        atflxsig = c.std()
        atflxmn = c.mean()
        # plot atlas flux histogram
        if self.config.instrument.plot_level >= 2:
            hist, edges = np.histogram(atspec, range=(low, upp),
                                       density=False, bins=50)
            p = figure(title='Atlas flux',
                       x_axis_label='flx', y_axis_label='N',
                       plot_width=self.config.instrument.plot_width,
                       plot_height=self.config.instrument.plot_height)
            p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
                   fill_color="navy", line_color="white", alpha=0.5)
            p.line([atflxmn-atflxsig, atflxmn-atflxsig], [0, np.max(hist)],
                   color='red',
                   legend_label="Sigma")
            p.line([atflxmn+atflxsig, atflxmn+atflxsig], [0, np.max(hist)],
                   color='red')
            p.line([atflxmn, atflxmn], [0, np.max(hist)], color='green',
                   legend_label='Mean')
            p.y_range.start = 0
            bokeh_plot(p, self.context.bokeh_session)
            input("Next? <cr>: ")
        # make sure atlas lines are significant w.r.t. background
        ampl_thresh = atflxmn + 2.5 * atflxsig
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
                fit, _ = curve_fit(gaus, xvec, yvec, p0=[spec_hgt[i], pk, 1.])
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
        self.action.args.at_flux = refas
        # plot results
        norm_fac = np.nanmax(atspec)
        if self.config.instrument.plot_level >= 1:
            p = figure(title=self.action.args.plotlabel +
                       "ATLAS LINES Ngood = %d, Nrej = %d" % (len(refws), nrej),
                       x_axis_label="Wavelength (A)",
                       y_axis_label="Normalized Flux",
                       plot_width=self.config.instrument.plot_width,
                       plot_height=self.config.instrument.plot_height)
            p.line(subwvals, subyvals / np.nanmax(subyvals),
                   legend_label='RefArc', color='lightgray')
            p.line(atwave, atspec / norm_fac, legend_label='Atlas',
                   color='blue')
            # Rejected: nearby neighbor
            p.diamond(rej_neigh_w, rej_neigh_y / norm_fac,
                      legend_label='NeighRej', color='cyan', size=8)
            # Rejected: fit failure
            p.diamond(rej_fit_w, rej_fit_y / norm_fac, legend_label='FitRej',
                      color='red', size=8)
            # Rejected: line parameter outside range
            p.diamond(rej_par_w, rej_par_a / norm_fac, legend_label='ParRej',
                      color='orange', size=8)
            p.diamond(refws, refas / norm_fac, legend_label='Kept',
                      color='green', size=10)
            p.x_range = Range1d(min(subwvals), max(subwvals))
            bokeh_plot(p, self.context.bokeh_session)
            if self.config.instrument.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.config.instrument.plot_pause)
            save_plot(p, filename="arc_%05d_%s_%s_%s_atlines.png" %
                                  (self.action.args.ccddata.header['FRAMENO'],
                                   self.action.args.illum,
                                   self.action.args.grating,
                                   self.action.args.ifuname))
        self.logger.info("Final atlas list has %d lines" % len(refws))

        log_string = GetAtlasLines.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class GetAtlasLines()
