from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot

import numpy as np
from scipy.signal.windows import boxcar
import scipy as sp
from scipy.optimize import curve_fit
from scipy.interpolate import interpolate
from scipy.stats import sigmaclip
from bokeh.plotting import figure
from bokeh.models import Range1d, LinearAxis
from bokeh.io import export_png
import time


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


class SolveArcs(BasePrimitive):
    """Solve the bar arc wavelengths"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        self.action.args.fincoeff = []
        self.action.args.xsvals = None
        self.action.args.av_bar_sig = []
        self.action.args.st_bar_sig = []
        self.action.args.av_bar_nls = []
        self.action.args.st_bar_nls = []

    def _perform(self):
        self.logger.info("Solving individual arc spectra")
        if self.config.instrument.plot_level >= 2:
            master_inter = True
        else:
            master_inter = False
        if self.config.instrument.plot_level >= 3:
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
        # convert list into ndarray
        at_wave = np.asarray(self.action.args.at_wave)
        at_flux = np.asarray(self.action.args.at_flux)
        # get x values starting at zero pixels
        self.action.args.xsvals = np.arange(0, len(
            self.context.arcs[self.config.instrument.REFBAR]))
        # loop over arcs and generate a wavelength solution for each
        for ib, b in enumerate(self.context.arcs):
            # print("")
            # self.logger.info("BAR %d" % ib)
            # print("")
            # Starting pascal shifted coeffs from fit_center()
            coeff = self.action.args.twkcoeff[ib]
            # get bar wavelengths
            bw = np.polyval(coeff, self.action.args.xsvals)
            # smooth spectrum according to slicer
            if 'Small' in self.action.args.ifuname:
                bspec = b
            else:
                if 'Large' in self.action.args.ifuname:
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
            at_flux_dat = []
            arc_pix_dat = []
            arc_int_dat = []
            rej_wave = []
            rej_flux = []
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
                        rej_flux.append(self.action.args.at_flux[iw])
                        nrej += 1
                        if verbose:
                            self.logger.info("Arc window rejected for line %.3f"
                                             % aw)
                        continue
                    # check if window no longer contains initial value
                    if minow > line_x > maxow:
                        rej_wave.append(aw)
                        rej_flux.append(self.action.args.at_flux[iw])
                        nrej += 1
                        if verbose:
                            self.logger.info(
                                "Arc window wandered off for line %.3f" % aw)
                        continue
                    yvec = bspec[minow:maxow + 1]
                    xvec = self.action.args.xsvals[minow:maxow + 1]
                    wvec = bw[minow:maxow + 1]
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
                        rej_flux.append(self.action.args.at_flux[iw])
                        nrej += 1
                        if verbose:
                            self.logger.info("Arc peak - cent offset rejected "
                                             "for line %.3f" % aw)
                        continue
                    # store data
                    arc_pix_dat.append(peak)
                    arc_int_dat.append(fit[0])  # plt_line[max_index])
                    at_wave_dat.append(aw)
                    at_flux_dat.append(self.action.args.at_flux[iw])
                    # plot, if requested
                    if do_inter:
                        ptitle = " Bar# %d - line %3d/%3d: xc = %.1f, " \
                                 "Wave = %9.2f" % \
                                 (ib, (iw + 1), len(self.action.args.at_wave),
                                  peak, aw)
                        atx0 = [i for i, v in enumerate(atwave)
                                if v >= min(wvec)][0]
                        atx1 = [i for i, v in enumerate(atwave)
                                if v >= max(wvec)][0]
                        atnorm = np.nanmax(yvec) / np.nanmax(atspec[atx0:atx1])
                        p = figure(
                            title=self.action.args.plotlabel +
                            "ATLAS/ARC LINE FITS" + ptitle,
                            x_axis_label="Wavelength (A)",
                            y_axis_label="Relative Flux",
                            plot_width=self.config.instrument.plot_width,
                            plot_height=self.config.instrument.plot_height)
                        ylim = [0, np.nanmax(yvec)]
                        p.line(atwave[atx0:atx1], atspec[atx0:atx1] * atnorm,
                               color='blue', legend_label='Atlas')
                        p.circle(atwave[atx0:atx1], atspec[atx0:atx1] * atnorm,
                                 color='green', legend_label='Atlas')
                        p.line([aw, aw], ylim, color='red',
                               legend_label='AtCntr')
                        p.x_range = Range1d(start=min(wvec), end=max(wvec))
                        p.extra_x_ranges = {"pix": Range1d(start=min(xvec),
                                                           end=max(xvec))}
                        p.add_layout(LinearAxis(x_range_name="pix",
                                                axis_label="CCD Y pix"),
                                     'above')
                        p.line(xplot, plt_line, color='black',
                               legend_label='Arc', x_range_name="pix")
                        p.circle(xvec, yvec, legend_label='Arc', color='red',
                                 x_range_name="pix")
                        ylim = [0, np.nanmax(plt_line)]
                        p.line([cent, cent], ylim, color='green',
                               legend_label='Cntr', line_dash='dashed',
                               x_range_name="pix")
                        p.line([sp_pk_x, sp_pk_x], ylim, color='magenta',
                               legend_label='Gpeak', line_dash='dashdot',
                               x_range_name="pix")
                        p.line([peak, peak], ylim, color='black',
                               legend_label='Peak', line_dash='dashdot',
                               x_range_name="pix")
                        p.y_range.start = 0
                        bokeh_plot(p, self.context.bokeh_session)

                        q = input(ptitle + " - Next? <cr>, q to quit: ")
                        if 'Q' in q.upper():
                            do_inter = False
                except IndexError:
                    if verbose:
                        self.logger.info(
                            "Atlas line not in observation: %.2f" % aw)
                    rej_wave.append(aw)
                    rej_flux.append(self.action.args.at_flux[iw])
                    nrej += 1
                    continue
                except ValueError:
                    if verbose:
                        self.logger.info(
                            "Interpolation error for line at %.2f" % aw)
                    rej_wave.append(aw)
                    rej_flux.append(self.action.args.at_flux[iw])
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
            rej_rsd_flux = []
            # Iteratively remove outliers
            for it in range(4):
                # self.logger.info("Iteration %d" % it)
                arc_dat = []
                arc_fdat = []
                at_dat = []
                at_fdat = []
                # Trim outliers
                for il, rsd in enumerate(resid):
                    if low < rsd < upp:
                        arc_dat.append(arc_pix_dat[il])
                        arc_fdat.append(arc_int_dat[il])
                        at_dat.append(at_wave_dat[il])
                        at_fdat.append(at_flux_dat[il])
                    else:
                        # self.logger.info("REJ: %d, %.2f, %.3f" %
                        #              (il, arc_pix_dat[il], at_wave_dat[il]))
                        rej_rsd_wave.append(at_wave_dat[il])
                        rej_rsd_flux.append(at_flux_dat[il])
                        rej_rsd.append(rsd)
                # refit
                arc_pix_dat = arc_dat.copy()
                arc_int_dat = arc_fdat.copy()
                at_wave_dat = at_dat.copy()
                at_flux_dat = at_fdat.copy()
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
            self.logger.info("RejRsd: %d, RejFit: %d" % (len(rej_rsd_wave),
                                                         len(rej_wave)))
            # print("")
            self.action.args.fincoeff.append(wfit)
            bar_sig.append(wsig)
            bar_nls.append(len(arc_pix_dat))
            # do plotting?
            if master_inter:
                # plot bar fit residuals
                ptitle = " for Bar %03d, Slice %02d, RMS = %.3f, N = %d" % \
                         (ib, int(ib / 5), wsig, len(arc_pix_dat))
                p = figure(title=self.action.args.plotlabel +
                           "RESIDUALS" + ptitle,
                           x_axis_label="Wavelength (A)",
                           y_axis_label="Fit - Inp (A)",
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.diamond(at_wave_dat, resid, legend_label='Rsd', size=8)
                if rej_rsd_wave:
                    p.diamond(rej_rsd_wave, rej_rsd, color='orange',
                              legend_label='Rej', size=8)
                xlim = [self.action.args.atminwave, self.action.args.atmaxwave]
                ylim = [np.nanmin(list(resid)+list(rej_rsd)),
                        np.nanmax(list(resid)+list(rej_rsd))]
                p.line(xlim, [0., 0.], color='black', line_dash='dotted')
                p.line(xlim, [wsig, wsig], color='gray', line_dash='dashdot')
                p.line(xlim, [-wsig, -wsig], color='gray', line_dash='dashdot')
                p.line([self.action.args.cwave, self.action.args.cwave],
                       ylim, legend_label='CWAV', color='magenta',
                       line_dash='dashdot')
                bokeh_plot(p, self.context.bokeh_session)
                input("Next? <cr>: ")

                # overplot atlas and bar using fit wavelengths
                p = figure(title=self.action.args.plotlabel +
                           "ATLAS/ARC FIT" + ptitle,
                           x_axis_label="Wavelength (A)",
                           y_axis_label="Flux",
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                bwav = pwfit(self.action.args.xsvals)
                p.line(bwav, b, color='darkgrey', legend_label='Arc')
                ylim = [np.nanmin(b), np.nanmax(b)]
                atnorm = np.nanmax(b) / np.nanmax(atspec)
                p.line(atwave, atspec * atnorm, color='blue',
                       legend_label='Atlas')
                p.line([self.action.args.cwave, self.action.args.cwave],
                       ylim, color='magenta', line_dash='dashdot',
                       legend_label='CWAV')
                p.diamond(at_wave, at_flux * atnorm, legend_label='Kept',
                          color='green', size=8)
                if rej_rsd_wave:
                    p.diamond(rej_rsd_wave, [rj*atnorm for rj in rej_rsd_flux],
                              color='orange', legend_label='RejRsd', size=6)
                p.diamond(rej_wave, [rj*atnorm for rj in rej_flux],
                          color='red', legend_label='RejFit', size=6)
                bokeh_plot(p, self.context.bokeh_session)
                q = input("Next? <cr>, q - quit: ")
                if 'Q' in q.upper():
                    master_inter = False
        # Plot final results
        # Plot fit sigmas
        self.action.args.av_bar_sig = float(np.nanmean(bar_sig))
        self.action.args.st_bar_sig = float(np.nanstd(bar_sig))
        self.logger.info("<STD>     = %.3f +- %.3f (A)" %
                         (self.action.args.av_bar_sig,
                          self.action.args.st_bar_sig))

        ptitle = self.action.args.plotlabel + \
            "FIT STATS <RMS> = %.3f +- %.3f" % (self.action.args.av_bar_sig,
                                                self.action.args.st_bar_sig)
        p = figure(title=ptitle, x_axis_label="Bar #", y_axis_label="RMS (A)",
                   plot_width=self.config.instrument.plot_width,
                   plot_height=self.config.instrument.plot_height)
        p.diamond(list(range(120)), bar_sig, size=8)
        xlim = [-1, 120]
        ylim = [np.nanmin(bar_sig), np.nanmax(bar_sig)]

        p.line(xlim, [self.action.args.av_bar_sig,
                      self.action.args.av_bar_sig], color='black')
        p.line(xlim, [(self.action.args.av_bar_sig -
                       self.action.args.st_bar_sig),
                      (self.action.args.av_bar_sig -
                       self.action.args.st_bar_sig)], color='black',
               line_dash='dotted')
        p.line(xlim, [(self.action.args.av_bar_sig +
                       self.action.args.st_bar_sig),
                      (self.action.args.av_bar_sig +
                       self.action.args.st_bar_sig)], color='black',
               line_dash='dotted')
        for ix in range(1, 24):
            sx = ix * 5 - 0.5
            p.line([sx, sx], ylim, color='black', line_dash='dashdot')
        if self.config.instrument.plot_level >= 1:
            bokeh_plot(p, self.context.bokeh_session)
            if self.config.instrument.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.config.instrument.plot_pause)
        export_png(p, filename="arc_%05d_resid_%s_%s_%s.png" %
                   (self.action.args.ccddata.header['FRAMENO'],
                    self.action.args.illum,
                    self.action.args.grating, self.action.args.ifuname))
        # Plot number of lines fit
        self.action.args.av_bar_nls = float(np.nanmean(bar_nls))
        self.action.args.st_bar_nls = float(np.nanstd(bar_nls))
        ptitle = self.action.args.plotlabel + \
            "FIT STATS <Nlns> = %.1f +- %.1f" % (self.action.args.av_bar_nls,
                                                 self.action.args.st_bar_nls)
        p = figure(title=ptitle, x_axis_label="Bar #", y_axis_label="N Lines",
                   plot_width=self.config.instrument.plot_width,
                   plot_height=self.config.instrument.plot_height)
        p.diamond(list(range(120)), bar_nls,  size=8)
        xlim = [-1, 120]
        ylim = [np.nanmin(bar_nls), np.nanmax(bar_nls)]
        self.logger.info("<N Lines> = %.1f +- %.1f" %
                         (self.action.args.av_bar_nls,
                          self.action.args.st_bar_nls))
        p.line(xlim, [self.action.args.av_bar_nls,
                      self.action.args.av_bar_nls], color='black')
        p.line(xlim, [(self.action.args.av_bar_nls -
                       self.action.args.st_bar_nls),
                      (self.action.args.av_bar_nls -
                       self.action.args.st_bar_nls)], color='black',
               line_dash='dotted')
        p.line(xlim, [(self.action.args.av_bar_nls +
                       self.action.args.st_bar_nls),
                      (self.action.args.av_bar_nls +
                       self.action.args.st_bar_nls)], color='black',
               line_dash='dotted')
        for ix in range(1, 24):
            sx = ix * 5 - 0.5
            p.line([sx, sx], ylim, color='black', line_dash='dashdot')

        if self.config.instrument.plot_level >= 1:
            bokeh_plot(p, self.context.bokeh_session)
            if self.config.instrument.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.config.instrument.plot_pause)
        export_png(p, filename="arc_%05d_nlines_%s_%s_%s.png" %
                   (self.action.args.ccddata.header['FRAMENO'],
                    self.action.args.illum,
                    self.action.args.grating, self.action.args.ifuname))
        # Plot coefs
        if self.config.instrument.plot_level >= 1:
            ylabs = ['Ang/px^4', 'Ang/px^3', 'Ang/px^2', 'Ang/px', 'Ang']
            for ic in reversed(range(len(self.action.args.fincoeff[0]))):
                ptitle = self.action.args.plotlabel + "COEF %d VALUES" % ic
                p = figure(title=ptitle, x_axis_label="Bar #",
                           y_axis_label="Coef %d (%s)" % (ic, ylabs[ic]),
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                coef = []
                for c in self.action.args.fincoeff:
                    coef.append(c[ic])
                p.diamond(list(range(120)), coef, size=8)
                ylim = [np.nanmin(coef), np.nanmax(coef)]
                for ix in range(1, 24):
                    sx = ix * 5 - 0.5
                    p.line([sx, sx], ylim, color='black', line_dash='dashdot')
                bokeh_plot(p, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)
                export_png(p, filename="arc_%05d_coef%d_%s_%s_%s.png" %
                           (self.action.args.ccddata.header['FRAMENO'], ic,
                            self.action.args.illum,
                            self.action.args.grating,
                            self.action.args.ifuname))

        log_string = SolveArcs.__module__ + "." + SolveArcs.__qualname__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class SolveArcs()
