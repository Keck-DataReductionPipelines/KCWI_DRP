from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_plotting import get_plot_lims, oplot_slices, \
    set_plot_lims, save_plot
from kcwidrp.primitives.GetAtlasLines import get_line_window, gaus

import numpy as np
from scipy.signal.windows import boxcar
import scipy as sp
from scipy.optimize import curve_fit
from scipy.interpolate import interpolate
from scipy.stats import sigmaclip
from bokeh.plotting import figure
from bokeh.models import Range1d, LinearAxis
import time


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
        """Solve individual arc bar spectra for wavelength"""
        self.logger.info("Solving individual arc spectra")
        # plot control booleans
        master_inter = (self.config.instrument.plot_level >= 2)
        do_inter = (self.config.instrument.plot_level >= 3)
        # output control
        verbose = (self.config.instrument.verbose > 1)

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
        next_bar_to_plot = 0
        poly_order = 4
        for ib, b in enumerate(self.context.arcs):
            # Starting with pascal shifted coeffs from fit_center()
            coeff = self.action.args.twkcoeff[ib]
            # get bar wavelengths
            bw = np.polyval(coeff, self.action.args.xsvals)
            # smooth spectrum according to slicer
            if 'Small' in self.action.args.ifuname:
                # no smoothing for Small slicer
                bspec = b
            else:
                if 'Large' in self.action.args.ifuname:
                    # max smoothing for Large slicer
                    win = boxcar(5)
                else:
                    # intermediate smoothing for Medium slicer
                    win = boxcar(3)
                # do the smoothing
                bspec = sp.signal.convolve(b, win, mode='same') / sum(win)
            # store values to fit
            at_wave_dat = []    # atlas line wavelengths
            at_flux_dat = []    # atlas line peak fluxes
            arc_pix_dat = []    # arc line pixel positions
            arc_int_dat = []    # arc line pixel intensities
            rej_wave = []       # rejected line wavelengths
            rej_flux = []       # rejected line fluxes
            gaus_sig = []
            nrej = 0
            # loop over lines
            for iw, aw in enumerate(self.action.args.at_wave):
                # get window for this line
                try:
                    # get arc line initial pixel position
                    line_x = [i for i, v in enumerate(bw) if v >= aw][0]
                    # get window for arc line
                    minow, maxow, count = get_line_window(
                        bspec, line_x, thresh=hgt,
                        logger=(self.logger if verbose else None))
                    # do we have enough points to fit?
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
                    # get data to fit
                    yvec = bspec[minow:maxow + 1]
                    xvec = self.action.args.xsvals[minow:maxow + 1]
                    wvec = bw[minow:maxow + 1]
                    f0 = max(yvec)
                    par_start = [f0, np.nanmean(xvec), 1.0]
                    par_bounds = ([f0*0.9, np.min(xvec), 0.5],
                                  [f0*1.1, np.max(xvec), 2.5])
                    # Gaussian fit
                    try:
                        fit, _ = curve_fit(gaus, xvec, yvec, p0=par_start)
                        #  bounds=par_bounds, method='trf')
                        sp_pk_x = fit[1]
                        gaus_sig.append(fit[2])
                    except (RuntimeError, ValueError):
                        rej_wave.append(aw)
                        rej_flux.append(self.action.args.at_flux[iw])
                        nrej += 1
                        if verbose:
                            self.logger.info(
                                "Arc Gaussian fit rejected for line %.3f" % aw)
                        # sp_pk_x = line_x
                        continue

                    # get interpolation of arc line
                    int_line = interpolate.interp1d(xvec, yvec, kind='cubic',
                                                    bounds_error=False,
                                                    fill_value='extrapolate')
                    # use very dense sampling
                    xplot = np.linspace(min(xvec), max(xvec), num=1000)
                    # re-sample line with dense sampling
                    plt_line = int_line(xplot)
                    # get peak position
                    max_index = plt_line.argmax()
                    peak = xplot[max_index]
                    # calculate centroid
                    cent = np.sum(xvec * yvec) / np.sum(yvec)
                    # how different is the centroid from the peak?
                    if abs(cent - peak) > 0.7:
                        # keep track of rejected line
                        rej_wave.append(aw)
                        rej_flux.append(self.action.args.at_flux[iw])
                        nrej += 1
                        if verbose:
                            self.logger.info("Arc peak - cent offset = %.2f "
                                             "rejected for line %.3f" %
                                             (abs(cent - peak), aw))
                        continue
                    if plt_line[max_index] < 100:
                        # keep track of rejected line
                        rej_wave.append(aw)
                        rej_flux.append(self.action.args.at_flux[iw])
                        nrej += 1
                        if verbose:
                            self.logger.info("Arc peak too low = %.2f "
                                             "rejected for line %.3f" %
                                             (plt_line[max_index], aw))
                        continue
                    # store surviving line data
                    arc_pix_dat.append(peak)
                    arc_int_dat.append(plt_line[max_index])
                    at_wave_dat.append(aw)
                    at_flux_dat.append(self.action.args.at_flux[iw])
                    # plot, if requested
                    if do_inter and ib == next_bar_to_plot:
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
            self.logger.info("")
            self.logger.info("Fitting wavelength solution starting with %d "
                             "lines after rejecting %d lines" %
                             (len(arc_pix_dat), nrej))
            # Fit wavelengths
            # Get poly order
            if self.action.args.dichroic_fraction <= 0.6:
                poly_order = 2
            elif 0.6 < self.action.args.dichroic_fraction < 0.75:
                poly_order = 3
            else:
                poly_order = 4
            self.logger.info("Fitting with polynomial order %d" % poly_order)
            # Initial fit
            wfit = np.polyfit(arc_pix_dat, at_wave_dat, poly_order)
            pwfit = np.poly1d(wfit)
            arc_wave_fit = pwfit(arc_pix_dat)
            # fit residuals
            resid = arc_wave_fit - at_wave_dat
            resid_c, low, upp = sigmaclip(resid, low=3., high=3.)
            wsig = resid_c.std()
            # maximum outlier
            max_resid = np.max(abs(resid))
            self.logger.info("wsig: %.3f, max_resid: %.3f" % (wsig, max_resid))
            # keep track of rejected lines
            rej_rsd = []        # rejected line residuals
            rej_rsd_wave = []   # rejected line wavelengths
            rej_rsd_flux = []   # rejected line fluxes
            # iteratively remove outliers
            it = 0
            while max_resid > 2.5 * wsig and it < 25:
                arc_dat = []    # arc line pixel values
                arc_fdat = []   # arc line flux data
                at_dat = []     # atlas line wavelength values
                at_fdat = []    # atlas line flux data
                # trim largest outlier
                for il, rsd in enumerate(resid):
                    if abs(rsd) < max_resid:
                        # append data for line that passed cut
                        arc_dat.append(arc_pix_dat[il])
                        arc_fdat.append(arc_int_dat[il])
                        at_dat.append(at_wave_dat[il])
                        at_fdat.append(at_flux_dat[il])
                    else:
                        if verbose:
                            self.logger.info("It%d REJ: %d, %.2f, %.3f, %.3f" %
                                             (it, il, arc_pix_dat[il],
                                              at_wave_dat[il], rsd))
                        # keep track of rejected lines
                        rej_rsd_wave.append(at_wave_dat[il])
                        rej_rsd_flux.append(at_flux_dat[il])
                        rej_rsd.append(rsd)
                # copy cleaned data back into input arrays
                arc_pix_dat = arc_dat.copy()
                arc_int_dat = arc_fdat.copy()
                at_wave_dat = at_dat.copy()
                at_flux_dat = at_fdat.copy()
                # refit cleaned data
                wfit = np.polyfit(arc_pix_dat, at_wave_dat, poly_order)
                # new wavelength function
                pwfit = np.poly1d(wfit)
                # new wavelengths for arc lines
                arc_wave_fit = pwfit(arc_pix_dat)
                # calculate residuals of arc lines
                resid = arc_wave_fit - at_wave_dat
                # get statistics
                resid_c, low, upp = sigmaclip(resid, low=3., high=3.)
                wsig = resid_c.std()
                # maximum outlier
                max_resid = np.max(abs(resid))
                # wsig = np.nanstd(resid)
                it += 1
            # END while max_resid > 3.5 * wsig and it < 5:
            # log arc bar results
            self.logger.info("")
            self.logger.info("BAR %03d, Slice = %02d, RMS = %.3f, N = %d" %
                             (ib, int(ib / 5), wsig, len(arc_pix_dat)))
            self.logger.info(
                "Nits: %d, wsig: %.3f, max_resid: %.3f" % (it, wsig, max_resid))
            self.logger.info("NRejRsd: %d, NRejFit: %d" % (len(rej_rsd_wave),
                                                           len(rej_wave)))
            self.logger.info("Line width median sigma: %.2f px" %
                             np.nanmedian(gaus_sig))
            self.logger.info("Coefs: " + ' '.join(['%.6g' % (c,)
                                                   for c in reversed(wfit)]))
            # store final fit coefficients
            self.action.args.fincoeff.append(wfit)
            # store statistics
            bar_sig.append(wsig)
            bar_nls.append(len(arc_pix_dat))
            # do plotting?
            if master_inter and ib == next_bar_to_plot:
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
                p.diamond(arc_wave_fit, arc_int_dat, color='darkgrey', size=8)
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
                q = input("Next? <int> or <cr>, q - quit: ")
                if 'Q' in q.upper():
                    master_inter = False
                else:
                    try:
                        next_bar_to_plot = int(q)
                    except ValueError:
                        next_bar_to_plot = ib + 1

        # Plot final results

        # plot output name stub
        pfname = "arc_%05d_%s_%s_%s_tf%02d" % (
            self.action.args.ccddata.header['FRAMENO'],
            self.action.args.illum, self.action.args.grating,
            self.action.args.ifuname, int(100*self.config.instrument.TAPERFRAC))

        # Plot coefs
        if self.config.instrument.plot_level >= 1:
            ylabs = ['Ang/px^4', 'Ang/px^3', 'Ang/px^2', 'Ang/px',
                     'Ang']
            ylabs = ylabs[-(poly_order+1):]
            for ic in reversed(
                    range(len(self.action.args.fincoeff[0]))):
                cn = poly_order - ic
                ptitle = self.action.args.plotlabel + "COEF %d VALUES" % cn
                p = figure(title=ptitle, x_axis_label="Bar #",
                           y_axis_label="Coef %d (%s)" % (cn, ylabs[ic]),
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                coef = []
                for c in self.action.args.fincoeff:
                    coef.append(c[ic])
                p.diamond(list(range(120)), coef, size=8)
                xlim = [-1, 120]
                ylim = get_plot_lims(coef)
                p.xgrid.grid_line_color = None
                oplot_slices(p, ylim)
                set_plot_lims(p, xlim=xlim, ylim=ylim)
                bokeh_plot(p, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)
                # save coefficients plot
                save_plot(p, filename=pfname + '_coef%d.png' % cn)

        # Plot number of lines fit
        self.action.args.av_bar_nls = float(np.nanmean(bar_nls))
        self.action.args.st_bar_nls = float(np.nanstd(bar_nls))
        ptitle = self.action.args.plotlabel + \
            "FIT STATS <Nlns> = %.1f +- %.1f" % (self.action.args.av_bar_nls,
                                                 self.action.args.st_bar_nls)
        p = figure(title=ptitle, x_axis_label="Bar #",
                   y_axis_label="N Lines",
                   plot_width=self.config.instrument.plot_width,
                   plot_height=self.config.instrument.plot_height)
        p.diamond(list(range(120)), bar_nls, size=8)
        xlim = [-1, 120]
        ylim = get_plot_lims(bar_nls)
        self.logger.info("<N Lines> = %.1f +- %.1f" %
                         (self.action.args.av_bar_nls,
                          self.action.args.st_bar_nls))
        p.line(xlim, [self.action.args.av_bar_nls,
                      self.action.args.av_bar_nls], color='red')
        p.line(xlim, [(self.action.args.av_bar_nls -
                       self.action.args.st_bar_nls),
                      (self.action.args.av_bar_nls -
                       self.action.args.st_bar_nls)], color='green',
               line_dash='dashed')
        p.line(xlim, [(self.action.args.av_bar_nls +
                       self.action.args.st_bar_nls),
                      (self.action.args.av_bar_nls +
                       self.action.args.st_bar_nls)], color='green',
               line_dash='dashed')
        p.xgrid.grid_line_color = None
        oplot_slices(p, ylim)
        set_plot_lims(p, xlim=xlim, ylim=ylim)
        if self.config.instrument.plot_level >= 1:
            bokeh_plot(p, self.context.bokeh_session)
            if self.config.instrument.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.config.instrument.plot_pause)
        # save N lines plot
        save_plot(p, filename=pfname + '_nlines.png')

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
        ylim = get_plot_lims(bar_sig)
        p.line(xlim, [self.action.args.av_bar_sig,
                      self.action.args.av_bar_sig], color='red')
        p.line(xlim, [(self.action.args.av_bar_sig -
                       self.action.args.st_bar_sig),
                      (self.action.args.av_bar_sig -
                       self.action.args.st_bar_sig)], color='green',
               line_dash='dashed')
        p.line(xlim, [(self.action.args.av_bar_sig +
                       self.action.args.st_bar_sig),
                      (self.action.args.av_bar_sig +
                       self.action.args.st_bar_sig)], color='green',
               line_dash='dashed')
        p.xgrid.grid_line_color = None
        oplot_slices(p, ylim)
        set_plot_lims(p, xlim=xlim, ylim=ylim)
        if self.config.instrument.plot_level >= 1:
            bokeh_plot(p, self.context.bokeh_session)
            if self.config.instrument.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.config.instrument.plot_pause)

        # save residual plot
        save_plot(p, filename=pfname + '_resid.png')

        log_string = SolveArcs.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: def _perform(self):
