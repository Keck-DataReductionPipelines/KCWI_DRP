from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_plotting import get_plot_lims, oplot_slices, \
    set_plot_lims, save_plot
from kcwidrp.primitives.GetAtlasLines import get_line_window, gaus
from kcwidrp.primitives.kcwi_file_primitives import plotlabel
from kcwidrp.core.kcwi_interactive import ask_next, pause, NEXT, QUIT

import numpy as np
from scipy.signal.windows import boxcar
import scipy as sp
from scipy.optimize import curve_fit
from scipy.interpolate import interpolate
from scipy.stats import sigmaclip
from kcwidrp.core.kcwi_plotting import figure
from kcwidrp.core import kcwi_plot_colors as pc
from kcwidrp.core.kcwi_plotting import oplot_cwave
from kcwidrp.core.kcwi_plotting import Range1d, LinearAxis


class SolveArcs(BasePrimitive):
    """
    Solve individual bar arc spectra for wavelength.

    For each bar, identifies the atlas lines found in GetAtlasLines.py in the
    observed spectrum and use the pixel positions along with the atlas
    wavelengths to determine the wavelength solution for the bar.

    Uses the following configuration parameters:

        * FRACMAX: fraction of line maximum to use for fitting peak.  Defaults to 0.5 for BLUE and 0.25 for RED.
        * LINETHRESH: the threshhold intensity for finding observed lines. Defaults to 100 for BLUE and 10 for RED (can also be set on the command line).

    Outputs diagnostic plots of the fitting and stores the coefficients for
    later use.

    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        self.action.args.fincoeff = []
        self.action.args.xsvals = None
        self.action.args.av_bar_sig = []
        self.action.args.st_bar_sig = []
        self.action.args.av_bar_nls = []
        self.action.args.st_bar_nls = []

    def _pre_condition(self):
        self.logger.info("Checking for master arc")
        if 'MARC' in self.action.args.ccddata.header['IMTYPE']:
            return True
        else:
            return False

    def _perform(self):
        """Solve individual arc bar spectra for wavelength"""
        self.logger.info("Solving individual arc spectra")
        # plot control booleans
        master_inter = (self.config.instrument.plot_level >= 2)
        do_inter = (self.config.instrument.plot_level >= 3)
        plab = plotlabel(self.action.args)
        # output control
        verbose = (self.config.instrument.verbose > 1)
        frac_max = self.config.instrument.FRACMAX
        self.logger.info("Finding line windows using fraction of line max of "
                         "%.2f" % frac_max)

        # Bar statistics
        bar_sig = []
        bar_nls = []
        # set thresh for finding lines
        hgt = self.config.instrument.LINETHRESH
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
        # Get poly order
        if self.action.args.dichroic_fraction <= 0.6:
            def_poly_order = 2
        elif 0.6 < self.action.args.dichroic_fraction < 0.75:
            def_poly_order = 3
        else:
            def_poly_order = 4
        poly_order = def_poly_order
        for ib, b in enumerate(self.context.arcs):
            self.logger.info("FITTING BAR %d" % ib)
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
                        logger=(self.logger if verbose else None),
                        frac_max=frac_max)
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
                    if line_x < minow or line_x > maxow:
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
                    # par_bounds = ([f0*0.9, np.min(xvec), 0.5],
                    #              [f0*1.1, np.max(xvec), 2.5])
                    # Gaussian fit
                    try:
                        fit, _ = curve_fit(gaus, xvec, yvec, p0=par_start,
                                           maxfev=5000)
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
                    if abs(cent - peak) > 0.8:
                        # keep track of rejected line
                        rej_wave.append(aw)
                        rej_flux.append(self.action.args.at_flux[iw])
                        nrej += 1
                        if verbose:
                            self.logger.info("Arc peak - cent offset = %.2f "
                                             "rejected for line %.3f" %
                                             (abs(cent - peak), aw))
                        continue
                    if plt_line[max_index] < hgt:
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
                            title=plab + "ATLAS/ARC LINE FITS" + ptitle,
                            x_axis_label="Wavelength (A)",
                            y_axis_label="Relative Flux",
                            plot_width=self.config.instrument.plot_width,
                            plot_height=self.config.instrument.plot_height)
                        ylim = [0, np.nanmax(yvec)]
                        p.line(atwave[atx0:atx1], atspec[atx0:atx1] * atnorm,
                               color=pc.REFERENCE, legend_label='Atlas')
                        p.circle(atwave[atx0:atx1], atspec[atx0:atx1] * atnorm,
                                 color=pc.REFERENCE, legend_label='Atlas')
                        p.line([aw, aw], ylim, color=pc.REFERENCE,
                               line_dash='dashed', legend_label='AtCntr')
                        p.x_range = Range1d(start=min(wvec), end=max(wvec))
                        p.extra_x_ranges = {"pix": Range1d(start=min(xvec),
                                                           end=max(xvec))}
                        p.add_layout(LinearAxis(x_range_name="pix",
                                                axis_label="CCD Y pix"),
                                     'above')
                        p.line(xplot, plt_line, color=pc.FIT,
                               legend_label='Arc', x_range_name="pix")
                        p.circle(xvec, yvec, legend_label='Arc', color=pc.DATA,
                                 x_range_name="pix")
                        ylim = [0, np.nanmax(plt_line)]
                        p.line([cent, cent], ylim, color=pc.FOUND,
                               legend_label='Cntr', line_dash='dashed',
                               x_range_name="pix")
                        p.line([sp_pk_x, sp_pk_x], ylim, color=pc.FOUND,
                               legend_label='Gpeak', line_dash='dotted',
                               x_range_name="pix")
                        p.line([peak, peak], ylim, color=pc.FOUND,
                               legend_label='Peak',
                               x_range_name="pix")
                        p.y_range.start = 0
                        bokeh_plot(p, self.context.bokeh_session)

                        q = ask_next(
                            self.logger,
                            what="Line fit: bar %d, line %d/%d, arc peak at %.1f px, atlas %.2f A." % (
                                ib, iw + 1, len(self.action.args.at_wave),
                                peak, aw),
                            help=("What you see:\n"
                                  "  Title \"ATLAS/ARC LINE FITS Bar# <n> - line <i>/<N> ...\". Two x axes:\n"
                                  "  wavelength (bottom, atlas) and CCD y pixel (top, arc), matched over\n"
                                  "  this window by the central wavelength fit.\n"
                                  "  Atlas:\n"
                                  "    Black line and circles: atlas spectrum.\n"
                                  "    Black dashed vertical: the atlas line.\n"
                                  "  Arc:\n"
                                  "    Blue circles: smoothed data in the fit window.\n"
                                  "    Red line: interpolation.\n"
                                  "    Magenta vertical: interpolated peak, the position used.\n"
                                  "    Magenta dashed: centroid.\n"
                                  "    Magenta dotted: Gaussian peak.\n"
                                  "What to look for:\n"
                                  "  Good: one clean arc line whose three magenta verticals coincide, at\n"
                                  "  the same place as the black atlas line.\n"
                                  "  Bad: blended or distorted lines, or an arc line far from the atlas\n"
                                  "  line (misidentification); the residual clipping may then reject it.\n"
                                  "Your options:\n"
                                  "  <cr>      show the next line of this bar.\n"
                                  "  q / quit  stop the line plots for all bars (bar plots continue).\n"
                                  "  Display only: no answer changes the reduction.\n"
                                  "Default:\n"
                                  "  <cr> shows the next line. Below plot_level 3 these are not shown."),
                            who="SolveArcs")
                        if q == QUIT:
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
            n_points = len(arc_pix_dat)
            self.logger.info("Fitting wavelength solution starting with %d "
                             "lines after rejecting %d lines" %
                             (n_points, nrej))
            # Fit wavelengths
            if n_points < 2:
                self.logger.warning("Not enough points for wavelength "
                                    "solution!  Using central coeffs")
                # store final fit coefficients
                self.action.args.fincoeff.append(coeff)
                # store statistics
                bar_sig.append(0.0)
                bar_nls.append(n_points)
                continue
            elif n_points < def_poly_order:
                poly_order = n_points - 1
            else:
                poly_order = def_poly_order
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
            # only reject if we have enough points
            if len(arc_pix_dat) > (poly_order + 1):
                while (max_resid > 2.5 * wsig or max_resid >= 10.) and it < 25:
                    arc_dat = []    # arc line pixel values
                    arc_fdat = []   # arc line flux data
                    at_dat = []     # atlas line wavelength values
                    at_fdat = []    # atlas line flux data
                    # trim largest outlier
                    for il, rsd in enumerate(resid):
                        if abs(rsd) < max_resid and abs(rsd) < 10.:
                            # append data for line that passed cut
                            arc_dat.append(arc_pix_dat[il])
                            arc_fdat.append(arc_int_dat[il])
                            at_dat.append(at_wave_dat[il])
                            at_fdat.append(at_flux_dat[il])
                        else:
                            if verbose:
                                self.logger.info("It%d REJ: %d, %.2f, %.3f, "
                                                 "%.3f" % (it, il,
                                                           arc_pix_dat[il],
                                                           at_wave_dat[il],
                                                           rsd))
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
            if poly_order < def_poly_order:
                nins = def_poly_order - poly_order
                wfit = np.insert(wfit, 0, np.zeros(nins, dtype=float))
            self.action.args.fincoeff.append(wfit)
            # store statistics
            bar_sig.append(wsig)
            bar_nls.append(len(arc_pix_dat))
            # do plotting?
            if master_inter and ib == next_bar_to_plot:
                # plot bar fit residuals
                # dispersion of this bar (A/px) for the RMS in px
                bar_disp = abs(wfit[-2])
                ptitle = " for Bar %03d, Slice %02d, RMS = %.3f A " \
                         "(%.2f px), N = %d" % \
                         (ib, int(ib / 5), wsig, wsig / bar_disp,
                          len(arc_pix_dat))
                p = figure(title=plab + "RESIDUALS" + ptitle,
                           x_axis_label="Wavelength (A)",
                           y_axis_label="Fit - Atlas (A)",
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.diamond(at_wave_dat, resid, color=pc.DATA,
                          legend_label='Kept', size=8)
                if rej_rsd_wave:
                    p.diamond(rej_rsd_wave, rej_rsd, color=pc.OUTLIER,
                              legend_label='RejRsd', size=8)
                if rej_wave:
                    p.scatter(rej_wave, np.zeros(len(rej_wave)), marker='x',
                              color=pc.MASKED, legend_label='RejFit', size=6)
                wlim = [self.action.args.atminwave, self.action.args.atmaxwave]
                xlim = get_plot_lims(wlim, padding=0.02)
                ylim = get_plot_lims(list(resid)+list(rej_rsd))
                p.line(wlim, [0., 0.], color=pc.THRESHOLD, line_dash='dotted')
                p.line(wlim, [wsig, wsig], color=pc.MEAN, line_dash='dashed',
                       legend_label='+-RMS')
                p.line(wlim, [-wsig, -wsig], color=pc.MEAN, line_dash='dashed')
                p.line([wlim[0], wlim[0]], ylim, color=pc.LIMIT,
                       legend_label='WavLim')
                p.line([wlim[1], wlim[1]], ylim, color=pc.LIMIT)
                oplot_cwave(p, self.action.args.cwave, ylim,
                            legend_label='CWAV')
                set_plot_lims(p, xlim=xlim, ylim=ylim)
                bokeh_plot(p, self.context.bokeh_session)
                pause(self.config.instrument, self.logger,
                      what="Wavelength-fit residuals of this bar, with kept and rejected lines.",
                      help=("What you see:\n"
                            "  Title \"RESIDUALS for Bar <n>, Slice <s>, RMS = <r> A (<p> px),\n"
                            "  N = <k>\": for each line, the wavelength the fit gives at the arc\n"
                            "  line's pixel minus the atlas wavelength (A). RMS is that of the kept\n"
                            "  lines; the px value divides it by this bar's dispersion.\n"
                            "  Blue diamonds: lines kept.\n"
                            "  Gold diamonds: lines rejected as outliers in the iterative clipping.\n"
                            "  Grey x (at 0): lines rejected before the fit (no clean window,\n"
                            "    failed Gaussian, centroid off).\n"
                            "  Grey dotted: zero.\n"
                            "  Red dashed: +- RMS.\n"
                            "  Purple verticals (\"WavLim\"): the wavelength limits of the atlas line\n"
                            "    list (as in P11); lines are only used between them.\n"
                            "  Orange dash-dot-dot: CWAVE.\n"
                            "  No WAVGOOD lines: WAVGOOD is only known after SolveGeom.\n"
                            "What to look for:\n"
                            "  Good: blue points scatter randomly about zero over the whole range;\n"
                            "  RMS (A) well below this bar's dispersion, so the px value in the\n"
                            "  title is well below 1.\n"
                            "  Bad: a trend or curve with wavelength (misidentified lines or too\n"
                            "  low a polynomial order), large gaps, many gold or grey lines (often\n"
                            "  a wrong atlas offset).\n"
                            "Your options:\n"
                            "  <cr>  show the fit overlay for this bar. Display only.\n"
                            "Default:\n"
                            "  <cr> continues. Below plot_level 2 this plot is not shown."),
                      who="SolveArcs")

                # overplot atlas and bar using fit wavelengths
                p = figure(title=plab + "ATLAS/ARC FIT" + ptitle,
                           x_axis_label="Wavelength (A)",
                           y_axis_label="Flux",
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                bwav = pwfit(self.action.args.xsvals)
                ylim = [np.nanmin(b), np.nanmax(b)]
                atnorm = np.nanmax(b) / np.nanmax(atspec)
                # atlas first, so the arc being inspected is drawn on top
                p.line(atwave, atspec * atnorm, color=pc.REFERENCE,
                       legend_label='Atlas')
                p.line(bwav, b, color=pc.DATA, legend_label='Arc')
                p.diamond(arc_wave_fit, arc_int_dat, color=pc.DATA, size=8)
                oplot_cwave(p, self.action.args.cwave, ylim,
                            legend_label='CWAV')
                p.diamond(at_wave, at_flux * atnorm, legend_label='Kept',
                          color=pc.FOUND, size=8)
                if rej_rsd_wave:
                    p.diamond(rej_rsd_wave, [rj*atnorm for rj in rej_rsd_flux],
                              color=pc.OUTLIER, legend_label='RejRsd', size=6)
                p.scatter(rej_wave, [rj*atnorm for rj in rej_flux], marker='x',
                          color=pc.MASKED, legend_label='RejFit', size=6)
                bokeh_plot(p, self.context.bokeh_session)
                q = ask_next(
                    self.logger,
                    what="Arc spectrum of this bar on its final wavelength solution, over the atlas.",
                    help=("What you see:\n"
                          "  Title \"ATLAS/ARC FIT for Bar <n> ...\".\n"
                          "  Blue line: this bar's arc on its final wavelength solution.\n"
                          "  Blue diamonds: lines used in the fit, at their fitted wavelengths.\n"
                          "  Black: atlas spectrum (scaled).\n"
                          "  Magenta diamonds (\"Kept\"): every line of the atlas line list.\n"
                          "  Gold: lines rejected as residual outliers.\n"
                          "  Grey x: lines rejected before the fit.\n"
                          "  Orange dash-dot-dot: CWAVE.\n"
                          "What to look for:\n"
                          "  Good: blue arc lines fall on the black atlas lines over the whole\n"
                          "  range, with blue diamonds on most magenta ones.\n"
                          "  Bad: blue lines shifted from black in part of the range: the\n"
                          "  solution for this bar is wrong there.\n"
                          "Your options:\n"
                          "  <cr>      show the next bar.\n"
                          "  <number>  jump to that bar (after the current one, up to 119); at\n"
                          "            plot_level 3 the line plots follow the jump.\n"
                          "  q / quit  stop the bar review (line plots stop too). The summary\n"
                          "            plots (coefficients, lines per bar, RMS) still follow.\n"
                          "  Display only: no answer changes the reduction.\n"
                          "Default:\n"
                          "  <cr> shows the next bar. Below plot_level 2 these are not shown."),
                    n_items=len(self.context.arcs), current=ib,
                    who="SolveArcs")
                if q == QUIT:
                    master_inter = False
                elif q == NEXT:
                    next_bar_to_plot = ib + 1
                else:
                    next_bar_to_plot = q

        # Plot final results
        nbars = self.config.instrument.NBARS
        # plot output name stub
        pfname = "arc_%05d_%s_%s_%s_tf%02d" % (
            self.action.args.ccddata.header['FRAMENO'],
            self.action.args.illum, self.action.args.grating,
            self.action.args.ifuname, int(100*self.config.instrument.TAPERFRAC))

        # Plot coefs
        if self.config.instrument.plot_level >= 1:
            ylabs = ['Ang/px^4', 'Ang/px^3', 'Ang/px^2', 'Ang/px',
                     'Ang']
            ylabs = ylabs[-(def_poly_order+1):]
            for ic in reversed(
                    range(len(self.action.args.fincoeff[0]))):
                # collect bar values for this particular coefficient
                coef = []
                for c in self.action.args.fincoeff:
                    try:
                        coef.append(c[ic])
                    except IndexError:
                        coef.append(0.0)
                # some stats
                cf_av = float(np.nanmean(coef))
                cf_st = float(np.nanstd(coef))
                cn = poly_order - ic
                if cn > 0:
                    ptitle = plab + "COEF %d OF C0-C%d VALUES <C%d> = " \
                                    "%.3g +- %.3g" % (cn, poly_order, cn,
                                                      cf_av, cf_st)
                else:
                    ptitle = plab + "COEF %d OF C0-C%d VALUES <C%d> = " \
                                    "%.2f +- %.2f" % (cn, poly_order, cn,
                                                      cf_av, cf_st)
                self.logger.info(ptitle)
                p = figure(title=ptitle, x_axis_label="Bar #",
                           y_axis_label="Coef %d (%s)" % (cn, ylabs[ic]),
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)

                p.diamond(list(range(nbars)), coef, color=pc.DATA, size=8)
                xlim = [-1, nbars]
                ylim = get_plot_lims(coef, clip=False)
                p.xgrid.grid_line_color = None
                oplot_slices(p, ylim)
                set_plot_lims(p, xlim=xlim, ylim=ylim)
                bokeh_plot(p, self.context.bokeh_session)
                pause(self.config.instrument, self.logger,
                      what="One coefficient of the wavelength solution for every bar.",
                      help=("What you see:\n"
                            "  Title \"COEF <k> OF C0-C<n> VALUES <Ck> = <mean> +- <std>\": coefficient\n"
                            "  k of each bar's wavelength polynomial in CCD pixel, against bar\n"
                            "  number. The polynomial has n + 1 coefficients, C0 to Cn: C0 is the\n"
                            "  wavelength (A) at the first pixel, C1 the dispersion (A/px), higher\n"
                            "  terms the curvature. One plot per coefficient.\n"
                            "  Blue diamonds: values.\n"
                            "  Black dash-dot lines: slice boundaries.\n"
                            "What to look for:\n"
                            "  Good: a smooth pattern within and between slices, no isolated\n"
                            "  points.\n"
                            "  Bad: a bar far off its neighbours has a poor solution (compare its\n"
                            "  N lines and RMS in the next plots).\n"
                            "Your options:\n"
                            "  <cr>  next coefficient (or the next plot). Display only.\n"
                            "Default:\n"
                            "  <cr> continues. At plot_level 1 the plot is shown for plot_pause\n"
                            "  seconds; at plot_level 0 it is not drawn."),
                      who="SolveArcs")
                # save coefficients plot
                save_plot(p, filename=pfname + '_coef%d.png' % cn)

        # Plot number of lines fit
        self.action.args.av_bar_nls = float(np.nanmean(bar_nls))
        self.action.args.st_bar_nls = float(np.nanstd(bar_nls))
        ptitle = plab + \
            "FIT STATS <Nlns> = %.1f +- %.1f" % (self.action.args.av_bar_nls,
                                                 self.action.args.st_bar_nls)
        p = figure(title=ptitle, x_axis_label="Bar #",
                   y_axis_label="N Lines",
                   plot_width=self.config.instrument.plot_width,
                   plot_height=self.config.instrument.plot_height)
        p.diamond(list(range(nbars)), bar_nls, color=pc.DATA, size=8,
                  legend_label='Bars')
        xlim = [-1, nbars]
        ylim = get_plot_lims(bar_nls)
        # mean of each slice's 5 bars
        for isl in range(nbars // 5):
            sl_av = float(np.nanmean(bar_nls[isl*5:isl*5+5]))
            p.line([isl*5 - 0.5, isl*5 + 4.5], [sl_av, sl_av],
                   color=pc.SUBSET_MEAN, line_width=3,
                   legend_label='Slice mean')
        self.logger.info("<N Lines> = %.1f +- %.1f" %
                         (self.action.args.av_bar_nls,
                          self.action.args.st_bar_nls))
        p.line(xlim, [self.action.args.av_bar_nls,
                      self.action.args.av_bar_nls], color=pc.MEAN,
               legend_label='Mean')
        p.line(xlim, [(self.action.args.av_bar_nls -
                       self.action.args.st_bar_nls),
                      (self.action.args.av_bar_nls -
                       self.action.args.st_bar_nls)], color=pc.MEAN,
               line_dash='dashed', legend_label='+-1 sigma')
        p.line(xlim, [(self.action.args.av_bar_nls +
                       self.action.args.st_bar_nls),
                      (self.action.args.av_bar_nls +
                       self.action.args.st_bar_nls)], color=pc.MEAN,
               line_dash='dashed')
        p.xgrid.grid_line_color = None
        oplot_slices(p, ylim)
        set_plot_lims(p, xlim=xlim, ylim=ylim)
        if self.config.instrument.plot_level >= 1:
            bokeh_plot(p, self.context.bokeh_session)
            pause(self.config.instrument, self.logger,
                  what="Number of arc lines used in the wavelength fit of every bar.",
                  help=("What you see:\n"
                        "  Title \"FIT STATS <Nlns> = <mean> +- <std>\": number of lines in each\n"
                        "  bar's final fit, against bar number.\n"
                        "  Blue diamonds (\"Bars\"): values.\n"
                        "  Red (\"Mean\"): mean of all bars.\n"
                        "  Red dashed: mean +- 1 standard deviation.\n"
                        "  Black thick, one per slice (\"Slice mean\"): mean of the slice's 5 bars.\n"
                        "  Black dash-dot: slice boundaries.\n"
                        "What to look for:\n"
                        "  Good: similar numbers for all bars and slices.\n"
                        "  Bad: bars or slices with far fewer lines than the rest; their\n"
                        "  solutions are less constrained, especially at the ends of the range.\n"
                        "  A red mean far below what this grating usually gives means the line\n"
                        "  identification failed: check the atlas offset (P39, P40) and the\n"
                        "  line list (P11).\n"
                        "Your options:\n"
                        "  <cr>  continue. Display only.\n"
                        "Default:\n"
                        "  <cr> continues. At plot_level 1 the plot is shown for plot_pause\n"
                        "  seconds; at plot_level 0 it is not drawn (the PNG is still saved)."),
                  who="SolveArcs")
        # save N lines plot
        save_plot(p, filename=pfname + '_nlines.png')

        # Plot fit sigmas
        self.action.args.av_bar_sig = float(np.nanmean(bar_sig))
        self.action.args.st_bar_sig = float(np.nanstd(bar_sig))
        self.logger.info("<STD>     = %.3f +- %.3f (A)" %
                         (self.action.args.av_bar_sig,
                          self.action.args.st_bar_sig))

        # mean dispersion (A/px) for the RMS in px
        av_disp = float(np.nanmean([abs(c[-2]) for c in
                                    self.action.args.fincoeff]))
        ptitle = plab + \
            "FIT STATS <RMS> = %.3f +- %.3f A (%.2f px)" % (
                self.action.args.av_bar_sig, self.action.args.st_bar_sig,
                self.action.args.av_bar_sig / av_disp)
        p = figure(title=ptitle, x_axis_label="Bar #", y_axis_label="RMS (A)",
                   plot_width=self.config.instrument.plot_width,
                   plot_height=self.config.instrument.plot_height)
        p.diamond(list(range(nbars)), bar_sig, color=pc.DATA, size=8,
                  legend_label='Bars')
        xlim = [-1, nbars]
        ylim = get_plot_lims(bar_sig)
        p.line(xlim, [self.action.args.av_bar_sig,
                      self.action.args.av_bar_sig], color=pc.MEAN,
               legend_label='Mean')
        p.line(xlim, [(self.action.args.av_bar_sig -
                       self.action.args.st_bar_sig),
                      (self.action.args.av_bar_sig -
                       self.action.args.st_bar_sig)], color=pc.MEAN,
               line_dash='dashed', legend_label='+-1 sigma')
        p.line(xlim, [(self.action.args.av_bar_sig +
                       self.action.args.st_bar_sig),
                      (self.action.args.av_bar_sig +
                       self.action.args.st_bar_sig)], color=pc.MEAN,
               line_dash='dashed')
        p.xgrid.grid_line_color = None
        oplot_slices(p, ylim)
        set_plot_lims(p, xlim=xlim, ylim=ylim)
        if self.config.instrument.plot_level >= 1:
            bokeh_plot(p, self.context.bokeh_session)
            pause(self.config.instrument, self.logger,
                  what="RMS of the wavelength fit of every bar (Angstrom).",
                  help=("What you see:\n"
                        "  Title \"FIT STATS <RMS> = <mean> +- <std> A (<p> px)\": RMS residual\n"
                        "  (A) of each bar's final fit, against bar number. A bar's RMS is that\n"
                        "  of its kept lines' residuals: the wavelength the fit gives at the arc\n"
                        "  line's pixel minus the atlas wavelength (as in P43). The px value\n"
                        "  divides the mean by the mean dispersion.\n"
                        "  Blue diamonds (\"Bars\"): values.\n"
                        "  Red (\"Mean\"): mean.\n"
                        "  Red dashed: mean +- 1 standard deviation.\n"
                        "  Black dash-dot: slice boundaries.\n"
                        "What to look for:\n"
                        "  Good: uniformly low RMS, well below the dispersion, so the px value\n"
                        "  in the title is well below 1.\n"
                        "  Bad: isolated bars or whole slices with high RMS: check them with\n"
                        "  the bar review (plot_level 2) before trusting those slices.\n"
                        "Your options:\n"
                        "  <cr>  continue. Display only.\n"
                        "Default:\n"
                        "  <cr> continues. At plot_level 1 the plot is shown for plot_pause\n"
                        "  seconds; at plot_level 0 it is not drawn (the PNG is still saved)."),
                  who="SolveArcs")

        # save residual plot
        save_plot(p, filename=pfname + '_resid.png')

        log_string = SolveArcs.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: def _perform(self):
