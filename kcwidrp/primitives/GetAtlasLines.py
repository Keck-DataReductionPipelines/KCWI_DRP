from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments


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
        self.action.args.at_flux = refas
        # plot results
        norm_fac = np.nanmax(atspec)
        p = figure(title=self.action.args.plotlabel +
                   "ATLAS LINES Ngood = %d, Nrej = %d" % (len(refws), nrej),
                   x_axis_label="Wavelength (A)",
                   y_axis_label="Normalized Flux",
                   plot_width=self.config.instrument.plot_width,
                   plot_height=self.config.instrument.plot_height)
        p.line(subwvals, subyvals / np.nanmax(subyvals), legend='RefArc',
               color='lightgray')
        p.line(atwave, atspec / norm_fac, legend='Atlas', color='blue')
        # Rejected: nearby neighbor
        p.diamond(rej_neigh_w, rej_neigh_y / norm_fac, legend='NeighRej',
                  color='cyan', size=8)
        # Rejected: fit failure
        p.diamond(rej_fit_w, rej_fit_y / norm_fac, legend='FitRej',
                  color='red', size=8)
        # Rejected: line parameter outside range
        p.diamond(rej_par_w, rej_par_a / norm_fac, legend='ParRej',
                  color='orange', size=8)
        p.diamond(refws, refas / norm_fac, legend='Kept', color='green',
                  size=10)
        p.x_range = Range1d(min(subwvals), max(subwvals))
        if self.config.instrument.plot_level >= 1:
            bokeh_plot(p)
            if self.config.instrument.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                pl.pause(self.config.instrument.plot_pause)
            export_png(p, "atlas_lines_%s_%s_%s_%05d.png" %
                       (self.action.args.illum, self.action.args.grating,
                        self.action.args.ifuname,
                        self.action.args.ccddata.header['FRAMENO']))
        self.logger.info("Final atlas list has %d lines" % len(refws))

        logstr = GetAtlasLines.__module__ + "." + GetAtlasLines.__qualname__
        self.action.args.ccddata.header['HISTORY'] = logstr
        self.logger.info(logstr)

        return self.action.args
    # END: class GetAtlasLines()


