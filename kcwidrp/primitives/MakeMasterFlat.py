from keckdrpframework.primitives.base_img import BaseImg
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer, strip_fname, plotlabel
from kcwidrp.core.kcwi_plotting import get_plot_lims
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_plotting import save_plot
from kcwidrp.core.bspline import Bspline
from kcwidrp.core.kcwi_plotting import figure
from kcwidrp.core.kcwi_plotting import Range1d
from kcwidrp.core.kcwi_interactive import ask_value, pause

import os
import numpy as np
from scipy.signal.windows import boxcar
import scipy as sp
from scipy.signal import find_peaks


def bm_ledge_position(cwave, dich):
    if dich:
        fit = [0.240742, 4060.56]
    else:
        fit = [0.240742, 4044.56]
    return fit[1] + fit[0] * cwave


class MakeMasterFlat(BaseImg):
    """
    Generate illumination correction from a stacked flat image.

    Uses b-spline fits along with geometry maps to generate a master image for
    illumination correction.  If the flat is internal, accounts for vignetting
    along one edge.  Also accounts for ledge seen in BM grating.

    Depending on the type of input stack, the following files are written out
    and entries are made in the proc file:

        * SFLAT - a \*_mflat.fits file and an MFLAT entry
        * SDOME - a \*_mdome.fits file and an MDOME entry
        * STWIF - a \*_mtwif.fits file and an MTWIF entry

"""

    def __init__(self, action, context):
        BaseImg.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if we can create a master flat based on the processing table
        """
        # check for flat stack to use in generating master flat
        self.logger.info("Checking precondition for MakeMasterFlat")
        tab = self.context.proctab.search_proctab(
            frame=self.action.args.ccddata,
            target_type=self.action.args.new_type, nearest=True)
        if len(tab) > 0:
            self.logger.info("already have %d master %s flats, "
                             " expecting 0" % (len(tab),
                                               self.action.args.new_type))
            return False
        else:
            self.logger.info("No %s master flat found, "
                             "checking for flat stack." %
                             self.action.args.new_type)
            self.stack_list = self.context.proctab.search_proctab(
                frame=self.action.args.ccddata,
                target_type=self.action.args.stack_type,
                target_group=self.action.args.groupid)
            self.logger.info(f"pre condition got {len(self.stack_list)},"
                             f" expecting {self.action.args.min_files}")
            # do we meet the criterion?
            if len(self.stack_list) >= 1:
                return True
            else:
                return False

    def _perform(self):
        """
        Returns an Argument() with the parameters that depends on this operation
        """
        self.logger.info("Creating master illumination correction")
        olab = plotlabel(self.action.args)

        suffix = self.action.args.new_type.lower()
        insuff = self.action.args.stack_type.lower()

        stack_list = list(self.stack_list['filename'])

        if len(stack_list) <= 0:
            self.logger.warning("No flat stack found!")
            return self.action.args

        # get root for maps
        tab = self.context.proctab.search_proctab(
            frame=self.action.args.ccddata, target_type='MARC',
            target_group=self.action.args.groupid)
        if len(tab) <= 0:
            self.logger.error("Geometry not solved!")
            return self.action.args

        mroot = strip_fname(tab['filename'][-1]) # ??????????????????

        # Wavelength map image
        wmf = mroot + '_wavemap.fits'
        self.logger.info("Reading image: %s" % wmf)
        wavemap = kcwi_fits_reader(
            os.path.join(self.config.instrument.cwd, 'redux', wmf))[0]

        # Slice map image
        slf = mroot + '_slicemap.fits'
        self.logger.info("Reading image: %s" % slf)
        slicemap = kcwi_fits_reader(os.path.join(
            self.config.instrument.cwd, 'redux', slf))[0]

        # Position map image
        pof = mroot + '_posmap.fits'
        self.logger.info("Reading image: %s" % pof)
        posmap = kcwi_fits_reader(os.path.join(
            self.config.instrument.cwd, 'redux', pof))[0]

        # Read in stacked flat image
        stname = strip_fname(stack_list[0]) + '_' + insuff + '.fits'

        self.logger.info("Reading image: %s" % stname)
        stacked = kcwi_fits_reader(os.path.join(
            self.config.instrument.cwd, 'redux', stname))[0]
        plab = " ".join(olab.split()[0:3]) + (
            " %d " % stacked.header['FRAMENO']
        ) + " ".join(olab.split()[4:])

        # get type of flat
        internal = ('SFLAT' in stacked.header['IMTYPE'])
        twiflat = ('STWIF' in stacked.header['IMTYPE'])
        domeflat = ('SDOME' in stacked.header['IMTYPE'])

        if internal:
            self.logger.info("Internal Flat")
        elif twiflat:
            self.logger.info("Twilight Flat")
        elif domeflat:
            self.logger.info("Dome Flat")
        else:
            self.logger.error("Flat of Unknown Type!")
            return self.action.args

        # knots per pixel
        knotspp = self.config.instrument.KNOTSPP

        # get image size
        ny = stacked.header['NAXIS2']

        # get binning
        xbin = self.action.args.xbinsize

        # Parameters for fitting

        if self.action.args.camera == 0:  # Blue
            # vignetted slice position range
            fitl = int(4/xbin)
            fitr = int(24/xbin)

            # flat fitting slice position range
            ffleft = int(10 / xbin)
            ffright = int(70 / xbin)

            corlim = 0

            # ymap not needed for Blue
            ymap = None
        else:   # Red
            # vignetted slice position range
            fitl = int(114 / xbin)
            fitr = int(134 / xbin)

            # flat fitting slice position range
            ffleft = int(70 / xbin)
            ffright = int(130 / xbin)

            corlim = int(140 / xbin)

            # Get ymap for trimming junk at ends
            ymap = posmap.copy()
            for i in range(ny):
                ymap.data[i, :] = float(i)

        # un-vignetted slice position range
        flatl = int(34 / xbin)
        flatr = int(72 / xbin)

        nrefx = int(ffright - ffleft)

        buffer = 6.0/float(xbin)

        # reference slice
        refslice = 9
        allidx = np.arange(int(140/xbin))
        newflat = stacked.data.copy()

        # dichroic fraction
        # try:
        #     dichroic_fraction = wavemap.header['DICHFRAC']
        # except KeyError:
        #     dichroic_fraction = 1.

        # get reference slice data
        q = [i for i, v in enumerate(slicemap.data.flat) if v == refslice]
        # get wavelength limits
        waves = wavemap.data.compress((wavemap.data > 0.).flat)
        waves = [waves.min(), waves.max()]
        self.logger.info("Wavelength limits: %.1f - %1.f" % (waves[0],
                                                             waves[1]))

        # correct vignetting if we are using internal flats
        if internal:
            self.logger.info("Internal flats require vignetting correction")
            # get good region for fitting
            if self.action.args.camera == 0:    # Blue
                wmin = waves[0]
                wmax = min([waves[1], 5620.])
            elif self.action.args.camera == 1:  # Red
                wmin = max([waves[0], 5580.])
                wmax = waves[1]
            else:
                self.logger.warning("Camera keyword not defined")
                wmin = waves[0]
                wmax = waves[1]
            # Use central fraction to calculate vignetting
            dw = (wmax - wmin) / 30.0   # fraction of full wavelength range
            wavemin = (wmin+wmax) / 2.0 - dw
            wavemax = (wmin+wmax) / 2.0 + dw
            self.logger.info("Using %.1f - %.1f A of slice %d" % (wavemin,
                                                                  wavemax,
                                                                  refslice))
            xflat = []
            yflat = []
            wflat = []
            qq = []
            for i in q:
                if wavemin < wavemap.data.flat[i] < wavemax:
                    xflat.append(posmap.data.flat[i])
                    yflat.append(stacked.data.flat[i])
                    wflat.append(wavemap.data.flat[i])
                    qq.append(i)
            # get un-vignetted portion
            qflat = [i for i, v in enumerate(xflat) if flatl <= v <= flatr]
            xflat = [xflat[i] for i in qflat]
            yflat = [yflat[i] for i in qflat]
            wflat = [wflat[i] for i in qflat]
            # sort on wavelength
            sw = np.argsort(wflat)
            ywflat = [yflat[i] for i in sw]
            wwflat = [wflat[i] for i in sw]
            ww0 = np.min(wwflat)
            # fit wavelength slope
            wavelinfit = np.polyfit(wwflat-ww0, ywflat, 2)
            wslfit = np.polyval(wavelinfit, wflat-ww0)
            # plot slope fit
            if self.config.instrument.plot_level >= 1:
                p = figure(title=plab + ' WAVE SLOPE FIT',
                           x_axis_label='wave px',
                           y_axis_label='counts',
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.circle(wwflat, ywflat, legend_label="Data")
                p.line(wflat, wslfit, line_color='red', line_width=3,
                       legend_label="Fit")
                bokeh_plot(p, self.context.bokeh_session)
                pause(self.config.instrument, self.logger,
                      what="Internal flat: lamp spectrum slope near the middle wavelength, reference slice.",
                      help=("What you see:\n"
                            '  Title "WAVE SLOPE FIT". Blue circles: counts of the un-vignetted part\n'
                            "  of the reference slice (slice 9) in a narrow window (1/15 of the\n"
                            "  range) around the middle wavelength; x is wavelength in A (the axis\n"
                            '  label says "wave px"). Red: quadratic fit, divided out before the\n'
                            "  vignetting is measured.\n"
                            "What to look for:\n"
                            "  Good: the red line follows the trend of the points.\n"
                            "  Bad: a fit pulled away by outliers, which then tilts the vignetting\n"
                            "  correction in the next plots.\n"
                            "Your options:\n"
                            "  <cr>  continue. Display only.\n"
                            "Default:\n"
                            "  <cr> continues. At plot_level 1 the plot is shown for plot_pause\n"
                            "  seconds; at plot_level 0 it is not drawn."),
                      who="MakeMasterFlat")
            # take out slope
            yflat = yflat / wslfit
            # now sort on slice position
            ss = np.argsort(xflat)
            xflat = [xflat[i] for i in ss]
            yflat = [yflat[i] for i in ss]
            # fit un-vignetted slope
            resflat = np.polyfit(xflat, yflat, 1)

            # select the points we will fit for the vignetting
            # get reference region
            xfit = [posmap.data.flat[i] for i in qq]
            yfit = [stacked.data.flat[i] for i in qq]
            wflat = [wavemap.data.flat[i] for i in qq]
            # take out wavelength slope
            yfit = yfit / np.polyval(wavelinfit, wflat-ww0)

            # select the vignetted region
            qfit = [i for i, v in enumerate(xfit) if fitl <= v <= fitr]
            xfit = [xfit[i] for i in qfit]
            yfit = [yfit[i] for i in qfit]
            # sort on slice position
            s = np.argsort(xfit)
            xfit = [xfit[i] for i in s]
            yfit = [yfit[i] for i in s]
            # fit vignetted slope
            resfit = np.polyfit(xfit, yfit, 1)
            # corrected data
            ycdata = stacked.data.flat[qq] / \
                np.polyval(wavelinfit, wavemap.data.flat[qq]-ww0)
            ycmin = 0.5     # np.min(ycdata)
            ycmax = 1.25    # np.max(ycdata)
            # compute the intersection
            xinter = -(resflat[1] - resfit[1]) / (resflat[0] - resfit[0])
            # plot slice profile and fits
            if self.config.instrument.plot_level >= 1:
                p = figure(title=plab + ' Vignetting',
                           x_axis_label='Slice Pos (px)',
                           y_axis_label='Ratio',
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.circle(posmap.data.flat[qq], ycdata, legend_label='Data')
                p.line(allidx, resfit[1] + resfit[0]*allidx,
                       line_color='purple', legend_label='Vign.')
                p.line(allidx, resflat[1] + resflat[0]*allidx, line_color='red',
                       legend_label='UnVign.')
                p.line([fitl, fitl], [ycmin, ycmax], line_color='blue')
                p.line([fitr, fitr], [ycmin, ycmax], line_color='blue')
                p.line([flatl, flatl], [ycmin, ycmax], line_color='green')
                p.line([flatr, flatr], [ycmin, ycmax], line_color='green')
                p.line([xinter-buffer, xinter-buffer], [ycmin, ycmax],
                       line_color='black')
                p.line([xinter + buffer, xinter + buffer], [ycmin, ycmax],
                       line_color='black')
                p.line([xinter, xinter], [ycmin, ycmax], line_color='red')
                p.y_range = Range1d(ycmin, ycmax)
                bokeh_plot(p, self.context.bokeh_session)
                pause(self.config.instrument, self.logger,
                      what="Internal flat: vignetting along the slice edge and its linear correction.",
                      help=("What you see:\n"
                            '  Title "Vignetting": slope-corrected flat ratio against position along\n'
                            "  the reference slice (px); y range fixed at 0.5 - 1.25.\n"
                            "  Blue circles: data. Purple: linear fit to the vignetted part (between\n"
                            "  the blue verticals). Red line: linear fit to the un-vignetted part\n"
                            "  (between the green verticals). Red vertical: where the two fits cross.\n"
                            "  Black verticals: buffer zone either side of the crossing.\n"
                            "What to look for:\n"
                            "  Good: each fit follows the data in its own region and they cross at\n"
                            "  the knee where the fall-off begins.\n"
                            "  Bad: fits that miss the data, or a crossing far from the knee: the\n"
                            "  internal flat will over- or under-correct that slice edge.\n"
                            "Your options:\n"
                            "  <cr>  continue. Display only.\n"
                            "Default:\n"
                            "  <cr> continues. At plot_level 1 the plot is shown for plot_pause\n"
                            "  seconds; at plot_level 0 it is not drawn."),
                      who="MakeMasterFlat")

            # figure out where the correction applies
            if self.action.args.camera == 0:
                qcor = [i for i, v in enumerate(posmap.data.flat)
                        if corlim <= v <= (xinter-buffer)]
            else:
                qcor = [i for i, v in enumerate(posmap.data.flat)
                        if (xinter + buffer) <= v <= corlim]
            # apply the correction!
            self.logger.info("Applying vignetting correction...")
            for i in qcor:
                newflat.flat[i] = (resflat[1]+resflat[0]*posmap.data.flat[i]) \
                                / (resfit[1]+resfit[0]*posmap.data.flat[i]) * \
                                stacked.data.flat[i]
            # now deal with the intermediate (buffer) region
            self.logger.info("Done, now handling buffer region")
            # get buffer points to fit in reference region
            qbff = [i for i in qq if (xinter-buffer) <=
                    posmap.data.flat[i] <= (xinter+buffer)]
            # get slice pos and data for buffer fitting
            xbuff = [posmap.data.flat[i] for i in qbff]
            ybuff = [stacked.data.flat[i] / np.polyval(wavelinfit,
                                                       wavemap.data.flat[i]-ww0)
                     for i in qbff]
            # sort on slice position
            ssp = np.argsort(xbuff)
            xbuff = [xbuff[i] for i in ssp]
            ybuff = [ybuff[i] for i in ssp]
            # fit buffer with low-order poly
            buffit = np.polyfit(xbuff, ybuff, 3)
            # plot buffer fit
            if self.config.instrument.plot_level >= 1:
                p = figure(title=plab + ' Buffer Region',
                           x_axis_label='Slice Pos (px)',
                           y_axis_label='Ratio',
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.circle(xbuff, ybuff)
                p.line(xbuff, np.polyval(buffit, xbuff), line_color='red')
                bokeh_plot(p, self.context.bokeh_session)
                pause(self.config.instrument, self.logger,
                      what="Internal flat: smooth fit across the buffer zone between corrected and uncorrected parts.",
                      help=("What you see:\n"
                            '  Title "Buffer Region": ratio data within the buffer zone around the\n'
                            "  vignetting knee, against slice position. Blue circles: data. Red:\n"
                            "  cubic fit used to join the corrected and uncorrected parts smoothly.\n"
                            "What to look for:\n"
                            "  Good: a smooth red curve through the points.\n"
                            "  Bad: a curve that swings away from the data at either end, which\n"
                            "  leaves a step in the flat at the edge of the buffer.\n"
                            "Your options:\n"
                            "  <cr>  continue. Display only.\n"
                            "Default:\n"
                            "  <cr> continues. At plot_level 1 the plot is shown for plot_pause\n"
                            "  seconds; at plot_level 0 it is not drawn."),
                      who="MakeMasterFlat")
            # get all buffer points in image
            qbuf = [i for i, v in enumerate(posmap.data.flat)
                    if (xinter-buffer) <= v <= (xinter+buffer)]
            # apply buffer correction to all buffer points in newflat
            for i in qbuf:
                newflat.flat[i] = \
                    (resflat[1] + resflat[0] * posmap.data.flat[i]) / \
                    np.polyval(buffit, posmap.data.flat[i]) * newflat.flat[i]
            self.logger.info("Vignetting correction complete.")

        self.logger.info("Fitting master illumination")
        # now fit master flat
        # get reference slice points
        if ymap is not None:
            self.logger.info("Using ymap to limit yrange of data in qref")
            qref = [i for i in q if ffleft <= posmap.data.flat[i] <= ffright and
                    50 <= ymap.data.flat[i] <= (ny-50)]
        else:
            qref = [i for i in q if ffleft <= posmap.data.flat[i] <= ffright]
        xfr = wavemap.data.flat[qref]
        yfr = newflat.flat[qref]
        # sort on wavelength
        s = np.argsort(xfr)
        xfr = xfr[s]
        yfr = yfr[s]

        wavegood0 = wavemap.header['WAVGOOD0']
        wavegood1 = wavemap.header['WAVGOOD1']

        # correction for BM where we see a ledge (except for twilights)
        if 'BM' in self.action.args.grating and not twiflat:
            ledge_wave = bm_ledge_position(self.action.args.cwave,
                                           self.action.args.dich)

            self.logger.info("BM ledge calculated wavelength "
                             "for ref slice = %.2f (A)" % ledge_wave)
            if wavegood0 <= ledge_wave <= wavegood1:
                self.logger.info("BM grating requires correction")
                qledge = [i for i, v in enumerate(xfr)
                          if ledge_wave-25 <= v <= ledge_wave+25]
                xledge = [xfr[i] for i in qledge]
                yledge = [yfr[i] for i in qledge]
                s = np.argsort(xledge)
                xledge = [xledge[i] for i in s]
                yledge = [yledge[i] for i in s]
                win = boxcar(250)
                smyledge = sp.signal.convolve(yledge,
                                              win, mode='same') / sum(win)
                ylmax = np.max(yledge)
                ylmin = np.min(yledge)
                fpoints = np.arange(0, 100) / 100. * 50 + (ledge_wave-25)
                ledgefit, ledgemsk = Bspline.iterfit(np.asarray(xledge),
                                                     smyledge, fullbkpt=fpoints,
                                                     upper=1, lower=1)
                ylfit, _ = ledgefit.value(np.asarray(fpoints))
                deriv = -(np.roll(ylfit, 1) - np.roll(ylfit, -1)) / 2.0
                # trim edges
                trm = int(len(deriv)/5)
                deriv = deriv[trm:-trm]
                xvals = fpoints[trm:-trm]
                peaks, properties = find_peaks(deriv, height=20)
                self.logger.info("%d Peak(s) found" % len(peaks))

                p = figure(title=plab +
                           ' Ledge', x_axis_label='Wavelength (A)',
                           y_axis_label='Value',
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.circle(xledge, smyledge, fill_color='green')
                p.line(fpoints, ylfit)
                ylim = get_plot_lims(ylfit)
                p.line([ledge_wave, ledge_wave], ylim, color='red')
                bokeh_plot(p, self.context.bokeh_session)
                pause(self.config.instrument, self.logger,
                      what="BM grating flat: the 'ledge' step near its predicted wavelength, reference slice.",
                      help=("What you see:\n"
                            '  Title "Ledge": reference-slice flat within +- 25 A of the wavelength\n'
                            "  where the BM grating's ledge (a step in the flat-field spectrum) is\n"
                            "  predicted from CWAVE. Green circles (blue outline): smoothed data.\n"
                            "  Blue line: spline fit. Red vertical: predicted ledge wavelength.\n"
                            "What to look for:\n"
                            "  Good: a clear step near the red line. Its position is found from the\n"
                            "  derivative in the next plot.\n"
                            "  Bad: no step, or several similar steps; then the peak choice in the\n"
                            "  next step needs care.\n"
                            "Your options:\n"
                            "  <cr>  continue. Display only.\n"
                            "Default:\n"
                            "  <cr> continues. Below plot_level 2 the plot is shown for plot_pause\n"
                            "  seconds (also at plot_level 0)."),
                      who="MakeMasterFlat")
                p = figure(title=plab + ' Deriv',
                           x_axis_label='px (Wavelength)',
                           y_axis_label='Value',
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                xx = list(range(len(deriv)))
                ylim = get_plot_lims(deriv)
                p.circle(xx, deriv)
                p.line([len(xx)/2, len(xx)/2], ylim, color='green')
                for pk in peaks:
                    p.line([pk, pk], ylim, color='red')
                bokeh_plot(p, self.context.bokeh_session)
                if len(peaks) == 0:
                    # No ledge detected: nothing to correct
                    self.logger.warning("No BM ledge peak found, "
                                        "skipping BM ledge correction")
                    ipk = None
                elif len(peaks) != 1:
                    self.logger.warning("Single peak not found!")
                    if self.config.instrument.plot_level >= 2:
                        auto_ipk = peaks[properties['peak_heights'].argmax()]

                        def parse_peak(text):
                            # today's test: str.isnumeric(), then int()
                            if not text.isnumeric():
                                raise ValueError("a non-negative integer "
                                                 "is needed")
                            return int(text)

                        def check_peak(value):
                            if value > len(xvals) - 1:
                                return "index must be from 0 to %d" % (
                                    len(xvals) - 1)
                            return None
                        ipk = ask_value(
                            self.logger,
                            "Peak? <int> (<cr> = highest peak, index %d)" %
                            auto_ipk,
                            help=("What you see:\n"
                                  '  Title "Deriv": derivative of the ledge fit, in 0.5 A steps, against\n'
                                  '  index 0-59 (axis "px (Wavelength)"; index i is at predicted ledge\n'
                                  "  wavelength - 15 A + 0.5 i A). Blue circles: derivative.\n"
                                  "  Green vertical (index 30): the predicted ledge wavelength.\n"
                                  "  Red verticals: peaks found (height >= 20); there is more than one.\n"
                                  "What to look for:\n"
                                  "  The ledge is a single sharp peak, usually the strongest one close to\n"
                                  "  the green line. Smaller or broader peaks come from noise or spectral\n"
                                  "  structure in the flat.\n"
                                  "Your options:\n"
                                  "  <index>  put the ledge at that index (read it off the x axis). The\n"
                                  "           step size is measured from linear fits 3-8 A either side,\n"
                                  "           and the reference-slice flat redward of it is divided by it.\n"
                                  "           This changes the master flat (*_mflat or *_mdome) and every\n"
                                  "           frame flat-fielded with it.\n"
                                  "  <cr>     use the automatic choice, the highest peak.\n"
                                  "Default:\n"
                                  "  <cr> takes the highest peak, as the reduction does below plot_level\n"
                                  "  2. The highest peak is not always the ledge: check the next plots\n"
                                  "  (BM Ledge Region) and rerun at plot_level 2 if the step remains."),
                            parse=parse_peak, default=auto_ipk,
                            effect="BM ledge peak set by user to index "
                                   "{value} (automatic choice: index %d)" %
                                   auto_ipk,
                            validate=check_peak, who="MakeMasterFlat",
                            what="Several BM ledge peaks found: enter the index (x value) of the peak that is the ledge.")
                    else:
                        # Pick the highest peak, using the peaks properties
                        self.logger.warning("Picking the highest peak. This may not be the correct peak!")
                        ipk = peaks[properties['peak_heights'].argmax()]
                else:
                    ipk = peaks[0]
                if ipk is not None:
                    apk = xvals[ipk]
                    if self.config.instrument.plot_level >= 3:
                        p = figure(
                            title=plab + ' Peak of ledge',
                            x_axis_label='Wave (A)',
                            y_axis_label='Value',
                            plot_width=self.config.instrument.plot_width,
                            plot_height=self.config.instrument.plot_height)
                        p.circle(xvals, deriv, legend_label='Data')
                        p.line([apk, apk], [-50, 200], line_color='red',
                               legend_label='Pk')
                        bokeh_plot(p, self.context.bokeh_session)
                        pause(self.config.instrument, self.logger,
                              what="BM ledge: derivative against wavelength with the chosen ledge position.",
                              help=("What you see:\n"
                                    '  Title "Peak of ledge": the derivative of the ledge fit against\n'
                                    '  wavelength (A). Blue circles: data. Red vertical ("Pk"): the ledge\n'
                                    "  position used.\n"
                                    "What to look for:\n"
                                    "  Good: the red line sits on the sharp peak you intended.\n"
                                    "  Bad: it sits on a noise peak; the correction will put a step into\n"
                                    "  the flat instead of removing one.\n"
                                    "Your options:\n"
                                    "  <cr>  continue. Display only.\n"
                                    "Default:\n"
                                    "  <cr> continues. Below plot_level 3 this plot is not shown."),
                              who="MakeMasterFlat")
                    xlow = apk - 3 - 5
                    xhi = apk - 3
                    zlow = apk + 3
                    zhi = apk + 3 + 5
                    qlow = [i for i, v in enumerate(fpoints) if xlow <= v <= xhi]
                    xlf = np.asarray([fpoints[i] for i in qlow])
                    ylf = np.asarray([ylfit[i] for i in qlow])
                    lowfit = np.polyfit(xlf, ylf, 1)
                    qhi = [i for i, v in enumerate(fpoints) if zlow <= v <= zhi]
                    xlf = np.asarray([fpoints[i] for i in qhi])
                    ylf = np.asarray([ylfit[i] for i in qhi])
                    hifit = np.polyfit(xlf, ylf, 1)
                    ratio = (hifit[1] + hifit[0] * apk) / \
                            (lowfit[1] + lowfit[0] * apk)
                    self.logger.info("BM ledge ratio: %.3f" % ratio)
                    # correct flat data
                    qcorr = [i for i, v in enumerate(xfr) if v >= apk]
                    for i in qcorr:
                        yfr[i] /= ratio
                    # plot BM ledge
                    if self.config.instrument.plot_level >= 1:
                        p = figure(
                            title=plab + ' BM Ledge Region',
                            x_axis_label='Wave (A)',
                            y_axis_label='Value',
                            plot_width=self.config.instrument.plot_width,
                            plot_height=self.config.instrument.plot_height)
                        # Input data
                        p.circle(xledge, yledge, fill_color='blue',
                                 legend_label='Data')
                        # correct input data
                        qcorrect = [i for i, v in enumerate(xledge) if v >= apk]
                        xplt = []
                        yplt = []
                        for i in qcorrect:
                            xplt.append(xledge[i])
                            yplt.append(yledge[i] / ratio)
                        p.circle(xplt, yplt, fill_color='orange',
                                 legend_label='Corrected')
                        p.line(fpoints, ylfit, line_color='red', legend_label='Fit')
                        p.line([xlow, xlow], [ylmin, ylmax], line_color='blue')
                        p.line([xhi, xhi], [ylmin, ylmax], line_color='blue')
                        p.line([zlow, zlow], [ylmin, ylmax], line_color='black')
                        p.line([zhi, zhi], [ylmin, ylmax], line_color='black')
                        p.line(fpoints, lowfit[1] + lowfit[0] * fpoints,
                               line_color='purple')
                        p.line(fpoints, hifit[1] + hifit[0] * fpoints,
                               line_color='green')
                        p.line([apk, apk], [ylmin, ylmax], line_color='green',
                               legend_label='Pk')
                        p.y_range = Range1d(ylmin, ylmax)
                        p.legend.location = 'top_left'
                        bokeh_plot(p, self.context.bokeh_session)
                        pause(self.config.instrument, self.logger,
                              what="BM ledge: flat data before and after the step correction.",
                              help=("What you see:\n"
                                    '  Title "BM Ledge Region", within +- 25 A of the predicted ledge.\n'
                                    "  Blue-filled circles: reference-slice flat data. Orange-filled\n"
                                    "  circles: the same data redward of the ledge after dividing by the\n"
                                    '  step ratio ("Corrected"). Red: spline fit to the smoothed data.\n'
                                    "  Blue verticals: window of the blueward linear fit (purple line);\n"
                                    "  black verticals: window of the redward linear fit (green line).\n"
                                    '  Green vertical ("Pk"): ledge position.\n'
                                    "What to look for:\n"
                                    "  Good: the orange points continue the level of the blue points below\n"
                                    "  the ledge, with no step left.\n"
                                    "  Bad: a step remains or is overcorrected: the wrong peak was taken.\n"
                                    "Your options:\n"
                                    "  <cr>  continue. Display only.\n"
                                    "Default:\n"
                                    "  <cr> continues. At plot_level 1 the plot is shown for plot_pause\n"
                                    "  seconds; at plot_level 0 it is not drawn."),
                              who="MakeMasterFlat")
        # END: handling BM grating ledge

        # if we are fitting a twilight flat, treat it like a sky image with a
        # larger number of knots
        if twiflat:
            knots = int(ny * knotspp)
        else:
            knots = 1000
        self.logger.info("Using %d knots for bspline fit" % knots)

        # generate a fit from ref slice points
        bkpt = np.min(xfr) + np.arange(knots+1) * \
            (np.max(xfr) - np.min(xfr)) / knots
        sftr, _ = Bspline.iterfit(xfr[nrefx:-nrefx], yfr[nrefx:-nrefx],
                                  fullbkpt=bkpt)
        yfitr, _ = sftr.value(xfr)

        # generate a blue slice spectrum bspline fit
        if self.action.args.camera == 0:  # Blue
            blueslice = 12
        else:
            blueslice = 11
        blueleft = 60 / xbin
        blueright = 80 / xbin
        qb = [i for i, v in enumerate(slicemap.data.flat) if v == blueslice]
        qblue = [i for i in qb if blueleft <= posmap.data.flat[i] <= blueright]
        xfb = wavemap.data.flat[qblue]
        yfb = newflat.flat[qblue]
        s = np.argsort(xfb)
        xfb = xfb[s]
        yfb = yfb[s]
        bkpt = np.min(xfb) + np.arange(knots+1) * \
            (np.max(xfb) - np.min(xfb)) / knots
        sftb, _ = Bspline.iterfit(xfb[nrefx:-nrefx], yfb[nrefx:-nrefx],
                                  fullbkpt=bkpt)
        yfitb, _ = sftb.value(xfb)

        # generate a red slice spectrum bspline fit
        if self.action.args.camera == 0:  # Blue
            redslice = 23
        else:
            redslice = 0
        redleft = 60 / xbin
        redright = 80 / xbin
        qr = [i for i, v in enumerate(slicemap.data.flat) if v == redslice]
        qred = [i for i in qr if redleft <= posmap.data.flat[i] <= redright]
        xfd = wavemap.data.flat[qred]
        yfd = newflat.flat[qred]
        s = np.argsort(xfd)
        xfd = xfd[s]
        yfd = yfd[s]
        bkpt = np.min(xfd) + np.arange(knots + 1) * \
            (np.max(xfd) - np.min(xfd)) / knots
        sftd, _ = Bspline.iterfit(xfd[nrefx:-nrefx], yfd[nrefx:-nrefx],
                                  fullbkpt=bkpt)
        yfitd, _ = sftd.value(xfd)

        # waves
        minwave = np.min(xfb)
        maxwave = np.max(xfd)
        # are we a twilight flat?
        if twiflat:
            nwaves = int(ny * knotspp)
        else:
            nwaves = 1000
        waves = minwave + (maxwave - minwave) * np.arange(nwaves+1) / nwaves
        if self.config.instrument.plot_level >= 1:
            # output filename stub
            rbfnam = "redblue_%05d_%s_%s_%s" % \
                      (stacked.header['FRAMENO'],
                       self.action.args.illum, self.action.args.grating,
                       self.action.args.ifuname)
            if xbin == 1:
                stride = int(len(xfr) / 8000.)
                if stride <= 0:
                    stride = 1
            else:
                stride = 1
            xrplt = xfr[::stride]
            yrplt = yfitr[::stride]
            yrplt_d = yfr[::stride]
            xbplt = xfb[::stride]
            ybplt = yfitb[::stride]
            ybplt_d = yfb[::stride]
            xdplt = xfd[::stride]
            ydplt = yfitd[::stride]
            ydplt_d = yfd[::stride]
            ymax = np.max(yrplt_d)
            ylim = [0, ymax + ymax * 0.1]
            p = figure(
                title=plab + ' Blue/Red fits',
                x_axis_label='Wave (A)',
                y_axis_label='Flux (e-)',
                y_range=ylim,
                plot_width=self.config.instrument.plot_width,
                plot_height=self.config.instrument.plot_height)
            p.line(xrplt, yrplt, line_color='black', legend_label='Ref')
            p.circle(xrplt, yrplt_d, size=1, line_alpha=0., fill_color='black')
            p.line(xbplt, ybplt, line_color='blue', legend_label='Blue')
            p.circle(xbplt, ybplt_d, size=1, line_alpha=0., fill_color='blue')
            p.line(xdplt, ydplt, line_color='red', legend_label='Red')
            p.circle(xdplt, ydplt_d, size=1, line_alpha=0., fill_color='red')
            bokeh_plot(p, self.context.bokeh_session)
            pause(self.config.instrument, self.logger,
                  what="Flat spectra of the reference slice and of the slices reaching furthest blue and red.",
                  help=("What you see:\n"
                        '  Title "Blue/Red fits", flux (e-) against wavelength.\n'
                        "  Black: reference slice (9), data (dots) and spline fit (line).\n"
                        "  Blue: the slice reaching furthest to the blue (12 for the blue\n"
                        "  camera, 11 for red). Red: the slice reaching furthest to the red\n"
                        "  (23 for the blue camera, 0 for red). These extend the illumination\n"
                        "  spectrum beyond the reference slice's wavelength range.\n"
                        "What to look for:\n"
                        "  Good: each line follows its dots; the curves have the same shape\n"
                        "  where they overlap.\n"
                        "  Bad: a fit that departs from its data, or very different shapes\n"
                        "  in the overlap.\n"
                        "Your options:\n"
                        "  <cr>  continue. Display only.\n"
                        "Default:\n"
                        "  <cr> continues. At plot_level 1 the plot is shown for plot_pause\n"
                        "  seconds; at plot_level 0 it is not drawn."),
                  who="MakeMasterFlat")
            save_plot(p, filename=rbfnam+'.png')

        wavebuffer = 0.1
        minrwave = np.min(xfr)
        maxrwave = np.max(xfr)
        wavebuffer2 = 0.05
        wlb0 = minrwave+(maxrwave-minrwave)*wavebuffer2
        wlb1 = minrwave+(maxrwave-minrwave)*wavebuffer
        wlr0 = minrwave+(maxrwave-minrwave)*(1.-wavebuffer)
        wlr1 = minrwave+(maxrwave-minrwave)*(1.-wavebuffer2)
        qbluefit = [i for i, v in enumerate(waves) if wlb0 < v < wlb1]
        qredfit = [i for i, v in enumerate(waves) if wlr0 < v < wlr1]

        nqb = len(qbluefit)
        nqr = len(qredfit)
        self.logger.info("Wavelength regions: blue = %.1f - %.1f, "
                         "red = %.1f - %.1f" % (wlb0, wlb1, wlr0, wlr1))
        self.logger.info("Fit points: blue = %d, red = %d" % (nqb, nqr))

        if nqb > 0:
            bluefit, _ = sftb.value(waves[qbluefit])
            refbluefit, _ = sftr.value(waves[qbluefit])
            blue_zero_cross = np.nanmin(bluefit) <= 0. or np.nanmin(
                refbluefit) <= 0.
            bluelinfit = np.polyfit(waves[qbluefit], refbluefit/bluefit, 1)
            bluelinfity = bluelinfit[1] + bluelinfit[0] * waves[qbluefit]
        else:
            bluefit = None
            blue_zero_cross = False
            refbluefit = None
            bluelinfit = None
            bluelinfity = None
        if nqr > 0:
            redfit, _ = sftd.value(waves[qredfit])
            refredfit, _ = sftr.value(waves[qredfit])
            red_zero_cross = np.nanmin(redfit) <= 0. or np.nanmin(
                refredfit) <= 0.
            redlinfit = np.polyfit(waves[qredfit], refredfit/redfit, 1)
            redlinfity = redlinfit[1] + redlinfit[0] * waves[qredfit]
        else:
            redfit = None
            red_zero_cross = False
            refredfit = None
            redlinfit = None
            redlinfity = None
        if blue_zero_cross:
            self.logger.info("Blue extension zero crossing detected")
        if red_zero_cross:
            self.logger.info("Red extension zero crossing detected")
        if self.config.instrument.plot_level >= 1:
            if nqb > 1:
                # plot blue fits
                p = figure(
                    title=plab + ' Blue fits',
                    x_axis_label='Wave (A)',
                    y_axis_label='Flux (e-)',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                p.line(waves[qbluefit], refbluefit, line_color='black',
                       legend_label='Ref')
                p.circle(waves[qbluefit], refbluefit, fill_color='black')
                p.line(waves[qbluefit], bluefit, line_color='blue',
                       legend_label='Blue')
                p.circle(waves[qbluefit], bluefit, fill_color='blue')
                bokeh_plot(p, self.context.bokeh_session)
                pause(self.config.instrument, self.logger,
                      what="Blue overlap: reference-slice and blue-slice flat fits.",
                      help=("What you see:\n"
                            '  Title "Blue fits": the spline fits in the overlap near the blue end\n'
                            "  of the reference slice (5-10 % into its range).\n"
                            "  Black line and filled circles: reference slice. Blue line and filled\n"
                            "  circles: the blue-end slice. (Circles have a blue outline.)\n"
                            "What to look for:\n"
                            "  Good: both curves smooth, with similar shapes; their ratio (next\n"
                            "  plot) scales the blue slice to join the reference.\n"
                            "  Bad: noisy or diverging curves: the blue extension of the flat will\n"
                            "  have a step or slope error.\n"
                            "Your options:\n"
                            "  <cr>  continue. Display only.\n"
                            "Default:\n"
                            "  <cr> continues. At plot_level 1 the plot is shown for plot_pause\n"
                            "  seconds; at plot_level 0 it is not drawn."),
                      who="MakeMasterFlat")
                # plot blue ratios
                p = figure(
                    title=plab + ' Blue ratios',
                    x_axis_label='Wave (A)',
                    y_axis_label='Ratio',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                p.line(waves[qbluefit], refbluefit/bluefit, line_color='black',
                       legend_label='Ref')
                p.circle(waves[qbluefit], refbluefit/bluefit,
                         fill_color='black')
                p.line(waves[qbluefit], bluelinfity,
                       line_color='blue', legend_label='Blue')
                p.circle(waves[qbluefit], bluelinfity,
                         fill_color='blue')
                bokeh_plot(p, self.context.bokeh_session)
                pause(self.config.instrument, self.logger,
                      what="Blue overlap: ratio of reference to blue-slice fit, and its linear fit.",
                      help=("What you see:\n"
                            '  Title "Blue ratios". Black line and filled circles: reference-slice\n'
                            "  fit divided by blue-slice fit. Blue line and filled circles: linear\n"
                            "  fit to that ratio, used to scale the blue-slice data blueward of the\n"
                            "  reference slice. (Circles have a blue outline.)\n"
                            "What to look for:\n"
                            "  Good: a smooth, nearly constant ratio, well followed by the line.\n"
                            "  Bad: a strongly curved or noisy ratio; the linear scaling will not\n"
                            "  match it. (If either fit crosses zero only an offset is applied.)\n"
                            "Your options:\n"
                            "  <cr>  continue. Display only.\n"
                            "Default:\n"
                            "  <cr> continues. At plot_level 1 the plot is shown for plot_pause\n"
                            "  seconds; at plot_level 0 it is not drawn."),
                      who="MakeMasterFlat")
            if nqr > 1:
                # plot red fits
                p = figure(
                    title=plab + ' Red fits',
                    x_axis_label='Wave (A)',
                    y_axis_label='Flux (e-)',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                p.line(waves[qredfit], refredfit, line_color='black',
                       legend_label='Ref')
                p.circle(waves[qredfit], refredfit, fill_color='black')
                p.line(waves[qredfit], redfit, line_color='red',
                       legend_label='Red')
                p.circle(waves[qredfit], redfit, fill_color='red')
                bokeh_plot(p, self.context.bokeh_session)
                pause(self.config.instrument, self.logger,
                      what="Red overlap: reference-slice and red-slice flat fits.",
                      help=("What you see:\n"
                            '  Title "Red fits": the spline fits in the overlap near the red end of\n'
                            "  the reference slice (90-95 % into its range).\n"
                            "  Black line and filled circles: reference slice. Red line and filled\n"
                            "  circles: the red-end slice. (Circles have a blue outline.)\n"
                            "What to look for:\n"
                            "  Good: both curves smooth, with similar shapes; their ratio (next\n"
                            "  plot) scales the red slice to join the reference.\n"
                            "  Bad: noisy or diverging curves: the red extension of the flat will\n"
                            "  have a step or slope error.\n"
                            "Your options:\n"
                            "  <cr>  continue. Display only.\n"
                            "Default:\n"
                            "  <cr> continues. At plot_level 1 the plot is shown for plot_pause\n"
                            "  seconds; at plot_level 0 it is not drawn."),
                      who="MakeMasterFlat")
                # plot red ratios
                p = figure(
                    title=plab + ' Red ratios',
                    x_axis_label='Wave (A)',
                    y_axis_label='Ratio',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                p.line(waves[qredfit], refredfit/redfit, line_color='black',
                       legend_label='Ref')
                p.circle(waves[qredfit], refredfit/redfit, fill_color='black')
                p.line(waves[qredfit], redlinfity,
                       line_color='red', legend_label='Red')
                p.circle(waves[qredfit], redlinfity, fill_color='red')
                bokeh_plot(p, self.context.bokeh_session)
                pause(self.config.instrument, self.logger,
                      what="Red overlap: ratio of reference to red-slice fit, and its linear fit.",
                      help=("What you see:\n"
                            '  Title "Red ratios". Black line and filled circles: reference-slice\n'
                            "  fit divided by red-slice fit. Red line and filled circles: linear fit\n"
                            "  to that ratio, used to scale the red-slice data redward of the\n"
                            "  reference slice. (Circles have a blue outline.)\n"
                            "What to look for:\n"
                            "  Good: a smooth, nearly constant ratio, well followed by the line.\n"
                            "  Bad: a strongly curved or noisy ratio; the linear scaling will not\n"
                            "  match it. (If either fit crosses zero only an offset is applied.)\n"
                            "Your options:\n"
                            "  <cr>  continue. Display only.\n"
                            "Default:\n"
                            "  <cr> continues. At plot_level 1 the plot is shown for plot_pause\n"
                            "  seconds; at plot_level 0 it is not drawn."),
                      who="MakeMasterFlat")

        # at this point we are going to try to merge the points
        self.logger.info("Correcting points outside %.1f - %.1f A"
                         % (minrwave, maxrwave))
        qselblue = [i for i, v in enumerate(xfb) if v <= minrwave]
        qselred = [i for i, v in enumerate(xfd) if v >= maxrwave]
        nqsb = len(qselblue)
        nqsr = len(qselred)
        blue_all_tie = yfitr[0]
        red_all_tie = yfitr[-1]
        self.logger.info("Blue/Red ref tie values: %.3f, %.3f"
                         % (blue_all_tie, red_all_tie))

        if nqsb > 0:
            self.logger.info("Blue ext tie value: %.3f" % yfb[qselblue[-1]])
            if blue_zero_cross:
                blue_offset = yfb[qselblue[-1]] - blue_all_tie
                bluefluxes = [yfb[i] - blue_offset for i in qselblue]
                self.logger.info("Blue zero crossing, only applying offset")
            else:
                blue_offset = yfb[qselblue[-1]] * \
                              (bluelinfit[1]+bluelinfit[0]*xfb[qselblue[-1]]) \
                              - blue_all_tie
                bluefluxes = [yfb[i] * (bluelinfit[1]+bluelinfit[0]*xfb[i])
                              - blue_offset for i in qselblue]
                self.logger.info("Blue linear ratio fit scaling applied")
            self.logger.info("Blue offset of %.2f applied" % blue_offset)
        else:
            bluefluxes = None
        if nqsr > 0:
            self.logger.info("Red ext tie value: %.3f" % yfd[qselred[0]])
            if red_zero_cross:
                red_offset = yfd[qselred[0]] - red_all_tie
                redfluxes = [yfd[i] - red_offset for i in qselred]
                self.logger.info("Red zero crossing, only applying offset")
            else:
                red_offset = yfd[qselred[0]] * \
                             (redlinfit[1]+redlinfit[0]*xfd[qselred[0]]) \
                             - red_all_tie
                redfluxes = [yfd[i] * (redlinfit[1]+redlinfit[0]*xfd[i])
                             - red_offset for i in qselred]
                self.logger.info("Red linear ratio fit scaling applied")
            self.logger.info("Red offset of %.2f applied" % red_offset)
        else:
            redfluxes = None
        allx = xfr
        allfx = xfr[nrefx:-nrefx]
        ally = yfr[nrefx:-nrefx]
        if nqsb > 0:
            bluex = xfb[qselblue]
            allx = np.append(bluex, allx)
            allfx = np.append(bluex[nrefx:], allfx)
            ally = np.append(bluefluxes[nrefx:], ally)
        if nqsr > 0:
            redx = xfd[qselred]
            allx = np.append(allx, redx)
            allfx = np.append(allfx, redx[:-nrefx])
            ally = np.append(ally, redfluxes[:-nrefx])
        s = np.argsort(allx)
        allx = allx[s]
        s = np.argsort(allfx)
        allfx = allfx[s]
        ally = ally[s]

        bkpt = np.min(allx) + np.arange(knots+1) * \
            (np.max(allx) - np.min(allx)) / knots
        sftall, _ = Bspline.iterfit(allfx, ally, fullbkpt=bkpt)
        yfitall, _ = sftall.value(allx)

        if self.config.instrument.plot_level >= 1:
            # output filename stub
            fltfnam = "flat_%05d_%s_%s_%s" % \
                      (stacked.header['FRAMENO'],
                       self.action.args.illum, self.action.args.grating,
                       self.action.args.ifuname)
            if xbin == 1:
                stride = int(len(allx) / 8000.)
            else:
                stride = 1
            xplt = allfx[::stride]
            yplt = ally[::stride]
            fxplt = allx[::stride]
            fplt = yfitall[::stride]
            yran = [0, np.nanmax(ally)]
            p = figure(
                title=plab + ' Master Illumination',
                x_axis_label='Wave (A)',
                y_axis_label='Flux (e-)',
                y_range=yran,
                plot_width=self.config.instrument.plot_width,
                plot_height=self.config.instrument.plot_height)
            p.circle(xplt, yplt, size=1, line_alpha=0., fill_color='black',
                     legend_label='Data')
            p.line(fxplt, fplt, line_color='red', legend_label='Fit',
                   line_width=2)
            p.line([minrwave, minrwave], yran, line_color='orange',
                   legend_label='Cor lim')
            p.line([maxrwave, maxrwave], yran, line_color='orange')
            p.legend.location = "top_left"
            bokeh_plot(p, self.context.bokeh_session)
            pause(self.config.instrument, self.logger,
                  what="Master illumination spectrum: combined flat data over the full range, with fit.",
                  help=("What you see:\n"
                        '  Title "Master Illumination", flux (e-) against wavelength.\n'
                        "  Black dots: reference-slice data plus the scaled blue- and red-end\n"
                        "  slice data. Red thick line: the spline fit, the illumination\n"
                        '  spectrum the flat is divided by. Orange verticals ("Cor lim"): the\n'
                        "  reference slice's range; outside it the end slices are used.\n"
                        "What to look for:\n"
                        "  Good: the red curve follows the data over the whole range, with no\n"
                        "  jump at the orange lines and no ringing at the ends.\n"
                        "  Bad: a jump or kink at an orange line (end-slice scaling), or wiggles\n"
                        "  where the data are smooth: these go straight into the master flat.\n"
                        "Your options:\n"
                        "  <cr>  continue. Display only.\n"
                        "Default:\n"
                        "  <cr> continues. At plot_level 1 the plot is shown for plot_pause\n"
                        "  seconds; at plot_level 0 it is not drawn."),
                  who="MakeMasterFlat")
            save_plot(p, filename=fltfnam+".png")

        # OK, Now we have extended to the full range... so... we are going to
        # make a ratio flat!
        comflat = np.zeros(newflat.shape, dtype=float)
        qz = [i for i, v in enumerate(wavemap.data.flat) if v >= 0]

        comvals = sftall.value(wavemap.data.flat[qz])

        comflat.flat[qz] = comvals
        ratio = np.zeros(newflat.shape, dtype=float)
        qzer = [i for i, v in enumerate(newflat.flat) if v != 0]
        ratio.flat[qzer] = comflat.flat[qzer] / newflat.flat[qzer]

        # trim negative points
        qq = [i for i, v in enumerate(ratio.flat) if v < 0]
        if len(qq) > 0:
            ratio.flat[qq] = 0.0

        # trim the high points near edges of slice
        qq = [i for i, v in enumerate(ratio.flat) if v >= 3. and
              (posmap.data.flat[i] <= 4/xbin or
               posmap.data.flat[i] >= 136/xbin)]
        if len(qq) > 0:
            ratio.flat[qq] = 0.0

        # don't correct low signal points
        qq = [i for i, v in enumerate(newflat.flat) if v < 30.]
        if len(qq) > 0:
            ratio.flat[qq] = 1.0

        # get master flat output name
        mfname = strip_fname(stack_list[0]) + '_' + suffix + '.fits'

        log_string = MakeMasterFlat.__module__
        stacked.header['IMTYPE'] = self.action.args.new_type
        stacked.header['HISTORY'] = log_string
        stacked.header['MASTFLAT'] = (True, 'master flat image?')
        stacked.header['WAVMAPF'] = wmf
        stacked.header['SLIMAPF'] = slf
        stacked.header['POSMAPF'] = pof

        # store flat in output frame
        stacked.data = ratio

        # output master flat
        kcwi_fits_writer(stacked, output_file=mfname,
                         output_dir=self.config.instrument.output_directory)
        self.context.proctab.update_proctab(frame=stacked, suffix=suffix,
                                            newtype=self.action.args.new_type,
                                            filename=stack_list[0]) ### HERE
        self.context.proctab.write_proctab(tfil=self.config.instrument.procfile)
        # self.action.args.name = stacked.header['OFNAME']
        # self.action.args.name = mfname

        self.logger.info(log_string)
        return self.action.args

    # END: class MakeMasterFlat()
