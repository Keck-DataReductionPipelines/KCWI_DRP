from keckdrpframework.primitives.base_img import BaseImg
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer
from kcwidrp.core.bokeh_plotting import bokeh_plot
from bokeh.plotting import figure
from bokeh.models import Range1d

import os
import time
import numpy as np
from scipy.signal.windows import boxcar
import scipy as sp
from scipy.signal import find_peaks
from pydl.pydlutils import bspline


def bm_ledge_position(cwave):
    fit = [0.240742, 4044.56]
    return fit[1] + fit[0] * cwave


class MakeMasterFlat(BaseImg):
    """Stack flat images and make master flat image"""

    def __init__(self, action, context):
        BaseImg.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if we can create a master flat based on the processing table
        :return:
        """
        # get list of master flats
        self.logger.info("Checking precondition for MakeMasterFlat")
        target_type = 'MFLAT'
        tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                             target_type=target_type,
                                             nearest=True)
        if len(tab) > 0:
            self.logger.info(f"already have {len(tab)} master flats, "
                             f" expecting 0")
            return False
        else:
            self.stack_list = self.context.proctab.n_proctab(
                frame=self.action.args.ccddata,
                target_type=self.action.args.stack_type,
                target_group=self.action.args.groupid)
            self.logger.info(f"pre condition got {len(self.stack_list)},"
                             f" expecting 1")
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

        suffix = self.action.args.new_type.lower()
        insuff = self.action.args.stack_type.lower()

        stack_list = list(self.stack_list['OFNAME'])

        if len(stack_list) <= 0:
            self.logger.warning("No flats found!")
            return self.action.args

        # get root for maps
        tab = self.context.proctab.n_proctab(
            frame=self.action.args.ccddata, target_type='ARCLAMP',
            target_group=self.action.args.groupid)
        if len(tab) <= 0:
            self.logger.error("Geometry not solved!")
            return self.action.args

        mroot = tab['OFNAME'][0].split('.fits')[0]

        # Wavelength map image
        wmf = mroot + '_wavemap.fits'
        self.logger.info("Reading image: %s" % wmf)
        wavemap = kcwi_fits_reader(
            os.path.join(os.path.dirname(self.action.args.name), 'redux',
                         wmf))[0]

        # Slice map image
        slf = mroot + '_slicemap.fits'
        self.logger.info("Reading image: %s" % slf)
        slicemap = kcwi_fits_reader(
            os.path.join(os.path.dirname(self.action.args.name), 'redux',
                         slf))[0]

        # Position map image
        pof = mroot + '_posmap.fits'
        self.logger.info("Reading image: %s" % pof)
        posmap = kcwi_fits_reader(
            os.path.join(os.path.dirname(self.action.args.name), 'redux',
                         pof))[0]

        # Read in stacked flat image
        stname = stack_list[0].split('.')[0] + '_' + insuff + '.fits'

        self.logger.info("Reading image: %s" % stname)
        stacked = kcwi_fits_reader(
            os.path.join(os.path.dirname(self.action.args.name), 'redux',
                         stname))[0]

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

        # vignetted slice position range
        fitl = int(4/xbin)
        fitr = int(24/xbin)

        # un-vignetted slice position range
        flatl = int(34/xbin)
        flatr = int(72/xbin)

        # flat fitting slice position range
        ffleft = int(10/xbin)
        ffright = int(70/xbin)

        buffer = 6.0/float(xbin)

        # reference slice
        refslice = 9
        allidx = np.arange(int(140/xbin))
        newflat = stacked.data.copy()

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
            dw = (waves[1] - waves[0]) / 30.0
            wavemin = (waves[0]+waves[1]) / 2.0 - dw
            wavemax = (waves[0]+waves[1]) / 2.0 + dw
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
                p = figure(title=self.action.args.plotlabel + ' WAVE SLOPE FIT',
                           x_axis_label='wave px',
                           y_axis_label='counts',
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.circle(wwflat, ywflat, legend_label="Data")
                p.line(wflat, wslfit, line_color='red', line_width=3,
                       legend_label="Fit")
                bokeh_plot(p, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)
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
                p = figure(title=self.action.args.plotlabel + ' Vignetting',
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
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)

            # figure out where the correction applies
            qcor = [i for i, v in enumerate(posmap.data.flat)
                    if 0 <= v <= (xinter-buffer)]
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
                p = figure(title=self.action.args.plotlabel + ' Buffer Region',
                           x_axis_label='Slice Pos (px)',
                           y_axis_label='Ratio',
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.circle(xbuff, ybuff)
                p.line(xbuff, np.polyval(buffit, xbuff), line_color='red')
                bokeh_plot(p, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)
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
        qref = [i for i in q if ffleft <= posmap.data.flat[i] <= ffright]
        xfr = wavemap.data.flat[qref]
        yfr = newflat.flat[qref]
        # sort on wavelength
        s = np.argsort(xfr)
        xfr = xfr[s]
        yfr = yfr[s]

        wavegood0 = wavemap.header['WAVGOOD0']
        wavegood1 = wavemap.header['WAVGOOD1']

        # correction for BM where we see a ledge
        if 'BM' in self.action.args.grating:
            ledge_wave = bm_ledge_position(self.action.args.cwave)

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
                win = boxcar(150)
                smyledge = sp.signal.convolve(yledge,
                                              win, mode='same') / sum(win)
                ylmax = np.max(yledge)
                ylmin = np.min(yledge)
                fpoints = np.arange(0, 100) / 100. * 50 + (ledge_wave-25)
                ledgefit, ledgemsk = bspline.iterfit(np.asarray(xledge),
                                                     smyledge, fullbkpt=fpoints,
                                                     upper=1, lower=1)
                ylfit, _ = ledgefit.value(np.asarray(fpoints))
                deriv = -(np.roll(ylfit, 1) - np.roll(ylfit, -1)) / 2.0
                deriv = deriv[4:-4]
                xvals = fpoints[4:-4]
                peaks, _ = find_peaks(deriv, height=100)
                if len(peaks) != 1:
                    raise ValueError
                ipk = peaks[0]
                apk = xvals[ipk]
                if self.config.instrument.plot_level >= 3:
                    p = figure(
                        title=self.action.args.plotlabel + ' Peak of ledge',
                        x_axis_label='Wave (A)',
                        y_axis_label='Value',
                        plot_width=self.config.instrument.plot_width,
                        plot_height=self.config.instrument.plot_height)
                    p.circle(xvals, deriv, legend_label='Data')
                    p.line([apk, apk], [-50, 200], line_color='red',
                           legend_label='Pk')
                    bokeh_plot(p, self.context.bokeh_session)
                    if self.config.instrument.plot_level >= 2:
                        input("Next? <cr>: ")
                    else:
                        time.sleep(self.config.instrument.plot_pause)
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
                        title=self.action.args.plotlabel + ' BM Ledge Region',
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
                    if self.config.instrument.plot_level >= 2:
                        input("Next? <cr>: ")
                    else:
                        time.sleep(self.config.instrument.plot_pause)
        # END: handling BM grating ledge

        # if we are fitting a twilight flat, treat it like a sky image with a
        # larger number of knots
        if twiflat:
            knots = int(ny * knotspp)
        else:
            knots = 100
        self.logger.info("Using %d knots for bspline fit" % knots)
        bkpt = np.min(xfr) + np.arange(knots+1) * \
            (np.max(xfr) - np.min(xfr)) / knots
        sftr, _ = bspline.iterfit(xfr, yfr, fullbkpt=bkpt)
        yfitr, _ = sftr.value(xfr)

        # generate a blue slice spectrum bspline fit
        blueslice = 12
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
        sftb, _ = bspline.iterfit(xfb, yfb, fullbkpt=bkpt)
        yfitb, _ = sftb.value(xfb)

        # generate a red slice spectrum bspline fit
        redslice = 23
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
        sftd, _ = bspline.iterfit(xfd, yfd, fullbkpt=bkpt)
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
            p = figure(
                title=self.action.args.plotlabel + ' Blue/Red fits',
                x_axis_label='Wave (A)',
                y_axis_label='Flux (e-)',
                plot_width=self.config.instrument.plot_width,
                plot_height=self.config.instrument.plot_height)
            p.line(xfr, yfitr, line_color='black', legend_label='Ref')
            p.line(xfb, yfitb, line_color='blue', legend_label='Blue')
            p.line(xfd, yfitd, line_color='red', legend_label='Red')
            bokeh_plot(p, self.context.bokeh_session)
            if self.config.instrument.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.config.instrument.plot_pause)

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
            bluelinfit = np.polyfit(waves[qbluefit], refbluefit/bluefit, 1)
        else:
            bluefit = None
            refbluefit = None
            bluelinfit = None
        if nqr > 0:
            redfit, _ = sftd.value(waves[qredfit])
            refredfit, _ = sftr.value(waves[qredfit])
            redlinfit = np.polyfit(waves[qredfit], refredfit/redfit, 1)
        else:
            redfit = None
            refredfit = None
            redlinfit = None
        if self.config.instrument.plot_level >= 1:
            if nqb > 1:
                # plot blue fits
                p = figure(
                    title=self.action.args.plotlabel + ' Blue fits',
                    x_axis_label='Wave (A)',
                    y_axis_label='Flux (e-)',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                p.line(waves[qbluefit], refbluefit, line_color='black',
                       legend_label='Ref')
                p.line(waves[qbluefit], bluefit, line_color='blue',
                       legend_label='Blue')
                bokeh_plot(p, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)
                # plot blue ratios
                p = figure(
                    title=self.action.args.plotlabel + ' Blue ratios',
                    x_axis_label='Wave (A)',
                    y_axis_label='Ratio',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                p.line(waves[qbluefit], refbluefit/bluefit, line_color='black',
                       legend_label='Ref')
                p.line(waves[qbluefit],
                       bluelinfit[1]+bluelinfit[0]*waves[qbluefit],
                       line_color='blue', legend_label='Blue')
                bokeh_plot(p, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)
            if nqr > 1:
                # plot red fits
                p = figure(
                    title=self.action.args.plotlabel + ' Red fits',
                    x_axis_label='Wave (A)',
                    y_axis_label='Flux (e-)',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                p.line(waves[qredfit], refredfit, line_color='black',
                       legend_label='Ref')
                p.line(waves[qredfit], redfit, line_color='red',
                       legend_label='Red')
                bokeh_plot(p, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)
                # plot red ratios
                p = figure(
                    title=self.action.args.plotlabel + ' Red ratios',
                    x_axis_label='Wave (A)',
                    y_axis_label='Ratio',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                p.line(waves[qredfit], refredfit/redfit, line_color='black',
                       legend_label='Ref')
                p.line(waves[qredfit],
                       redlinfit[1]+redlinfit[0]*waves[qredfit],
                       line_color='red', legend_label='Red')
                bokeh_plot(p, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)

        # at this point we are going to try to merge the points
        qselred = [i for i, v in enumerate(xfd) if v >= maxrwave]
        qselblue = [i for i, v in enumerate(xfb) if v <= minrwave]
        nqsr = len(qselred)
        nqsb = len(qselblue)

        if nqsr > 0:
            redfluxes = [yfd[i] * (redlinfit[1]+redlinfit[0]*xfd[i])
                         for i in qselred]
        else:
            redfluxes = None
        if nqsb > 0:
            bluefluxes = [yfb[i] * (bluelinfit[1]+bluelinfit[0]*xfb[i])
                          for i in qselblue]
        else:
            bluefluxes = None
        allx = xfr
        ally = yfr
        if nqsb > 0:
            allx = np.append(xfb[qselblue], allx)
            ally = np.append(bluefluxes, ally)
        if nqsr > 0:
            allx = np.append(allx, xfd[qselred])
            ally = np.append(ally, redfluxes)
        s = np.argsort(allx)
        allx = allx[s]
        ally = ally[s]

        bkpt = np.min(allx) + np.arange(knots+1) * \
            (np.max(allx) - np.min(allx)) / knots
        sftall, _ = bspline.iterfit(allx, ally, fullbkpt=bkpt)
        yfitall, _ = sftall.value(allx)

        if self.config.instrument.plot_level >= 1:
            p = figure(
                title=self.action.args.plotlabel + ' Master Illumination',
                x_axis_label='Wave (A)',
                y_axis_label='Flux (e-)',
                plot_width=self.config.instrument.plot_width,
                plot_height=self.config.instrument.plot_height)
            p.circle(allx, ally, size=1, line_alpha=0., fill_color='purple',
                     legend_label='Data')
            p.line(allx, yfitall, line_color='red', legend_label='Fit')
            p.circle(xfr, yfr, size=1, line_alpha=0., fill_color='black',
                     legend_label='Ref Data')
            p.line(xfr, yfitr, line_color='green', legend_label='Ref Fit')
            bokeh_plot(p, self.context.bokeh_session)
            if self.config.instrument.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.config.instrument.plot_pause)

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
        mfname = stack_list[0].split('.fits')[0] + '_' + suffix + '.fits'

        log_string = MakeMasterFlat.__module__ + "." + \
            MakeMasterFlat.__qualname__
        stacked.header['IMTYPE'] = self.action.args.new_type
        stacked.header['HISTORY'] = log_string
        stacked.header['MASTFLAT'] = (True, 'master flat image?')
        stacked.header['WAVMAPF'] = wmf
        stacked.header['SLIMAPF'] = slf
        stacked.header['POSMAPF'] = pof

        stacked.data = ratio

        # output master flat
        kcwi_fits_writer(stacked, output_file=mfname)
        self.context.proctab.update_proctab(frame=stacked, suffix=suffix,
                                            newtype=self.action.args.new_type)
        self.context.proctab.write_proctab()

        self.logger.info(log_string)
        return self.action.args

    # END: class MakeMasterFlat()
