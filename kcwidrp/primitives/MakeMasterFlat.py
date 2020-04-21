from keckdrpframework.primitives.base_img import BaseImg
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer
from kcwidrp.core.bokeh_plotting import bokeh_plot
from bokeh.plotting import figure
import os
import time
import numpy as np
from scipy.interpolate import CubicSpline


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

        # get image size
        nx = stacked.header['NAXIS1']
        ny = stacked.header['NAXIS2']

        # get binning
        xbin = self.action.args.xbinsize
        ybin = self.action.args.ybinsize

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
        fflice = refslice
        ffslice2 = refslice
        sm = 25
        allidx = np.arange(int(140/xbin))
        newflat = stacked.data.copy()

        # correct vignetting if we are using internal flats
        if internal:
            self.logger.info("Internal flats require vignetting correction")
            # get good region for fitting
            waves = wavemap.data.compress((wavemap.data > 0.).flat)
            waves = [waves.min(), waves.max()]
            dw = (waves[1] - waves[0]) / 30.0
            wavemin = (waves[0]+waves[1]) / 2.0 - dw
            wavemax = (waves[0]+waves[1]) / 2.0 + dw
            self.logger.info("Using %.1f - %.1f A of slice %d" % (wavemin,
                                                                  wavemax,
                                                                  refslice))
            # get reference slice data
            q = [i for i, v in enumerate(slicemap.data.flat) if v == refslice]
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
            # fit wavelength slope
            wavelinfit = np.polyfit(wwflat, ywflat, 2)
            wslfit = np.polyval(wavelinfit, wflat)
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
            yfit = yfit / np.polyval(wavelinfit, wflat)

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
            # compute the intersection
            xinter = -(resflat[1] - resfit[1]) / (resflat[0] - resfit[0])
            # plot slice profile and fits
            if self.config.instrument.plot_level >= 1:
                p = figure(title=self.action.args.plotlabel + ' Vignetting',
                           x_axis_label='Slice Pos (px)',
                           y_axis_label='Ratio',
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.circle(posmap.data.flat[qq],
                         stacked.data.flat[qq] / np.polyval(
                             wavelinfit, wavemap.data.flat[qq]))
                p.line(allidx, resfit[1] + resfit[0]*allidx,
                       line_color='purple')
                p.line(allidx, resflat[1] + resflat[0]*allidx, line_color='red')
                p.line([fitl, fitl], [0.5, 1.5], line_color='blue')
                p.line([fitr, fitr], [0.5, 1.5], line_color='blue')
                p.line([flatl, flatl], [0.5, 1.5], line_color='green')
                p.line([flatr, flatr], [0.5, 1.5], line_color='green')
                p.line([xinter-buffer, xinter-buffer], [0.5, 1.5],
                       line_color='black')
                p.line([xinter + buffer, xinter + buffer], [0.5, 1.5],
                       line_color='black')
                p.line([xinter, xinter], [0.5, 1.5], line_color='red')
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
                                                       wavemap.data.flat[i])
                     for i in qbff]
            # sort on slice position
            sp = np.argsort(xbuff)
            xbuff = [xbuff[i] for i in sp]
            ybuff = [ybuff[i] for i in sp]
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

        stacked.data = newflat
        # get master flat output name
        mfname = stack_list[0].split('.fits')[0] + '_' + suffix + '.fits'

        log_string = MakeMasterFlat.__module__ + "." + \
            MakeMasterFlat.__qualname__
        stacked.header['IMTYPE'] = self.action.args.new_type
        stacked.header['HISTORY'] = log_string
        stacked.header['WAVMAPF'] = wmf
        stacked.header['SLIMAPF'] = slf
        stacked.header['POSMAPF'] = pof

        # output master flat
        kcwi_fits_writer(stacked, output_file=mfname)
        self.context.proctab.update_proctab(frame=stacked, suffix=suffix,
                                            newtype=self.action.args.new_type)
        self.context.proctab.write_proctab()

        self.logger.info(log_string)
        return self.action.args

    # END: class MakeMasterFlat()
