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

        buffer = 5.0/float(xbin)

        # reference slice
        refslice = 9
        fflice = refslice
        ffslice2 = refslice
        sm = 25
        allidx = np.arange(int(140/xbin))
        newflat = stacked.data.copy()

        # correct vignetting if we are using internal flats
        if internal:
            # get good region for fitting
            waves = wavemap.data.compress((wavemap.data > 0.).flat)
            waves = [waves.min(), waves.max()]
            dw = (waves[1] - waves[0]) / 30.0
            wavemin = (waves[0]+waves[1]) / 2.0 - dw
            wavemax = (waves[0]+waves[1]) / 2.0 + dw
            print(wavemin, wavemax)

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

            qflat = [i for i, v in enumerate(xflat) if flatl <= v <= flatr]
            xflat = [xflat[i] for i in qflat]
            yflat = [yflat[i] for i in qflat]
            wflat = [wflat[i] for i in qflat]
            print(np.min(xflat), np.max(xflat))
            sw = np.argsort(wflat)
            ywflat = [yflat[i] for i in sw]
            wwflat = [wflat[i] for i in sw]
            wavelinfit = np.polyfit(wwflat, ywflat, 2)
            print(wavelinfit)
            yfit = np.polyval(wavelinfit, wflat)

            if self.config.instrument.plot_level >= 1:
                p = figure(title=self.action.args.plotlabel + ' WAVE SLOPE FIT',
                           x_axis_label='wave px',
                           y_axis_label='counts',
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.circle(wwflat, ywflat, legend_label="Data")
                p.line(wflat, yfit, line_color='red', line_width=3,
                       legend_label="Fit")
                bokeh_plot(p, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)

            yflat = yflat / np.polyval(wavelinfit, wflat)
            ss = np.argsort(xflat)
            xflat = [xflat[i] for i in ss]
            yflat = [yflat[i] for i in ss]
            wflat = [wflat[i] for i in ss]
            resflat = np.polyfit(xflat, yflat, 1)

            # select the point we will fit for the vignetting
            xfit = [posmap.data.flat[i] for i in qq]
            yfit = [stacked.data.flat[i] for i in qq]
            wflat = [wavemap.data.flat[i] for i in qq]
            yfit = yfit / np.polyval(wavelinfit, wflat)

            # fit the vignetted region
            qfit = [i for i, v in enumerate(xfit) if fitl <= v <= fitr]
            xfit = [xfit[i] for i in qfit]
            yfit = [yfit[i] for i in qfit]
            s = np.argsort(xfit)
            xfit = [xfit[i] for i in s]
            yfit = [yfit[i] for i in s]
            resfit = np.polyfit(xfit, yfit, 1)

            if self.config.instrument.plot_level >= 1:
                p = figure(title=self.action.args.plotlabel + ' Vignetting',
                           x_axis_label='Slice Pos (px)',
                           y_axis_label='Ratio',
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.circle(posmap.data.flat[qq],
                         stacked.data.flat[qq] / np.polyval(
                             wavelinfit, wavemap.data.flat[qq]))
                p.line(allidx, resfit[1] + resfit[0]*allidx, line_color='purple')
                p.line(allidx, resflat[1] + resflat[0]*allidx, line_color='red')
                p.line([fitl, fitl], [0.5, 1.5], line_color='blue')
                p.line([fitr, fitr], [0.5, 1.5], line_color='blue')
                p.line([flatl, flatl], [0.5, 1.5], line_color='green')
                p.line([flatr, flatr], [0.5, 1.5], line_color='green')
                bokeh_plot(p, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)
            # compute the intersection
            xinter = -(resflat[1]-resfit[1])/(resflat[0]-resfit[0])
            # figure out where the correction applies
            qinter = [i for i, v in enumerate(posmap.data.flat)
                      if 0 <= v <= (xinter-buffer)]
            # apply the correction!
            for i in qinter:
                newflat.flat[i] = (resflat[1]+resflat[0]*posmap.data.flat[i]) \
                                / (resfit[1]+resfit[0]*posmap.data.flat[i]) * \
                                stacked.data.flat[i]
            # now deal with the intermediate (buffer) region
            qspline = [i for i, v in enumerate(posmap.data.flat)
                       if (xinter-buffer) <= v <= (xinter+buffer)]
            xspline = [posmap.data.flat[i] for i in qspline]
            posmin = np.min(xspline)
            posmax = np.max(xspline)
            valuemin = (resflat[1]+resflat[0]*posmin) / \
                       (resfit[1]+resfit[0]*posmin)
            valuemax = 1
            cs = CubicSpline([posmin, posmax], [valuemin, valuemax])
            for i in qspline:
                newflat.flat[i] = cs(posmap.data.flat[i]) * newflat.flat[i]

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
