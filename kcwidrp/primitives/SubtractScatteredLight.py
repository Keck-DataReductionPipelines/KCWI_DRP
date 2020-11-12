from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_plotting import save_plot

from bokeh.plotting import figure
import numpy as np
from scipy.signal import savgol_filter
import time


class SubtractScatteredLight(BasePrimitive):
    """Subtract scattered light between slices"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):

        # Header keyword to update
        key = 'SCATSUB'
        keycom = "was scattered light subtracted?"

        # Skip if nod-and-shuffle
        if self.action.args.nasmask:
            self.logger.info("NAS Mask: skipping scattered light subtraction")
            self.action.args.ccddata.header[key] = (False, keycom)
        elif self.config.instrument.skipscat:
            self.logger.info("Skipping scattered light subtraction by request")
            self.action.args.ccddata.header[key] = (False, keycom)
        else:
            # Binning
            ybin = self.action.args.ybinsize
            # Get size of image
            siz = self.action.args.ccddata.data.shape
            # Get x range for scattered light
            x0 = int(siz[1] / 2 - 180 / self.action.args.xbinsize)
            x1 = int(siz[1] / 2 + 180 / self.action.args.xbinsize)
            # Get y limits
            y0 = 0
            # y1 = int(siz[0] / 2 - 1)
            # y2 = y1 + 1
            y3 = siz[0]
            # print("x limits: %d, %d, y limits: %d, %d" % (x0, x1, y0, y3))
            # Y data values
            yvals = np.nanmedian(self.action.args.ccddata.data[y0:y3, x0:x1],
                                 axis=1)
            # Fix extreme values
            yvals[0] = np.nanmedian(yvals[1:10])
            yvals[-1] = np.nanmedian(yvals[-11:-2])
            # X data values
            xvals = np.arange(len(yvals), dtype=np.float)
            # filter window
            fwin = 151
            scat = savgol_filter(yvals, fwin, 3)
            signal_to_noise = np.mean(scat) / np.nanstd(yvals - scat)
            if signal_to_noise < 25.:
                if signal_to_noise < 5:
                    fwin = 501
                else:
                    fwin = 303
                scat = savgol_filter(yvals, fwin, 3)
                signal_to_noise = np.mean(scat) / np.nanstd(yvals - scat)
            self.logger.info("Smoothing scattered light with window of %d px"
                             % fwin)
            self.logger.info("Mean signal to noise = %.2f" % signal_to_noise)
            if self.config.instrument.plot_level >= 1:
                # output filename stub
                scfnam = "scat_%05d_%s_%s_%s" % \
                         (self.action.args.ccddata.header['FRAMENO'],
                          self.action.args.illum, self.action.args.grating,
                          self.action.args.ifuname)
                # plot
                p = figure(title=self.action.args.plotlabel +
                           "SCATTERED LIGHT FIT",
                           x_axis_label='y pixel', y_axis_label='e-',
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.circle(xvals, yvals, legend_label="Scat")
                p.line(xvals, scat, color='red', line_width=3,
                       legend_label="Smooth")
                p.legend.location = "top_left"
                bokeh_plot(p, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)
                save_plot(p, filename=scfnam+".png")
            # Subtract scattered light
            self.logger.info("Starting scattered light subtraction")
            for ix in range(0, siz[1]):
                self.action.args.ccddata.data[y0:y3, ix] = \
                    self.action.args.ccddata.data[y0:y3, ix] - scat
            self.action.args.ccddata.header[key] = (True, keycom)

        log_string = SubtractScatteredLight.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        # write out int image
        kcwi_fits_writer(self.action.args.ccddata,
                         table=self.action.args.table,
                         output_file=self.action.args.name,
                         output_dir=self.config.instrument.output_directory,
                         suffix="intd")
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix="intd")
        self.context.proctab.write_proctab()

        return self.action.args
    # END: SubtractScatteredLight()
