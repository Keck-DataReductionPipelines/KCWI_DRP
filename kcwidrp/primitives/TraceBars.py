from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import write_table, strip_fname,\
    plotlabel
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_plotting import save_plot

from bokeh.plotting import figure
import numpy as np
import os
import time


class TraceBars(BasePrimitive):
    """
    Derive bar traces

    Using the output from FindBars.py, trace each bar to the top and bottom of
    the image.

    Uses the following configuration parameter:

        * saveintims: if set to ``True`` write out a warped version of bars image in \*_warped.fits.  Defaults to ``False``.

    Writes out resulting spatial control points as a FITS table in
    \*_trace.fits.

    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        self.logger.info("Checking for master contbars")
        if 'MCBARS' in self.action.args.ccddata.header['IMTYPE']:
            return True
        else:
            return False

    def _perform(self):
        self.logger.info("Tracing continuum bars")
        if self.config.instrument.plot_level >= 1:
            do_plot = True
        else:
            do_plot = False
        if len(self.action.args.middle_centers) < 1:
            self.logger.error("No bars found")
        elif not self.action.args.bar_avg:
            self.logger.error("No threshold for tracing")
        else:
            # initialize
            plab = plotlabel(self.action.args)
            samp = int(80 / self.action.args.ybinsize)
            win = self.action.args.window
            bar_thresh = self.action.args.bar_avg
            self.logger.info("Tracing bars with threshold of %.1f" % bar_thresh)
            xi = []     # x input
            xo = []     # x output
            yi = []     # y input (and output)
            barid = []  # bar id number
            slid = []   # slice id number
            # loop over bars
            for barn, barx in enumerate(self.action.args.middle_centers):
                # nearest pixel to bar center
                barxi = int(barx + 0.5)
                self.logger.info("bar number %d is at %.3f" % (barn, barx))
                # middle row data
                xi.append(barx)
                xo.append(barx)
                yi.append(self.action.args.middle_row)
                barid.append(barn)
                slid.append(int(barn/5))
                # trace up
                samy = self.action.args.middle_row + samp
                done = False
                while samy < (self.action.args.ccddata.data.shape[0] - win) \
                        and not done:
                    ys = np.median(
                        self.action.args.ccddata.data[(samy - win):
                                                      (samy + win + 1),
                                                      (barxi - win):
                                                      (barxi + win + 1)],
                        axis=0)
                    ys = ys - np.nanmin(ys)
                    if np.nanmax(ys) > bar_thresh and np.nansum(ys) > 0:
                        xs = list(range(barxi - win, barxi + win + 1))
                        xc = np.nansum(xs * ys) / np.nansum(ys)
                        xi.append(xc)
                        xo.append(barx)
                        yi.append(samy)
                        barid.append(barn)
                        slid.append(int(barn/5))
                        barxi = int(xc)
                    else:
                        done = True
                    samy += samp
                # trace down
                # nearest pixel to bar center
                barxi = int(barx + 0.5)
                samy = self.action.args.middle_row - samp
                done = False
                while samy >= win and not done:
                    ys = np.median(
                        self.action.args.ccddata.data[(samy - win):
                                                      (samy + win + 1),
                                                      (barxi - win):
                                                      (barxi + win + 1)],
                        axis=0)
                    ys = ys - np.nanmin(ys)
                    if np.nanmax(ys) > bar_thresh and np.nansum(ys) > 0:
                        xs = list(range(barxi - win, barxi + win + 1))
                        xc = np.sum(xs * ys) / np.sum(ys)
                        xi.append(xc)
                        xo.append(barx)
                        yi.append(samy)
                        barid.append(barn)
                        slid.append(int(barn / 5))
                        barxi = int(xc)
                    else:
                        done = True
                    samy -= samp
            # end loop over bars
            # create source and destination coords
            yo = yi
            dst = np.column_stack((xi, yi))
            src = np.column_stack((xo, yo))
            if do_plot:
                # output filename stub
                nbars = self.config.instrument.NBARS
                trcfnam = "bars_%05d_%s_%s_%s" % \
                          (self.action.args.ccddata.header['FRAMENO'],
                           self.action.args.illum, self.action.args.grating,
                           self.action.args.ifuname)
                # plot them
                p = figure(title=plab + 'SPATIAL CONTROL POINTS',
                           x_axis_label="CCD X (px)", y_axis_label="CCD Y (px)",
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.scatter(xi, yi, marker='x', size=2, color='blue')
                p.scatter(self.action.args.middle_centers,
                          [self.action.args.middle_row]*nbars, color='red')
                bokeh_plot(p, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)
                save_plot(p, filename=trcfnam+".png")
            trace = {
                'src': src,
                'dst': dst,
                'barid': barid,
                'slid': slid,
                'MIDROW': self.action.args.middle_row,
                'WINDOW': self.action.args.window,
                'REFDELX': self.action.args.reference_delta_x,
                'CBARSNO': self.action.args.contbar_image_number,
                'CBARSFL': self.action.args.contbar_image}

            # in this line we pass the trace information to an argument
            # instead of writing it to a table
            self.context.trace = trace
            name = strip_fname(self.action.args.contbar_image) + \
                "_trace.fits"
            write_table(table=[src, dst, barid, slid],
                        names=('src', 'dst', 'barid', 'slid'),
                        output_dir=os.path.join(
                            self.config.instrument.cwd,
                            self.config.instrument.output_directory),
                        output_name=name,
                        clobber=self.config.instrument.clobber,
                        comment=['Source and destination fiducial points',
                                 'Derived from KCWI continuum bars images',
                                 'For defining spatial transformation'],
                        keywords={'MIDROW': (self.action.args.middle_row,
                                             "Middle Row of image"),
                                  'WINDOW': (self.action.args.window,
                                             "Window for bar"),
                                  'REFDELX': (
                                      self.action.args.reference_delta_x,
                                      "Reference bar sep in px"),
                                  'CBARSNO': (
                                      self.action.args.contbar_image_number,
                                      "Cont. bars image number"),
                                  'CBARSFL': (self.action.args.contbar_image,
                                              "Cont. bars image")})

            if self.config.instrument.saveintims:
                from kcwidrp.primitives.kcwi_file_primitives import \
                    kcwi_fits_writer
                from skimage import transform as tf
                # fit transform
                # NOTE: we do not need an asymmetric polynomial for this
                # global fit (only for per slice fitting: see SolveGeom.py)
                self.logger.info("Fitting spatial control points")
                tform = tf.estimate_transform('polynomial', src, dst, order=3)
                self.logger.info("Transforming bars image")
                warped = tf.warp(self.action.args.ccddata.data, tform)
                # write out warped image
                self.action.args.ccddata.data = warped
                kcwi_fits_writer(
                    self.action.args.ccddata, output_file=self.action.args.name,
                    output_dir=self.config.instrument.output_directory,
                    suffix='warped')
                self.logger.info("Transformed bars produced")

            log_string = TraceBars.__module__
            self.action.args.ccddata.header['HISTORY'] = log_string
            self.logger.info(log_string)

            return self.action.args
    # END: class TraceBars()
