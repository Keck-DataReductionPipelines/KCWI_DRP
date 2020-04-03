from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments
from kcwidrp.primitives.kcwi_file_primitives import write_table
from kcwidrp.core.bokeh_plotting import bokeh_plot, bokeh_save

from bokeh.plotting import figure, show
from bokeh.models import Range1d, LinearAxis
import numpy as np
import os
import time


class TraceBars(BasePrimitive):
    """Derive bar traces"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Tracing continuum bars")
        if self.config.instrument.plot_level >= 1:
            do_plot = True
        else:
            do_plot = False
        if len(self.action.args.midcntr) < 1:
            self.logger.error("No bars found")
        else:
            # initialize
            samp = int(80 / self.action.args.ybinsize)
            win = self.action.args.win
            xi = []     # x input
            xo = []     # x output
            yi = []     # y input (and output)
            barid = []  # bar id number
            slid = []   # slice id number
            # loop over bars
            for barn, barx in enumerate(self.action.args.midcntr):
                # nearest pixel to bar center
                barxi = int(barx + 0.5)
                self.logger.info("bar number %d is at %.3f" % (barn, barx))
                # middle row data
                xi.append(barx)
                xo.append(barx)
                yi.append(self.action.args.midrow)
                barid.append(barn)
                slid.append(int(barn/5))
                # trace up
                samy = self.action.args.midrow + samp
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
                    xs = list(range(barxi - win, barxi + win + 1))
                    xc = np.sum(xs * ys) / np.sum(ys)
                    if np.nanmax(ys) > 255:
                        xi.append(xc)
                        xo.append(barx)
                        yi.append(samy)
                        barid.append(barn)
                        slid.append(int(barn/5))
                    else:
                        done = True
                    samy += samp
                # trace down
                samy = self.action.args.midrow - samp
                done = False
                while samy >= win and not done:
                    ys = np.median(
                        self.action.args.ccddata.data[(samy - win):
                                                      (samy + win + 1),
                                                      (barxi - win):
                                                      (barxi + win + 1)],
                        axis=0)
                    ys = ys - np.nanmin(ys)
                    xs = list(range(barxi - win, barxi + win + 1))
                    xc = np.sum(xs * ys) / np.sum(ys)
                    if np.nanmax(ys) > 255:
                        xi.append(xc)
                        xo.append(barx)
                        yi.append(samy)
                        barid.append(barn)
                        slid.append(int(barn / 5))
                    else:
                        done = True
                    # disable for now
                    samy -= samp
            # end loop over bars
            # create source and destination coords
            yo = yi
            dst = np.column_stack((xi, yi))
            src = np.column_stack((xo, yo))
            if do_plot:
                # plot them
                p = figure(title=self.action.args.plotlabel +
                           'SPATIAL CONTROL POINTS',
                           x_axis_label="CCD X (px)", y_axis_label="CCD Y (px)",
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.scatter(xi, yi, marker='x', size=2, color='blue')
                p.scatter(self.action.args.midcntr,
                          [self.action.args.midrow]*120, color='red')
                bokeh_plot(p, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)
            trace = {
                'src': src,
                'dst': dst,
                'barid': barid,
                'slid': slid,
                'MIDROW': self.action.args.midrow,
                'WINDOW': self.action.args.win,
                'REFDELX': self.action.args.refdelx,
                'CBARSNO': self.action.args.cbarsno,
                'CBARSFL': self.action.args.cbarsfl}

            # in this line we pass the trace information to an argument
            # instead of writing it to a table
            self.context.trace = trace
            ofname = self.action.args.cbarsfl.split('.')[0] + "_trace.fits"
            write_table(table=[src, dst, barid, slid],
                        names=('src', 'dst', 'barid', 'slid'),
                        output_dir=os.path.dirname(self.action.args.name),
                        output_name=ofname,
                        comment=['Source and destination fiducial points',
                                 'Derived from KCWI continuum bars images',
                                 'For defining spatial transformation'],
                        keywords={'MIDROW': (self.action.args.midrow,
                                             "Middle Row of image"),
                                  'WINDOW': (self.action.args.win,
                                             "Window for bar"),
                                  'REFDELX': (self.action.args.refdelx,
                                              "Reference bar sep in px"),
                                  'CBARSNO': (self.action.args.cbarsno,
                                              "Cont. bars image number"),
                                  'CBARSFL': (self.action.args.cbarsfl,
                                              "Cont. bars image")})

            if self.config.instrument.saveintims:
                # fit transform
                self.logger.info("Fitting spatial control points")
                tform = tf.estimate_transform('polynomial', src, dst, order=3)
                self.logger.info("Transforming bars image")
                warped = tf.warp(self.action.args.ccddata.data, tform)
                # write out warped image
                self.action.args.ccddata.data = warped
                kcwi_fits_writer(self.action.args.ccddata,
                                 self.action.args.table,
                                 output_file=self.action.args.name,
                                 suffix='warped')
                self.logger.info("Transformed bars produced")

            logstr = TraceBars.__module__ + "." + TraceBars.__qualname__
            self.action.args.ccddata.header['HISTORY'] = logstr
            self.logger.info(logstr)

            return self.action.args
    # END: class TraceBars()


