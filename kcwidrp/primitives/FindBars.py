from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments
from kcwidrp.core.bokeh_plotting import bokeh_plot, bokeh_save

import numpy as np
import logging
from bokeh.util.logconfig import basicConfig, bokeh_logger as bl
from bokeh.plotting import figure, show
from scipy.signal import find_peaks
import time


class FindBars(BasePrimitive):
    """Find bars in middle row of cont bars image"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger
        basicConfig(level=logging.ERROR)

    def _perform(self):
        self.logger.info("Finding continuum bars")
        # initialize
        refbar = self.config.instrument.REFBAR
        midcntr = []
        # get image dimensions
        nx = self.action.args.ccddata.data.shape[1]
        ny = self.action.args.ccddata.data.shape[0]
        # get binning
        ybin = self.action.args.ybinsize
        win = int(10 / ybin)
        # select from center rows of image
        midy = int(ny / 2)
        midvec = np.median(
            self.action.args.ccddata.data[(midy-win):(midy+win+1), :], axis=0)
        # set threshold for peak finding
        midavg = np.average(midvec)
        self.logger.info("peak threshold = %f" % midavg)
        # find peaks above threshold
        midpeaks, _ = find_peaks(midvec, height=midavg)
        # do we have the requisite number?
        if len(midpeaks) != self.config.instrument.NBARS:
            self.logger.error("Did not find %d peaks: n peaks = %d" %
                              (self.config.instrument.NBARS, len(midpeaks)))
        else:
            self.logger.info("found %d bars" % len(midpeaks))

            if self.config.instrument.plot_level >= 1:
                # plot the peak positions
                x = np.arange(len(midvec))
                p = figure(
                    title=self.action.args.plotlabel +
                    "BARS MID TRACE Thresh = %.2f" % midavg,
                    x_axis_label='CCD X (px)', y_axis_label='e-',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                p.line(x, midvec, color='blue', legend_label="MidTrace")
                p.scatter(midpeaks, midvec[midpeaks], marker='x', color='red',
                          legend_label="FoundBar")
                p.line([0, nx], [midavg, midavg], color='grey',
                       line_dash='dashed')
                p.legend.location = "bottom_center"
                bokeh_save(p)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)
                # calculate the bar centroids

            # Do we plot?
            if self.config.instrument.plot_level >= 2:
                do_inter = True
            else:
                do_inter = False
            for ip, peak in enumerate(midpeaks):
                xs = list(range(peak-win, peak+win+1))
                ys = midvec[xs] - np.nanmin(midvec[xs])
                xc = np.sum(xs*ys) / np.sum(ys)
                midcntr.append(xc)
                if do_inter:
                    p = figure(
                        title=self.action.args.plotlabel +
                        "BAR %d CENTRD = %.2f" % (ip, xc),
                        x_axis_label='CCD X (px)', y_axis_label='e-',
                        plot_width=self.config.instrument.plot_width,
                        plot_height=self.config.instrument.plot_height)
                    p.line(xs, ys, color='blue', legend_label='Bar Trace')
                    p.circle(xs, ys, color='red', legend_label='Bar Trace')
                    p.line([xc, xc], [midavg, midvec[peak]], color='green',
                           legend_label='Cntrd')
                    bokeh_save(p)
                    if do_inter:
                        q = input("Next? <cr>, q - quit: ")
                        if 'Q' in q.upper():
                            do_inter = False
                    else:
                        time.sleep(self.config.instrument.plot_pause)
            self.logger.info("Found middle centroids for continuum bars")
        # store peaks
        self.action.args.midcntr = midcntr
        # store the row where we got them
        self.action.args.midrow = midy
        self.action.args.win = win
        # calculate reference delta x based on refbar
        self.action.args.refdelx = 0.
        for ib in range(refbar-1, refbar+3):
            self.action.args.refdelx += (midcntr[ib] - midcntr[ib-1])
        self.action.args.refdelx /= 4.
        # store image info
        self.action.args.cbarsno = self.action.args.ccddata.header['FRAMENO']
        self.action.args.cbarsfl = self.action.args.ccddata.header['OFNAME']

        logstr = FindBars.__module__ + "." + FindBars.__qualname__
        self.action.args.ccddata.header['HISTORY'] = logstr
        self.logger.info(logstr)

        return self.action.args
    # END: class FindBars()

