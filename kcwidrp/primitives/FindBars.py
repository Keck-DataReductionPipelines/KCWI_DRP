from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot

import numpy as np
import logging
from bokeh.util.logconfig import basicConfig
from bokeh.plotting import figure
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
        reference_bar = self.config.instrument.REFBAR
        middle_centers = []
        # get image dimensions
        x_size = self.action.args.ccddata.data.shape[1]
        y_size = self.action.args.ccddata.data.shape[0]
        # get binning
        y_binning = self.action.args.ybinsize
        # get camera
        camera = self.action.args.camera
        window = int(10 / y_binning)
        # select from center rows of image
        div_fac = 2
        middle_y_row = int(y_size / div_fac)
        # handle dichroic edge by moving middle until we get enough peaks
        n_tries = 0
        peaks_found = 0
        middle_vector = None
        bar_thresh = None
        average_value_middle_vector = None
        stdev_value_middle_vector = None
        peaks_in_middle_vector = None
        while peaks_found != self.config.instrument.NBARS and n_tries < 5:
            self.logger.info("Middle row for bars finding: %d" % middle_y_row)
            middle_vector = np.median(
                self.action.args.ccddata.data[
                    (middle_y_row-window):(middle_y_row+window+1), :], axis=0)
            # set threshold for peak finding
            average_value_middle_vector = np.average(middle_vector)
            stdev_value_middle_vector = np.nanstd(middle_vector)
            bar_thresh = average_value_middle_vector + \
                0.5 * stdev_value_middle_vector
            self.logger.info("peak threshold = %f" % bar_thresh)
            # find peaks above threshold
            peaks_in_middle_vector, _ = find_peaks(
                middle_vector, height=bar_thresh)
            peaks_found = len(peaks_in_middle_vector)
            n_tries += 1
            # do we have the requisite number?
            if peaks_found != self.config.instrument.NBARS:
                self.logger.warning("Did not find %d peaks: n peaks = %d" %
                                    (self.config.instrument.NBARS, peaks_found))
                self.logger.info("Calculating new middle row and retrying")
                if camera == 0:     # BLUE
                    div_fac += 0.5
                else:               # RED
                    div_fac -= 0.5
                middle_y_row = int(y_size / div_fac)
        if peaks_found != self.config.instrument.NBARS:
            self.logger.error("Did not find %d peaks: n peaks = %d" %
                              (self.config.instrument.NBARS, peaks_found))
        else:
            self.logger.info("found %d bars" % peaks_found)

            if self.config.instrument.plot_level >= 1:
                # plot the peak positions
                x = np.arange(len(middle_vector))
                p = figure(
                    title=self.action.args.plotlabel +
                    "BARS MID TRACE Thresh = %.2f" % bar_thresh,
                    x_axis_label='CCD X (px)', y_axis_label='e-',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                p.line(x, middle_vector, color='blue', legend_label="MidTrace")
                p.scatter(peaks_in_middle_vector,
                          middle_vector[peaks_in_middle_vector], marker='x',
                          color='red', legend_label="FoundBar")
                p.line([0, x_size], [bar_thresh, bar_thresh],
                       color='grey', line_dash='dashed')
                p.legend.location = "bottom_center"
                bokeh_plot(p, self.context.bokeh_session)
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
            for ip, peak in enumerate(peaks_in_middle_vector):
                xs = list(range(peak-window, peak+window+1))
                ys = middle_vector[xs] - np.nanmin(middle_vector[xs])
                xc = np.sum(xs*ys) / np.sum(ys)
                middle_centers.append(xc)
                if do_inter:
                    p = figure(
                        title=self.action.args.plotlabel +
                        "BAR %d CENTRD = %.2f" % (ip, xc),
                        x_axis_label='CCD X (px)', y_axis_label='e-',
                        plot_width=self.config.instrument.plot_width,
                        plot_height=self.config.instrument.plot_height)
                    p.line(xs, ys, color='blue', legend_label='Bar Trace')
                    p.circle(xs, ys, color='red', legend_label='Bar Trace')
                    p.line([xc, xc], [bar_thresh, middle_vector[peak]],
                           color='green', legend_label='Cntrd')
                    bokeh_plot(p, self.context.bokeh_session)
                    if do_inter:
                        q = input("Next? <cr>, q - quit: ")
                        if 'Q' in q.upper():
                            do_inter = False
                    else:
                        time.sleep(self.config.instrument.plot_pause)
            self.logger.info("Found middle centroids for continuum bars")
        # store threshold
        self.action.args.bar_thresh = bar_thresh
        self.action.args.bar_avg = average_value_middle_vector
        self.action.args.bar_std = stdev_value_middle_vector
        # store peaks
        self.action.args.middle_centers = middle_centers
        # store the row where we got them
        self.action.args.middle_row = middle_y_row
        self.action.args.window = window
        # calculate reference delta x based on refbar
        self.action.args.reference_delta_x = 0.
        for ib in range(reference_bar-1, reference_bar+3):
            self.action.args.reference_delta_x += \
                (middle_centers[ib] - middle_centers[ib-1])
        self.action.args.reference_delta_x /= 4.
        # store image info
        self.action.args.contbar_image_number = \
            self.action.args.ccddata.header['FRAMENO']
        self.action.args.contbar_image = \
            self.action.args.ccddata.header['OFNAME']

        log_string = FindBars.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class FindBars()
