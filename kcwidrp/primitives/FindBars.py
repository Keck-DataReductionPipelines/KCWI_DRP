from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.primitives.kcwi_file_primitives import plotlabel

import numpy as np
import logging
from kcwidrp.core.kcwi_plotting import figure
from kcwidrp.core.kcwi_interactive import pause, ask_next, QUIT
from scipy.signal import find_peaks
import time


class FindBars(BasePrimitive):
    """
    Find bars in middle row of cont bars image

    Uses the middle row with surrounding rows to find the individual continuum
    bar locations in the image.  Checks to be sure the correct number of bars
    has been found (see NBARS in kcwidrp/configs/kcwi.cfg).

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
        # get grating
        grating = self.action.args.grating
        window = int(10 / y_binning)
        # get plot label
        plab = plotlabel(self.action.args)
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
        peaks_vector = None
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
            peaks_vector = []
            for pk in peaks_in_middle_vector:
                if 5 < pk < (self.action.args.ccddata.header['NAXIS1'] - 5):
                    peaks_vector.append(pk)
            peaks_found = len(peaks_vector)
            if self.config.instrument.plot_level >= 3:
                # plot the peak positions
                x = np.arange(len(middle_vector))
                p = figure(
                    title=plab + "BARS TRACE at = %d, TRY %d" % (middle_y_row,
                                                                 n_tries),
                    x_axis_label='CCD X (px)', y_axis_label='e-',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                p.line(x, middle_vector, color='blue', legend_label="Trace")
                p.scatter(peaks_vector,
                          middle_vector[peaks_vector], marker='x',
                          color='red', legend_label="FoundBar")
                p.line([0, x_size], [bar_thresh, bar_thresh],
                       color='grey', line_dash='dashed')
                p.legend.location = "bottom_center"
                bokeh_plot(p, self.context.bokeh_session)
                pause(self.config.instrument, self.logger,
                    what="Continuum-bar profile across one CCD row, with the peaks found above the threshold.",
                    help=(
                        "What you see:\n"
                        "  Title \"BARS TRACE at = <row>, TRY <n>\": the continuum-bars image\n"
                        "  across CCD x, median of a few rows around <row>.\n"
                        "  Blue: the profile. Red x: peaks found above the threshold.\n"
                        "  Grey dashed: threshold (mean + 0.5 x standard deviation).\n"
                        "What to look for:\n"
                        "  Good: exactly one red x on each of the 120 bars (24 groups of 5), none\n"
                        "  on noise or between bars.\n"
                        "  Bad: missing or extra peaks. Then the DRP moves to another row and\n"
                        "  tries again (up to 5 tries; one plot per try). If no try finds 120\n"
                        "  bars, FindBars logs an error and the geometry cannot be solved.\n"
                        "Your options:\n"
                        "  <cr>  continue. Display only: the reduction is not changed.\n"
                        "Default:\n"
                        "  <cr> continues. Below plot_level 3 this plot is not shown."),
                    who="FindBars")

            n_tries += 1
            # do we have the requisite number?
            if peaks_found != self.config.instrument.NBARS:
                self.logger.warning("Did not find %d peaks: n peaks = %d" %
                                    (self.config.instrument.NBARS, peaks_found))
                self.logger.info("Calculating new middle row and retrying")
                if camera == 0:     # BLUE
                    div_fac += 0.5
                else:               # RED
                    if 'RH4' in grating:
                        div_fac += 0.5
                    else:
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
                    title=plab + "BARS MID TRACE Thresh = %.2f" % bar_thresh,
                    x_axis_label='CCD X (px)', y_axis_label='e-',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                p.line(x, middle_vector, color='blue', legend_label="MidTrace")
                p.scatter(peaks_vector,
                          middle_vector[peaks_vector], marker='x',
                          color='red', legend_label="FoundBar")
                p.line([0, x_size], [bar_thresh, bar_thresh],
                       color='grey', line_dash='dashed')
                p.legend.location = "bottom_center"
                bokeh_plot(p, self.context.bokeh_session)
                pause(self.config.instrument, self.logger,
                    what="Continuum-bar profile at the row used for bar finding: all bars found.",
                    help=(
                        "What you see:\n"
                        "  Title \"BARS MID TRACE Thresh = <value>\": the continuum-bars profile\n"
                        "  across CCD x at the row where all bars were found.\n"
                        "  Blue: profile. Red x: the bar peaks. Grey dashed: peak threshold.\n"
                        "What to look for:\n"
                        "  Good: one red x on top of every bar, 24 groups of 5, regular spacing.\n"
                        "  Bad: a red x on a noise spike or between two bars; this shifts the\n"
                        "  bar numbering and every later step that uses it.\n"
                        "Your options:\n"
                        "  <cr>  continue. Display only: the reduction is not changed.\n"
                        "Default:\n"
                        "  <cr> continues. At plot_level 1 the plot is shown for plot_pause\n"
                        "  seconds; at plot_level 0 it is not drawn."),
                    who="FindBars")
                # calculate the bar centroids

            # Do we plot?
            if self.config.instrument.plot_level >= 3:
                do_inter = True
            else:
                do_inter = False
            for ip, peak in enumerate(peaks_vector):
                xs = list(range(peak-window, peak+window+1))
                ys = middle_vector[xs] - np.nanmin(middle_vector[xs])
                xc = np.sum(xs*ys) / np.sum(ys)
                middle_centers.append(xc)
                if do_inter:
                    p = figure(
                        title=plab + "BAR %d CENTRD = %.2f" % (ip, xc),
                        x_axis_label='CCD X (px)', y_axis_label='e-',
                        plot_width=self.config.instrument.plot_width,
                        plot_height=self.config.instrument.plot_height)
                    p.line(xs, ys, color='blue', legend_label='Bar Trace')
                    p.circle(xs, ys, color='red', legend_label='Bar Trace')
                    p.line([xc, xc], [bar_thresh, middle_vector[peak]],
                           color='green', legend_label='Cntrd')
                    bokeh_plot(p, self.context.bokeh_session)
                    if do_inter:
                        q = ask_next(self.logger,
                            what="Profile of this bar across CCD x and its centroid.",
                            help=(
                                "What you see:\n"
                                "  Title \"BAR <n> CENTRD = <x>\": bar n across CCD x, minimum subtracted.\n"
                                "  Blue line with red circles: the profile (both labelled \"Bar Trace\").\n"
                                "  Green vertical: the centroid, the x position used for this bar (its\n"
                                "  height spans threshold to peak of the raw profile, so it need not\n"
                                "  touch the blue curve).\n"
                                "What to look for:\n"
                                "  Good: one clean, roughly symmetric peak with the green line at its\n"
                                "  centre.\n"
                                "  Bad: a double or truncated peak, or a flat profile: the centroid, and\n"
                                "  the trace that starts from it, will be off.\n"
                                "Your options:\n"
                                "  <cr>      show the next bar.\n"
                                "  q / quit  stop showing centroids (all are still computed).\n"
                                "  Display only: no answer changes the reduction.\n"
                                "Default:\n"
                                "  <cr> shows the next bar. Below plot_level 3 these plots are not shown."),
                            who="FindBars")
                        if q == QUIT:
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
        try:
            if (reference_bar-1) > 0 and \
                    (reference_bar+3) < self.config.instrument.NBARS:
                for ib in range(reference_bar-1, reference_bar+3):
                    self.action.args.reference_delta_x += \
                        (middle_centers[ib] - middle_centers[ib-1])
                ndiv = 4.
            else:
                for ib in range(reference_bar, reference_bar+1):
                    self.action.args.reference_delta_x += \
                        (middle_centers[ib] - middle_centers[ib-1])
                ndiv = 2.
        except IndexError:
            self.logger.warning(
                "Not enough bars per slice to determine ref x sep")
            self.action.args.reference_delta_x = 1.
            ndiv = 1.
        self.action.args.reference_delta_x /= ndiv
        # store image info
        self.action.args.contbar_image_number = \
            self.action.args.ccddata.header['FRAMENO']
        self.action.args.contbar_image = \
            self.action.args.name
            # self.action.args.ccddata.header['OFNAME']

        log_string = FindBars.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class FindBars()
