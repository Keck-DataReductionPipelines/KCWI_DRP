from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_plotting import get_plot_lims, oplot_slices, \
    set_plot_lims, save_plot
from kcwidrp.primitives.kcwi_file_primitives import plotlabel

import time

import numpy as np
from bokeh.plotting import figure
from scipy import signal


class ArcOffsets(BasePrimitive):
    """
    Derive offset of each bar relative to reference bar.

    Using cross correlation techniques, this routine calculates the relative
    offsets between the reference bar (loaded as a parameter as
    config.instrument.REFBAR).  The other important config parameters is
    TUKEYALPHA, which controls the roll-off of the edges of the spectra that
    are cross-correlated.  If a bright line near the edge is throwing off the
    cross-correlations, you can increase this parameter to lessen its impact.

    See kcwidrp/configs/kcwi.cfg

    Arcs must be available as context.arcs.
    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        self.logger.info("Checking for master arc")
        if 'MARC' in self.action.args.ccddata.header['IMTYPE']:
            return True
        else:
            return False

    def _perform(self):
        self.logger.info("Finding inter-bar offsets")
        nbars = self.config.instrument.NBARS
        tkalpha = self.config.instrument.TUKEYALPHA
        self.logger.info("Using Tukey Alpha value of %.2f" % tkalpha)
        arcs = self.context.arcs
        if arcs is not None:
            # Do we plot?
            do_plot = (self.config.instrument.plot_level >= 2)
            plab = plotlabel(self.action.args)
            # Compare with reference arc
            reference_arc = arcs[self.config.instrument.REFBAR][:]
            tkwgt = signal.windows.tukey(len(reference_arc), alpha=tkalpha)
            reference_arc *= tkwgt
            # number of cross-correlation samples (avoiding ends)
            number_of_samples = len(reference_arc[10:-10])
            # possible offsets
            offsets_array = np.arange(1-number_of_samples, number_of_samples)
            # Collect offsets
            offsets = []
            next_bar_to_plot = 0
            for arc_number, arc in enumerate(arcs):
                arc *= tkwgt
                # Cross-correlate, avoiding junk on the ends
                cross_correlation = np.correlate(reference_arc[10:-10],
                                                 arc[10:-10], mode='full')
                # Calculate offset
                ncross = len(cross_correlation)
                cr0 = int(ncross * 0.25)
                cr1 = int(ncross * 0.75)
                central_cross = cross_correlation[cr0:cr1]
                central_offs = offsets_array[cr0:cr1]
                offset = central_offs[central_cross.argmax()]
                offsets.append(offset)
                self.logger.info("Arc %d Slice %d XCorr shift = %d" %
                                 (arc_number, int(arc_number/5), offset))
                # display if requested
                if do_plot and arc_number == next_bar_to_plot:
                    p = figure(title=plab +
                               "BAR OFFSET for Arc: %d Slice: %d = %d" %
                               (arc_number, int(arc_number/5), offset),
                               x_axis_label="CCD y (px)", y_axis_label="e-",
                               plot_width=self.config.instrument.plot_width,
                               plot_height=self.config.instrument.plot_height)
                    x = range(len(reference_arc))
                    p.line(x, reference_arc, color='green',
                           legend_label='ref bar (%d)' %
                           self.config.instrument.REFBAR)
                    p.line(x, np.roll(arc, offset), color='red',
                           legend_label='bar %d' % arc_number)
                    bokeh_plot(p, self.context.bokeh_session)
                    q = input("Next? <int> or <cr>, q to quit: ")
                    if 'Q' in q.upper():
                        do_plot = False
                    else:
                        try:
                            next_bar_to_plot = int(q)
                        except ValueError:
                            next_bar_to_plot = arc_number + 1
            # plot output name stub
            pfname = "arc_%05d_%s_%s_%s_ta%03d" % (
                self.action.args.ccddata.header['FRAMENO'],
                self.action.args.illum, self.action.args.grating,
                self.action.args.ifuname,
                int(100 * self.config.instrument.TUKEYALPHA))
            # plot offsets
            if self.config.instrument.plot_level >= 1:
                p = figure(title=plab + "BAR OFFSETS ",
                           x_axis_label="Bar #",
                           y_axis_label="Wave offset (px)",
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.diamond(list(range(nbars)), offsets, size=8)
                xlim = [-1, nbars]
                ylim = get_plot_lims(offsets)
                p.xgrid.grid_line_color = None
                oplot_slices(p, ylim)
                set_plot_lims(p, xlim=xlim, ylim=ylim)
                bokeh_plot(p, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)
                # save offset plot
                save_plot(p, filename=pfname + '_baroffs.png')
            self.context.bar_offsets = offsets
        else:
            self.logger.error("No extracted arcs found")

        log_string = ArcOffsets.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class ArcOffsets()
