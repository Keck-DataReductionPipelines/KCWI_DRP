from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot

import numpy as np
from bokeh.plotting import figure


class ArcOffsets(BasePrimitive):
    """Derive offset of each bar relative to reference bar"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Finding inter-bar offsets")
        arcs = self.context.arcs
        if arcs is not None:
            # Do we plot?
            if self.config.instrument.plot_level >= 2:
                do_plot = True
            else:
                do_plot = False
            # Compare with reference arc
            reference_arc = arcs[self.config.instrument.REFBAR][:]
            # number of cross-correlation samples (avoiding ends)
            number_of_samples = len(reference_arc[10:-10])
            # possible offsets
            offsets_array = np.arange(1-number_of_samples, number_of_samples)
            # Collect offsets
            offsets = []
            for arc_number, arc in enumerate(arcs):
                # Cross-correlate, avoiding junk on the ends
                cross_correlation = np.correlate(reference_arc[10:-10],
                                                 arc[10:-10], mode='full')
                # Calculate offset
                offset = offsets_array[cross_correlation.argmax()]
                offsets.append(offset)
                self.logger.info("Arc %d Slice %d XCorr shift = %d" %
                                 (arc_number, int(arc_number/5), offset))
                # display if requested
                if do_plot:
                    p = figure(title=self.action.args.plotlabel +
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
                    q = input("Next? <cr>, q to quit: ")
                    if 'Q' in q.upper():
                        do_plot = False
            self.context.bar_offsets = offsets
        else:
            self.logger.error("No extracted arcs found")

        log_string = ArcOffsets.__module__ + "." + ArcOffsets.__qualname__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class ArcOffsets()
