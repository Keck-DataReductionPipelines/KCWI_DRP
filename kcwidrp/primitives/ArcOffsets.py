from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments
from kcwidrp.core.bokeh_plotting import bokeh_plot

import numpy as np
from bokeh.plotting import figure, show

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
            refarc = arcs[self.config.instrument.REFBAR][:]
            # number of cross-correlation samples (avoiding ends)
            nsamp = len(refarc[10:-10])
            # possible offsets
            offar = np.arange(1-nsamp, nsamp)
            # Collect offsets
            offsets = []
            for na, arc in enumerate(arcs):
                # Cross-correlate, avoiding junk on the ends
                xcorr = np.correlate(refarc[10:-10], arc[10:-10], mode='full')
                # Calculate offset
                offset = offar[xcorr.argmax()]
                offsets.append(offset)
                self.logger.info("Arc %d Slice %d XCorr shift = %d" %
                                 (na, int(na/5), offset))
                # display if requested
                if do_plot:
                    p = figure(title=self.action.args.plotlabel +
                               "BAR OFFSET for Arc: %d Slice: %d = %d" %
                               (na, int(na/5), offset),
                               x_axis_label="CCD y (px)", y_axis_label="e-",
                               plot_width=self.config.instrument.plot_width,
                               plot_height=self.config.instrument.plot_height)
                    x = range(len(refarc))
                    p.line(x, refarc, color='green', legend='ref bar (%d)' %
                           self.config.instrument.REFBAR)
                    p.line(x, np.roll(arc, offset), color='red',
                           legend='bar %d' % na)
                    bokeh_plot(p)
                    q = input("Next? <cr>, q to quit: ")
                    if 'Q' in q.upper():
                        do_plot = False
            self.context.baroffs = offsets
        else:
            self.logger.error("No extracted arcs found")

        logstr = ArcOffsets.__module__ + "." + ArcOffsets.__qualname__
        self.action.args.ccddata.header['HISTORY'] = logstr
        self.logger.info(logstr)

        return self.action.args
    # END: class ArcOffsets()


