from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot

from bokeh.plotting import figure
# from bokeh.layouts import gridplot
import numpy as np
import math
import time


class SubtractOverscan(BasePrimitive):
    """Fit overscan region and subtract result from image"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        # image sections for each amp
        bsec, dsec, tsec, direc = self.action.args.map_ccd
        namps = len(bsec)
        # polynomial fit order
        if namps == 4:
            porder = 2
        else:
            porder = 7
        # header keyword to update
        key = 'OSCANSUB'
        keycom = 'Overscan subtracted?'
        # is it performed?
        performed = False
        # loop over amps
        # plts = []   # plots for each amp

        for ia in range(namps):
            # get gain
            gain = self.action.args.ccddata.header['GAIN%d' % (ia + 1)]
            # check if we have enough data to fit
            if (bsec[ia][3] - bsec[ia][2]) > self.config.instrument.minoscanpix:
                # pull out an overscan vector
                x0 = bsec[ia][2] + self.config.instrument.oscanbuf
                x1 = bsec[ia][3] - self.config.instrument.oscanbuf
                y0 = bsec[ia][0]
                y1 = bsec[ia][1] + 1
                osvec = np.nanmedian(
                    self.action.args.ccddata.data[y0:y1, x0:x1], axis=1)
                nsam = x1 - x0
                xx = np.arange(len(osvec), dtype=np.float)
                # fit it, avoiding first 50 px
                if direc[ia]:
                    # forward read skips first 50 px
                    oscoef = np.polyfit(xx[50:], osvec[50:], porder)
                else:
                    # reverse read skips last 50 px
                    oscoef = np.polyfit(xx[:-50], osvec[:-50], porder)
                # generate fitted overscan vector for full range
                osfit = np.polyval(oscoef, xx)
                # calculate residuals
                resid = (osvec - osfit) * math.sqrt(nsam) * gain / 1.414
                sdrs = float("%.3f" % np.std(resid))
                self.logger.info("Amp%d Read noise from oscan in e-: %.3f" %
                                 ((ia + 1), sdrs))
                self.action.args.ccddata.header['OSCNRN%d' % (ia + 1)] = \
                    (sdrs, "amp%d RN in e- from oscan" % (ia + 1))

                if self.config.instrument.plot_level >= 1:
                    x = np.arange(len(osvec))
                    p = figure(title=self.action.args.plotlabel +
                               'OSCAN amp %d' % (ia+1),
                               x_axis_label='overscan px',
                               y_axis_label='counts',
                               plot_width=self.config.instrument.plot_width,
                               plot_height=self.config.instrument.plot_height)
                    p.circle(x, osvec, legend_label="Data")
                    p.line(x, osfit, line_color='red', line_width=3,
                           legend_label="Fit")
                    bokeh_plot(p, self.context.bokeh_session)
                    # plts.append(p)
                    if self.config.instrument.plot_level >= 2:
                        input("Next? <cr>: ")
                    else:
                        time.sleep(self.config.instrument.plot_pause)
                # subtract it
                for ix in range(dsec[ia][2], dsec[ia][3] + 1):
                    self.action.args.ccddata.data[y0:y1, ix] = \
                        self.action.args.ccddata.data[y0:y1, ix] - osfit
                performed = True
            else:
                self.logger.info("not enough overscan px to fit amp %d")
        # if self.config.instrument.plot_level >= 1 and len(plts) > 0:
        #    bokeh_plot(gridplot(plts, ncols=(2 if namps > 2 else 1),
        #                        plot_width=500, plot_height=300,
        #                        toolbar_location=None),
        #                        self.context.bokeh_session)
        #    if self.config.instrument.plot_level >= 2:
        #        input("Next? <cr>: ")
        #    else:
        #        time.sleep(self.config.instrument.plot_pause)

        self.action.args.ccddata.header[key] = (performed, keycom)

        log_string = SubtractOverscan.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class SubtractOverscan()
