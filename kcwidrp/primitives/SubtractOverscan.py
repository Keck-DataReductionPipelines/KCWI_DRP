from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_plotting import save_plot
from kcwidrp.core.bokeh_plotting import bokeh_clear

from bokeh.plotting import figure
from bokeh.layouts import gridplot
import numpy as np
import math
import time


class SubtractOverscan(BasePrimitive):
    """
    Determines overscan offset and subtracts it from the image.

    Uses the BIASSEC header keyword to determine where to calculate the overscan
    offset.  Subtracts the overscan offset and records the value in the header.
    In addition, performs a polynomial fit and uses the residuals to determine
    the read noise in the overscan.  Records the overscan readnoise in the
    header as well.

    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        # image sections for each amp
        bsec, dsec, tsec, direc, amps, aoff = self.action.args.map_ccd
        namps = len(amps)
        # polynomial fit order
        if namps == 4:
            porder = 2
        else:
            porder = 7
        minoscanpix = self.config.instrument.minoscanpix
        oscanbuf = self.config.instrument.oscanbuf
        frameno = self.action.args.ccddata.header['FRAMENO']
        # header keyword to update
        key = 'OSCANSUB'
        keycom = 'Overscan subtracted?'
        # is it performed?
        performed = False
        # loop over amps
        plts = []   # plots for each amp

        for ia in amps:
            # bias correct amp number for indexing python arrays
            iac = ia - aoff
            # get gain
            gain = self.action.args.ccddata.header['GAIN%d' % ia]
            # check if we have enough data to fit
            if (bsec[iac][3] - bsec[iac][2]) > minoscanpix:
                # pull out an overscan vector
                x0 = bsec[iac][2] + oscanbuf
                x1 = bsec[iac][3] - oscanbuf
                y0 = bsec[iac][0]
                y1 = bsec[iac][1] + 1
                # get overscan value to subtract
                osval = int(np.nanmedian(
                    self.action.args.ccddata.data[y0:y1, x0:x1]))
                # vector to fit for determining read noise
                osvec = np.nanmedian(
                    self.action.args.ccddata.data[y0:y1, x0:x1], axis=1)
                nsam = x1 - x0
                xx = np.arange(len(osvec), dtype=np.float32)
                # fit it, avoiding first 50 px
                if direc[iac][0]:
                    # forward read skips first 50 px
                    oscoef = np.polyfit(xx[50:], osvec[50:], porder)
                    # generate fitted overscan vector for full range
                    osfit = np.polyval(oscoef, xx)
                    # calculate residuals
                    resid = (osvec[50:] - osfit[50:]) * math.sqrt(nsam) * \
                        gain / 1.414
                else:
                    # reverse read skips last 50 px
                    oscoef = np.polyfit(xx[:-50], osvec[:-50], porder)
                    # generate fitted overscan vector for full range
                    osfit = np.polyval(oscoef, xx)
                    # calculate residuals
                    resid = (osvec[:-50] - osfit[:-50]) * math.sqrt(nsam) * \
                        gain / 1.414

                sdrs = float("%.3f" % np.std(resid))
                self.logger.info("Img # %05d, Amp %d [%d:%d, %d:%d]" % (frameno,
                                                                        ia,
                                                                        x0, x1,
                                                                        y0, y1))
                self.logger.info("Amp%d oscan counts (DN): %d" %
                                 (ia, osval))
                self.logger.info("Amp%d Read noise from oscan in e-: %.3f" %
                                 (ia, sdrs))
                self.action.args.ccddata.header['OSCNRN%d' % ia] = \
                    (sdrs, "amp%d RN in e- from oscan" % ia)
                self.action.args.ccddata.header['OSCNVAL%d' % ia] = \
                    (osval, "amp%d oscan counts (DN)" % ia)

                if self.config.instrument.plot_level >= 2:
                    x = np.arange(len(osvec))
                    p = figure(title='Img # %05d OSCAN [%d:%d, %d:%d] '
                                     'amp %d, noise: %.3f e-/px' %
                                     (frameno, x0, x1, y0, y1, ia, sdrs),
                               x_axis_label='overscan px',
                               y_axis_label='counts',
                               plot_width=self.config.instrument.plot_width,
                               plot_height=self.config.instrument.plot_height)
                    p.circle(x, osvec, legend_label="Data")
                    p.line(x, osfit, line_color='red', line_width=3,
                           legend_label="Fit")
                    bokeh_plot(p, self.context.bokeh_session)
                    plts.append(p)
                    if self.config.instrument.plot_level >= 3:
                        input("Next? <cr>: ")
                    else:
                        time.sleep(self.config.instrument.plot_pause)
                # subtract osval
                xx0 = dsec[iac][2]
                xx1 = dsec[iac][3] + 1
                self.action.args.ccddata.data[y0:y1, xx0:xx1] -= osval
                # for ix in range(dsec[iac][2], dsec[iac][3] + 1):
                #     self.action.args.ccddata.data[y0:y1, ix] = \
                #         self.action.args.ccddata.data[y0:y1, ix] - osfit
                performed = True
            else:
                self.logger.info("not enough overscan px to fit amp %d" % ia)
                performed = False
        if self.config.instrument.plot_level >= 3 and len(plts) > 0:
            bokeh_plot(gridplot(plts, ncols=(2 if namps > 2 else 1),
                                plot_width=500, plot_height=300,
                                toolbar_location=None),
                       self.context.bokeh_session)
            save_plot(gridplot(plts, ncols=(2 if namps > 2 else 1),
                               plot_width=500, plot_height=300,
                               toolbar_location=None),
                      filename="oscan_%05d.png" % frameno)
            if self.config.instrument.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.config.instrument.plot_pause)
            bokeh_clear(self.context.bokeh_session)

        self.action.args.ccddata.header[key] = (performed, keycom)

        log_string = SubtractOverscan.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class SubtractOverscan()
