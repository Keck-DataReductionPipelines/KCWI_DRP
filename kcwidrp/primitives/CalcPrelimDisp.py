from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments

import math


class CalcPrelimDisp(BasePrimitive):
    """Calculate dispersion based on configuration parameters"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        # get binning
        ybin = self.action.args.ybinsize
        # 0 - compute alpha
        prelim_alpha = self.action.args.grangle - 13.0 - self.action.args.adjang
        # 1 - compute preliminary angle of diffraction
        prelim_beta = self.action.args.camangle - prelim_alpha
        # 2 - compute preliminary dispersion
        prelim_disp = math.cos(prelim_beta/math.degrees(1.)) / \
            self.action.args.rho / self.config.instrument.FCAM * \
            (self.config.instrument.PIX*ybin) * 1.e4
        prelim_disp *= math.cos(
            self.config.instrument.GAMMA/math.degrees(1.))
        self.logger.info("Initial alpha, beta (deg): %.3f, %.3f" %
                         (prelim_alpha, prelim_beta))
        self.logger.info("Initial calculated dispersion (A/binned pix): %.3f" %
                         prelim_disp)
        self.context.prelim_disp = prelim_disp

        logstr = CalcPrelimDisp.__module__ + "." + CalcPrelimDisp.__qualname__
        self.action.args.ccddata.header['HISTORY'] = logstr
        self.logger.info(logstr)

        return self.action.args
    # END: class CalcPrelimDisp()


