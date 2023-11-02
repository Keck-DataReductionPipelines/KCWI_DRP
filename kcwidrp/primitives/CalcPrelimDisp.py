from keckdrpframework.primitives.base_primitive import BasePrimitive

import math


class CalcPrelimDisp(BasePrimitive):
    """
    Calculate dispersion based on configuration parameters.

    The parameters of the grating equation are calculated as:

    * BLUE

      * alpha = grating_angle - 13 - adjustment_angle
      * adjustment_angle = (180 for BH and 0 for all other gratings)

    * RED

      * alpha = 155.892 - grating_angle

    * BLUE & RED

      * beta = camera_angle - alpha

    dispersion = cos(beta)/rho/focal_length x (pixel_scale x binning) * 1.e4

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
        # get binning
        y_binning = self.action.args.ybinsize
        # 0 - compute alpha
        if 'BLUE' in self.action.args.ccddata.header['CAMERA'].upper():
            preliminary_alpha = self.action.args.grangle - 13.0 - \
                self.action.args.adjang
        else:
            # red_grat_norm_angle = 156.1748047  # Caltech AIT value
            red_grat_norm_angle = 155.892
            preliminary_alpha = red_grat_norm_angle - self.action.args.grangle
        # 1 - compute preliminary angle of diffraction
        preliminary_beta = self.action.args.camangle - preliminary_alpha
        # 2 - compute preliminary dispersion
        preliminary_dispersion = math.cos(preliminary_beta/math.degrees(1.)) / \
            self.action.args.rho / self.config.instrument.FCAM * \
            (self.config.instrument.PIX*y_binning) * 1.e4
        preliminary_dispersion *= math.cos(
            self.config.instrument.GAMMA/math.degrees(1.))
        self.logger.info("Initial alpha, beta (deg): %.3f, %.3f" %
                         (preliminary_alpha, preliminary_beta))
        self.logger.info("Initial calculated dispersion (A/binned pix): %.3f" %
                         preliminary_dispersion)
        self.context.prelim_disp = preliminary_dispersion

        log_string = CalcPrelimDisp.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class CalcPrelimDisp()
