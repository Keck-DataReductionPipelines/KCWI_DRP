from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer

import numpy as np


class FlagSaturation(BasePrimitive):
    """Remove known bad columns"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Flagging saturated pixels")

        # Header keyword to update
        key = 'SATFLAG'
        keycom = 'flagged saturated pixels?'

        # Create flags for saturated pixels
        flags = np.zeros(self.action.args.ccddata.data.shape, dtype=np.uint8)

        sat = np.where(self.action.args.ccddata.data > 60000)
        flags[sat] += 8
        number_of_sat_pixels = np.sum(sat)

        sat = np.where(self.action.args.ccddata.data == 0)
        flags[sat] += 8
        number_of_sat_pixels += np.sum(sat)

        self.action.args.ccddata.header[key] = (True, keycom)

        self.logger.info("Flagged %d saturated pixels" % number_of_sat_pixels)
        self.action.args.ccddata.header['NSATFLAG'] = \
            (number_of_sat_pixels, 'number of saturated pixels flagged')

        log_string = FlagSaturation.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        # add flags array
        # DN 2023-may-28: commenting out mask update because
        # it causes bad things later on
        # self.action.args.ccddata.mask = flags
        self.action.args.ccddata.flags = flags

        if self.config.instrument.saveintims:
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name,
                             output_dir=self.config.instrument.output_directory,
                             suffix="fsat")

        return self.action.args
    # END: class CorrectDefects()
