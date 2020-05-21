from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import parse_imsec

from astropy.nddata import StdDevUncertainty

import numpy as np


class CreateUncertaintyImage(BasePrimitive):
    """Generate a variance image based on Poisson noise plus readnoise"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        """Assumes units of image are electron"""

        # Header keyword to update
        key = 'UNCSTD'
        keycom = 'stddev uncertainty created?'

        self.logger.info("Create uncertainty image")
        # start with Poisson noise
        self.action.args.ccddata.uncertainty = StdDevUncertainty(
            np.sqrt(np.abs(self.action.args.ccddata.data)), copy=True)
        # add readnoise, if known
        if 'BIASRN1' in self.action.args.ccddata.header:
            number_of_amplifiers = self.action.args.ccddata.header['NVIDINP']
            for amplifier in range(number_of_amplifiers):
                # get amp parameters
                bias_readnoise = self.action.args.ccddata.header[
                    'BIASRN%d' % (amplifier + 1)]
                section = self.action.args.ccddata.header['ATSEC%d' %
                                                          (amplifier + 1)]
                parsed_section, read_forward = parse_imsec(section)
                self.action.args.ccddata.uncertainty.array[
                    parsed_section[0]:(parsed_section[1]+1),
                    parsed_section[2]:(parsed_section[3]+1)] = \
                    np.sqrt(
                    self.action.args.ccddata.uncertainty.array[
                        parsed_section[0]:(parsed_section[1] + 1),
                        parsed_section[2]:(parsed_section[3] + 1)] ** 2 +
                        bias_readnoise ** 2)
        else:
            self.logger.warn("Readnoise undefined, uncertainty Poisson only")
        # check for flags and mask arrays
        if self.action.args.ccddata.flags is None:
            self.action.args.ccddata.flags = np.zeros(
                self.action.args.ccddata.data.shape, dtype=np.uint8)
        if self.action.args.ccddata.mask is None:
            self.action.args.ccddata.mask = np.zeros(
                self.action.args.ccddata.data.shape, dtype=np.uint8)
        # document variance image creation
        self.action.args.ccddata.header[key] = (True, keycom)

        log_string = CreateUncertaintyImage.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class CreateUncertaintyImage()
