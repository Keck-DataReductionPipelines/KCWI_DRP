from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer

import numpy as np


class RectifyImage(BasePrimitive):
    """Ensure output image has a consistent orientation"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):

        # Header keyword to update
        key = 'IMGRECT'
        keycom = 'Image rectified?'

        # get amp mode
        ampmode = self.action.args.ccddata.header['AMPMODE'].strip().upper()

        if '__B' in ampmode or '__G' in ampmode:
            newimg = np.rot90(self.action.args.ccddata.data, 2)
            self.action.args.ccddata.data = newimg
            if self.action.args.ccddata.uncertainty:
                newunc = np.rot90(self.action.args.ccddata.uncertainty.array, 2)
                self.action.args.ccddata.uncertainty.array = newunc
            if hasattr(self.action.args.ccddata, 'mask'):
                newmask = np.rot90(self.action.args.ccddata.mask, 2)
                self.action.args.ccddata.mask = newmask
            else:
                self.logger.info("No mask data to rectify")
        elif '__D' in ampmode or '__F' in ampmode:
            newimg = np.fliplr(self.action.args.ccddata.data)
            self.action.args.ccddata.data = newimg
            if self.action.args.ccddata.uncertainty:
                newunc = np.fliplr(self.action.args.ccddata.uncertainty.array)
                self.action.args.ccddata.uncertainty.array = newunc
            if hasattr(self.action.args.ccddata, 'mask'):
                newmask = np.fliplr(self.action.args.ccddata.mask)
                self.action.args.ccddata.mask = newmask
            else:
                self.logger.info("No mask data to rectify")
        elif '__A' in ampmode or '__H' in ampmode or 'TUP' in ampmode:
            newimg = np.flipud(self.action.args.ccddata.data)
            self.action.args.ccddata.data = newimg
            if self.action.args.ccddata.uncertainty:
                newunc = np.flipud(self.action.args.ccddata.uncertainty.array)
                self.action.args.ccddata.uncertainty.array = newunc
            if hasattr(self.action.args.ccddata, 'mask'):
                newmask = np.flipud(self.action.args.ccddata.mask)
                self.action.args.ccddata.mask = newmask
            else:
                self.logger.info("No mask data to rectify")

        self.action.args.ccddata.header[key] = (True, keycom)

        logstr = RectifyImage.__module__ + "." + RectifyImage.__qualname__
        self.action.args.ccddata.header['HISTORY'] = logstr
        self.logger.info(logstr)

        # write out int image
        kcwi_fits_writer(self.action.args.ccddata,
                         table=self.action.args.table,
                         output_file=self.action.args.name, suffix="int")
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix="int")
        self.context.proctab.write_proctab()
        return self.action.args
    # END: class RectifyImage()


