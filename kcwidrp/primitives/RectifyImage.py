from keckdrpframework.primitives.base_primitive import BasePrimitive
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
            mask = getattr(self.action.args.ccddata, "mask", None)
            if mask is not None:
                newmask = np.rot90(mask, 2)
                self.action.args.ccddata.mask = newmask
            else:
                self.logger.info("No mask data to rectify")
            flags = getattr(self.action.args.ccddata, "flags", None)
            if flags is not None:
                newflags = np.rot90(flags, 2)
                self.action.args.ccddata.flags = newflags
            else:
                self.logger.info("No flags data to rectify")
        elif '__D' in ampmode or '__F' in ampmode:
            newimg = np.fliplr(self.action.args.ccddata.data)
            self.action.args.ccddata.data = newimg
            if self.action.args.ccddata.uncertainty:
                newunc = np.fliplr(self.action.args.ccddata.uncertainty.array)
                self.action.args.ccddata.uncertainty.array = newunc
            mask = getattr(self.action.args.ccddata, "mask", None)
            if mask is not None:
                newmask = np.fliplr(mask)
                self.action.args.ccddata.mask = newmask
            else:
                self.logger.info("No mask data to rectify")
            flags = getattr(self.action.args.ccddata, "flags", None)
            if flags is not None:
                newflags = np.fliplr(flags)
                self.action.args.ccddata.flags = newflags
            else:
                self.logger.info("No flags data to rectify")
        elif '__A' in ampmode or '__H' in ampmode or 'TUP' in ampmode:
            newimg = np.flipud(self.action.args.ccddata.data)
            self.action.args.ccddata.data = newimg
            if self.action.args.ccddata.uncertainty:
                newunc = np.flipud(self.action.args.ccddata.uncertainty.array)
                self.action.args.ccddata.uncertainty.array = newunc
            mask = getattr(self.action.args.ccddata, "mask", None)
            if mask is not None:
                newmask = np.flipud(mask)
                self.action.args.ccddata.mask = newmask
            else:
                self.logger.info("No mask data to rectify")
            flags = getattr(self.action.args.ccddata, "flags", None)
            if flags is not None:
                newflags = np.flipud(flags)
                self.action.args.ccddata.flags = newflags
            else:
                self.logger.info("No flags data to rectify")

        self.action.args.ccddata.header[key] = (True, keycom)

        log_string = RectifyImage.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        # write out int image
        kcwi_fits_writer(self.action.args.ccddata,
                         table=self.action.args.table,
                         output_file=self.action.args.name,
                         output_dir=self.config.instrument.output_directory,
                         suffix="int")
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix="int")
        self.context.proctab.write_proctab()
        return self.action.args
    # END: class RectifyImage()
