from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer

import os


class CorrectIllumination(BasePrimitive):
    """Subtract master bias frame"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if we can correct illumination based on the processing table
        :return:
        """
        self.logger.info("Checking precondition for CorrectIllumination")
        target_type = 'MFLAT'
        tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                             target_type=target_type,
                                             nearest=True)
        self.logger.info("pre condition got %d master flats, expected 1"
                         % len(tab))
        if len(tab) <= 0:
            return False
        else:
            return True

    def _perform(self):

        # Header keyword to update
        key = 'FLATCOR'
        keycom = 'flat corrected?'
        target_type = 'MFLAT'

        self.logger.info("Correcting Illumination")
        tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                             target_type=target_type,
                                             nearest=True)
        self.logger.info("%d master flat frames found" % len(tab))

        if len(tab) > 0:
            mfname = tab['OFNAME'][0].split('.')[0] + '_' + \
                     target_type.lower() + ".fits"
            self.logger.info("Reading image: %s" % mfname)
            mflat = kcwi_fits_reader(
                os.path.join(os.path.dirname(self.action.args.name), 'redux',
                             mfname))[0]

            # do the subtraction
            self.action.args.ccddata.data *= mflat.data

            # update header keywords
            self.action.args.ccddata.header[key] = (True, keycom)
            self.action.args.ccddata.header['MFFILE'] = (mfname,
                                                         "Master flat filename")
        else:

            self.action.args.ccddata.header[key] = (False, keycom)

        log_string = CorrectIllumination.__module__ + "." + \
            CorrectIllumination.__qualname__
        self.action.args.ccddata.header['HISTORY'] = log_string

        # write out int image
        kcwi_fits_writer(self.action.args.ccddata,
                         table=self.action.args.table,
                         output_file=self.action.args.name, suffix="intf")
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix="intf")
        self.context.proctab.write_proctab()

        self.logger.info(log_string)

        return self.action.args
    # END: class CorrectIllumination()
