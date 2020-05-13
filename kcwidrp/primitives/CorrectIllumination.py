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
        self.action.args.master_flat = None
        self.logger.info("Checking precondition for CorrectIllumination")
        # first check for internal flat
        target_type = 'MFLAT'
        tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                             target_type=target_type,
                                             nearest=True)
        if len(tab) <= 0:
            # next look for twilight flat
            target_type = 'MTWIF'
            tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                                 target_type=target_type,
                                                 nearest=True)
            if len(tab) <= 0:
                # finally look for dome flat
                target_type = 'MDOME'
                tab = self.context.proctab.n_proctab(
                    frame=self.action.args.ccddata,
                    target_type=target_type,
                    nearest=True)
                if len(tab) <= 0:
                    precondition = False
                else:
                    precondition = True
            else:
                precondition = True
        else:
            precondition = True

        self.logger.info("pre condition got %d %s flats, expected 1"
                         % (len(tab), target_type))
        if precondition:
            self.action.args.master_flat = tab['OFNAME'][0].split('.')[0] + \
                                           '_' + target_type.lower() + ".fits"
        return precondition

    def _perform(self):

        # Header keyword to update
        key = 'FLATCOR'
        keycom = 'flat corrected?'

        self.logger.info("Correcting Illumination")
        if self.action.args.master_flat:
            mflat = kcwi_fits_reader(
                os.path.join(os.path.dirname(self.action.args.name), 'redux',
                             self.action.args.master_flat))[0]

            # do the subtraction
            self.action.args.ccddata.data *= mflat.data

            # update header keywords
            self.action.args.ccddata.header[key] = (True, keycom)
            self.action.args.ccddata.header['MFFILE'] = (
                self.action.args.master_flat, "Master flat filename")
        else:
            self.logger.error("No master flat found, "
                              "cannot correct illumination.")
            self.action.args.ccddata.header[key] = (False, keycom)

        log_string = CorrectIllumination.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string

        # write out intf image
        kcwi_fits_writer(self.action.args.ccddata,
                         table=self.action.args.table,
                         output_file=self.action.args.name,
                         output_dir=self.config.instrument.output_directory,
                         suffix="intf")
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix="intf")
        self.context.proctab.write_proctab()

        self.logger.info(log_string)

        return self.action.args
    # END: class CorrectIllumination()
