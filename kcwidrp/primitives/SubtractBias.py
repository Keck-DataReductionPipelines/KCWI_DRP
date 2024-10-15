from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    get_master_name, get_unique_CCD_master_name

import os


class SubtractBias(BasePrimitive):
    """
    Subtract the master bias frame.

    Reads in the master bias created by MakeMasterBias.py and performs the
    subtraction (after verifying amplifier configuration agreement).  Records
    the processing in the header.

    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):

        # Header keyword to update
        key = 'BIASSUB'
        keycom = 'master bias subtracted?'
        target_type = 'MBIAS'

        self.logger.info("Subtracting master bias")
        tab = self.context.proctab.search_proctab(
            frame=self.action.args.ccddata, target_type=target_type,
            nearest=True)
        self.logger.info("%d master bias frames found" % len(tab))

        if len(tab) > 0:
            # mbname = get_master_name(tab, target_type)
            mbname = get_unique_CCD_master_name(self.action.args.ccddata)

            # mbname = master_bias_name(self.action.args.ccddata)
            self.logger.info("Reading image: %s" % mbname)
            mbias = kcwi_fits_reader(
                os.path.join(self.context.config.instrument.cwd, 'redux',
                             mbname))[0]

            # do the subtraction
            self.action.args.ccddata.data -= mbias.data
            bsec, dsec, tsec, direc, amps, aoff = self.action.args.map_ccd

            # transfer bias read noise
            namps = self.action.args.ccddata.header['NVIDINP']
            if len(amps) != namps:
                self.logger.warning("Amp count disagreement!")
            for ia in amps:
                self.action.args.ccddata.header['BIASRN%d' % ia] = \
                    mbias.header['BIASRN%d' % ia]

            self.action.args.ccddata.header[key] = (True, keycom)
            self.action.args.ccddata.header['MBFILE'] = (mbname,
                                                         "Master bias filename")
        else:

            self.action.args.ccddata.header[key] = (False, keycom)

        log_string = SubtractBias.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class SubtractBias()
