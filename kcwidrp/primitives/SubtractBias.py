from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader

import os


class SubtractBias(BasePrimitive):
    """Subtract master bias frame"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):

        # Header keyword to update
        key = 'BIASSUB'
        keycom = 'master bias subtracted?'
        target_type = 'MBIAS'

        self.logger.info("Subtracting master bias")
        tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                             target_type=target_type,
                                             nearest=True)
        self.logger.info("%d master bias frames found" % len(tab))

        if len(tab) > 0:
            mbname = tab['OFNAME'][0].split('.')[0] + '_' + \
                     target_type.lower() + ".fits"
            self.logger.info("Reading image: %s" % mbname)
            mbias = kcwi_fits_reader(
                os.path.join(os.path.dirname(self.action.args.name), 'redux',
                             mbname))[0]

            # do the subtraction
            self.action.args.ccddata.data -= mbias.data

            # transfer bias read noise
            namps = self.action.args.ccddata.header['NVIDINP']
            for ia in range(namps):
                self.action.args.ccddata.header['BIASRN%d' % (ia + 1)] = \
                    mbias.header['BIASRN%d' % (ia + 1)]

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
