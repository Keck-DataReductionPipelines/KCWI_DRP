from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments
from keckdrpframework.primitives.base_img import BaseImg
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer

import os
import ccdproc


class MakeMasterFlat(BaseImg):
    """Stack flat images and make master flat image"""

    def __init__(self, action, context):
        BaseImg.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if we can create a master flat based on the processing table
        :return:
        """
        # get list of master flats
        self.logger.info("Checking precondition for MakeMasterFlat")
        self.stack_list = self.context.proctab.n_proctab(
            frame=self.action.args.ccddata, target_type='SFLAT',
            target_group=self.action.args.groupid)
        self.logger.info(f"pre condition got {len(self.stack_list)},"
                         f" expecting 1")
        # do we meet the criterion?
        if len(self.stack_list) >= 1:
            return True
        else:
            return False

    def _perform(self):
        """
        Returns an Argument() with the parameters that depends on this operation
        """
        args = self.action.args
        suffix = args.new_type.lower()

        stack_list = list(self.stack_list['OFNAME'])

        stname = stack_list[0].split('.')[0] + "_sflat.fits"

        self.logger.info("Reading image: %s" % stname)
        stacked = kcwi_fits_reader(
            os.path.join(os.path.dirname(self.action.args.name), 'redux',
                         stname))[0]

        # get master flat output name
        mfname = stack_list[0].split('.fits')[0] + '_' + suffix + '.fits'

        log_string = MakeMasterFlat.__module__ + "." + \
                     MakeMasterFlat.__qualname__
        stacked.header['HISTORY'] = log_string

        # output master flat
        kcwi_fits_writer(stacked, output_file=mfname)
        self.context.proctab.update_proctab(frame=stacked, suffix=suffix,
                                            newtype=args.new_type)
        self.context.proctab.write_proctab()

        self.logger.info(log_string)
        return self.action.args

    # END: class MakeMasterFlat()
