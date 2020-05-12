from keckdrpframework.primitives.base_img import BaseImg
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer

import os
import ccdproc


class MakeMasterDark(BaseImg):
    """Stack dark frames into master dark"""

    def __init__(self, action, context):
        BaseImg.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if we can build a stacked frame based on the processing table
        :return:
        """
        # get list of dark frames
        self.logger.info("Checking precondition for stack_darks")
        self.combine_list = self.context.proctab.n_proctab(
            frame=self.action.args.ccddata, target_type='DARK',
            target_group=self.action.args.groupid)
        self.logger.info(f"pre condition got {len(self.combine_list)},"
                         f" expecting {self.action.args.min_files}")
        # do we meet the criterion?
        if len(self.combine_list) >= self.action.args.min_files:
            return True
        else:
            return False

    def _perform(self):
        """
        Returns an Argument() with the parameters that depends on this operation
        """
        args = self.action.args
        method = 'average'
        suffix = args.new_type.lower()

        combine_list = list(self.combine_list['OFNAME'])
        # get master dark output name
        mdname = combine_list[0].split('.fits')[0] + '_' + suffix + '.fits'
        stack = []
        stackf = []
        for dark in combine_list:
            # get dark intensity (int) image file name in redux directory
            stackf.append(dark.split('.fits')[0] + '_int.fits')
            darkfn = os.path.join(args.in_directory, stackf[-1])
            # using [0] gets just the image data
            stack.append(kcwi_fits_reader(darkfn)[0])

        stacked = ccdproc.combine(stack, method=method, sigma_clip=True,
                                  sigma_clip_low_thresh=None,
                                  sigma_clip_high_thresh=2.0)
        stacked.unit = stack[0].unit
        stacked.header.IMTYPE = args.new_type
        stacked.header['NSTACK'] = (len(combine_list),
                                    'number of images stacked')
        stacked.header['STCKMETH'] = (method, 'method used for stacking')
        for ii, fname in enumerate(stackf):
            stacked.header['STACKF%d' % (ii + 1)] = (fname, "stack input file")

        log_string = MakeMasterDark.__module__
        stacked.header['HISTORY'] = log_string
        self.logger.info(log_string)

        kcwi_fits_writer(stacked, output_file=mdname,
                         output_dir=self.config.instrument.output_directory)
        self.context.proctab.update_proctab(frame=stacked, suffix=suffix,
                                            newtype=args.new_type)
        self.context.proctab.write_proctab()
        return self.action.args
    # END: class StackDarks()
