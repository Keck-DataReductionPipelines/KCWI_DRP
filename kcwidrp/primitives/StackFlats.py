from keckdrpframework.primitives.base_img import BaseImg
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer

import os
import ccdproc


class StackFlats(BaseImg):
    """Stack flat images"""

    def __init__(self, action, context):
        BaseImg.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if we can build a stacked frame based on the processing table
        :return:
        """
        # get list of input flats
        self.logger.info("Checking precondition for StackFlats")
        self.stacked_list = self.context.proctab.n_proctab(
            frame=self.action.args.ccddata,
            target_type=self.action.args.stack_type,
            target_group=self.action.args.groupid)
        if len(self.stacked_list) > 0:
            self.logger.info(f"already have {len(self.stacked_list)},"
                             f" stacked flats, expecting 0")
            return False
        else:
            self.combine_list = self.context.proctab.n_proctab(
                frame=self.action.args.ccddata,
                target_type=self.action.args.want_type,
                target_group=self.action.args.groupid)
            self.logger.info(f"pre condition got {len(self.combine_list)},"
                             f" expecting {self.action.args.min_files}")
            # do we meet the criterion?
            if len(self.combine_list) >= self.action.args.min_files:
                return True
            else:
                self.action.new_event = None
                return False

    def _perform(self):
        """
        Returns an Argument() with the parameters that depends on this operation
        """
        method = 'average'
        suffix = self.action.args.stack_type.lower()

        self.logger.info("Stacking flats using method %s" % method)

        combine_list = list(self.combine_list['OFNAME'])
        # get flat stack output name
        stname = combine_list[0].split('.fits')[0] + '_' + suffix + '.fits'
        stack = []
        stackf = []
        for flat in combine_list:
            # get flat intensity (int) image file name in redux directory
            stackf.append(flat.split('.fits')[0] + '_intd.fits')
            flatfn = os.path.join(self.action.args.in_directory, stackf[-1])
            # using [0] gets just the image data
            stack.append(kcwi_fits_reader(flatfn)[0])

        stacked = ccdproc.combine(stack, method=method, sigma_clip=True,
                                  sigma_clip_low_thresh=None,
                                  sigma_clip_high_thresh=2.0)
        stacked.header['IMTYPE'] = self.action.args.stack_type
        stacked.header['NSTACK'] = (len(combine_list),
                                    'number of images stacked')
        stacked.header['STCKMETH'] = (method, 'method used for stacking')
        for ii, fname in enumerate(stackf):
            stacked.header['STACKF%d' % (ii + 1)] = (fname, "stack input file")

        log_string = StackFlats.__module__
        stacked.header['HISTORY'] = log_string

        # output stacked flat
        kcwi_fits_writer(stacked, output_file=stname,
                         output_dir=self.config.instrument.output_directory)

        self.context.proctab.update_proctab(frame=stacked, suffix=suffix,
                                            newtype=self.action.args.stack_type)
        self.context.proctab.write_proctab()

        self.logger.info(log_string)
        return self.action.args

    # END: class StackFlats()
