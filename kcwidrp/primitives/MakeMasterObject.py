from keckdrpframework.primitives.base_img import BaseImg
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer, strip_fname, get_master_name

import os
import ccdproc
import numpy as np
from astropy.stats import mad_std


class MakeMasterObject(BaseImg):
    """
    Stack object frames into master object

    Generate a master object frame based on the instrument config parameter
    object_min_nframes, which defaults to 1 for the BLUE channel and 1 for the
    RED channel.  The combine method 'median' will be used to mitigate cosmic
    rays (especially for the RED channel) if object_min_nframes is 3 or less.
    If larger than 3, then the method 'average' will be used.  A high sigma
    clipping of 2.0 is used to help with the CRs.

    Uses the ccdproc.combine routine to peform the stacking.

    Writes out a \*_mobj.fits file and records a master object frame in the proc
    table, no matter how many frames are combined.
    """

    def __init__(self, action, context):
        BaseImg.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if we can build a stacked frame based on the processing table
        :return:
        """
        # get list of arc frames
        self.logger.info("Checking precondition for MakeMasterObject")
        self.combine_list = self.context.proctab.search_proctab(
            frame=self.action.args.ccddata, target_type='OBJECT',
            target_group=self.action.args.groupid)
        self.logger.info(f"pre condition got {len(self.combine_list)},"
                         f" expecting {self.action.args.min_files}")

        # do we meet the criterion?
        if len(self.combine_list) >= self.action.args.min_files > 1:
            return True
        else:
            if self.action.args.min_files <= 1:
                self.logger.info(f"skipping stacking because only one combined")
            else:
                self.action.new_event = None
            return False

    def _perform(self):
        """
        Returns an Argument() with the parameters that depends on this operation
        """
        args = self.action.args
        method = 'median'   # default for 3 or fewer images
        sig_up = 2.0        # default upper sigma rejection limit
        suffix = args.new_type.lower()
        log_string = MakeMasterObject.__module__

        combine_list = list(self.combine_list['filename'])

        # more than three and we can average
        nstack = len(combine_list)
        if nstack > 3:
            method = 'average'

        self.logger.info("Combining Master Object with method %s" % method)

        # get master arc output name
        # maname = strip_fname(combine_list[0]) + '_' + suffix + '.fits'
        if self.action.args.min_files > 1:
            stack = []
            stacko = []
            for obj in combine_list:
                # get object intensity (int) image file name in redux directory
                stackf = strip_fname(obj) + '_intf.fits'
                objfn = os.path.join(self.config.instrument.cwd,
                self.config.instrument.output_directory, stackf)
                stacko.append(stackf)
                # using [0] gets just the image data
                stack.append(kcwi_fits_reader(objfn)[0])

            stacked = ccdproc.combine(stack, method=method, sigma_clip=True,
                                      sigma_clip_low_thresh=None,
                                      sigma_clip_high_thresh=sig_up,
                                      sigma_clip_func=np.ma.median,
                                      sigma_clip_dev_func=mad_std)
            stacked.unit = stack[0].unit
            stacked.header['IMTYPE'] = args.new_type
            stacked.header['NSTACK'] = (nstack, 'number of images stacked')
            stacked.header['STCKMETH'] = (method, 'method used for stacking')
            stacked.header['STCKSIGU'] = (sig_up,
                                          'Upper sigma rejection for stacking')
            for ii, fname in enumerate(stacko):
                stacked.header['STACKF%d' % (ii + 1)] = (fname,
                                                         "stack input file")
            stacked.header['HISTORY'] = log_string
            self.action.args.ccddata = stacked
            mobj_name = get_master_name(combine_list, "mobj")
            kcwi_fits_writer(stacked, output_file=mobj_name,
                             output_dir=self.config.instrument.output_directory)
            self.context.proctab.update_proctab(frame=stacked, suffix=suffix,
                                                newtype=args.new_type,
                                                filename=self.action.args.name) ### HERE
            # self.action.args.name = mobj_name
            # self.action.args.name = stacked.header['OFNAME']
        else:
            mobj_name = get_master_name(combine_list, "mobj")
            self.action.args.ccddata.header['IMTYPE'] = args.new_type
            self.action.args.ccddata.header['HISTORY'] = log_string
            kcwi_fits_writer(self.action.args.ccddata, output_file=mobj_name,
                             output_dir=self.config.instrument.output_directory)
            self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                                suffix=suffix,
                                                newtype=args.new_type,
                                                filename=self.action.args.name)
        self.logger.info(log_string)

        self.context.proctab.write_proctab(tfil=self.config.instrument.procfile)
        return self.action.args
    # END: class MakeMasterObject()
