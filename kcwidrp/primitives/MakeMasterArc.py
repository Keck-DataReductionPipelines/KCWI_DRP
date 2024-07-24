from keckdrpframework.primitives.base_img import BaseImg
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer, strip_fname , get_unique_STATEID_master_name

import os
import ccdproc
import numpy as np
from astropy.stats import mad_std


class MakeMasterArc(BaseImg):
    """
    Stack arc frames into master arc

    Generate a master arc frame based on the instrument config parameter
    arc_min_nframes, which defaults to 1 for the BLUE channel and 3 for the RED
    channel.  It is assumed that each frame is well-exposed and the combine
    method 'median' will be used to mitigate cosmic rays (especially for the RED
    channel).  A high sigma clipping of 2.0 is used to help with the CRs.

    Uses the ccdproc.combine routine to peform the stacking.

    Writes out a \*_marc.fits file and records a master arc frame in the proc
    table, no matter how many frames are combined.

    """

    def __init__(self, action, context):
        BaseImg.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if we can build a stacked frame based on the processing table
        """
        # get list of arc frames
        self.logger.info("Checking precondition for stack_arcs")
        self.combine_list = self.context.proctab.search_proctab(
            frame=self.action.args.ccddata, target_type='ARCLAMP',
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
        method = 'median'
        sig_up = 2.0
        suffix = args.new_type.lower()
        log_string = MakeMasterArc.__module__

        combine_list = list(self.combine_list['filename'])
        # get master arc output name
        # maname = strip_fname(combine_list[0]) + '_' + suffix + '.fits'
        if self.action.args.min_files > 1:
            stack = []
            stackf = []
            for arc in combine_list:
                # get arc intensity (int) image file name in redux directory
                stackf.append(strip_fname(arc) + '_int.fits')
                arcfn = os.path.join(self.config.instrument.cwd,
                self.config.instrument.output_directory, stackf[-1])
                # using [0] gets just the image data
                stack.append(kcwi_fits_reader(arcfn)[0])

            stacked = ccdproc.combine(stack, method=method, sigma_clip=True,
                                      sigma_clip_low_thresh=None,
                                      sigma_clip_high_thresh=sig_up,
                                      sigma_clip_func=np.ma.median,
                                      sigma_clip_dev_func=mad_std)
            stacked.unit = stack[0].unit
            stacked.header['IMTYPE'] = args.new_type
            stacked.header['NSTACK'] = (len(combine_list),
                                        'number of images stacked')
            stacked.header['STCKMETH'] = (method, 'method used for stacking')
            stacked.header['STCKSIGU'] = (sig_up,
                                          'Upper sigma rejection for stacking')
            for ii, fname in enumerate(stackf):
                stacked.header['STACKF%d' % (ii + 1)] = (fname,
                                                         "stack input file")
            stacked.header['HISTORY'] = log_string
            self.action.args.ccddata = stacked
            # maname = get_unique_STATEID_master_name(stacked, suffix)
            maname = strip_fname(combine_list[0]) + '_' + suffix + '.fits'

            kcwi_fits_writer(stacked, output_file=maname,
                             output_dir=self.config.instrument.output_directory)
            self.context.proctab.update_proctab(frame=stacked, suffix=suffix,
                                                newtype=args.new_type,
                                                filename=self.action.args.name) ### HERE
            # self.action.args.name = maname
            # self.action.args.name = stacked.header['OFNAME']
        else:
            maname = strip_fname(combine_list[0]) + '_' + suffix + '.fits'
            self.action.args.ccddata.header['IMTYPE'] = args.new_type
            self.action.args.ccddata.header['HISTORY'] = log_string
            kcwi_fits_writer(self.action.args.ccddata, output_file=maname,
                             output_dir=self.config.instrument.output_directory)
            self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                                suffix=suffix,
                                                newtype=args.new_type,
                                                filename=self.action.args.name)
        self.logger.info(log_string)

        self.context.proctab.write_proctab(tfil=self.config.instrument.procfile)
        return self.action.args
    # END: class MakeMasterArc()
