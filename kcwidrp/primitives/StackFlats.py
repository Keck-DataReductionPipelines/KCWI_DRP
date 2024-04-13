from keckdrpframework.primitives.base_img import BaseImg
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer, strip_fname, get_unique_STATEID_master_name

import os
import ccdproc
import numpy as np
from astropy.stats import mad_std


def inv_median(a):
    return 1.0/np.median(a)


class StackFlats(BaseImg):
    """
    Stack flat images

    Generate a stacked flat frame based on one of the following instrument
    config parameters:

        * flat_min_nframes: defaults to 6
        * dome_min_nframes: defaults to 3
        * twiflat_min_nframes: defaults to 1

    It is assumed that each frame is well-exposed and the combine
    method 'average' will be used.  A high sigma clipping of 2.0 is used to
    help with the CRs.

    Uses the ccdproc.combine routine to peform the stacking.

    Writes out the following files and adds entries in the proc table depending
    on the input IMTYPE:

        * FLATLAMP - writes out \*_sflat.fits and puts SFLAT in proc table
        * DOMEFLAT - writes out \*_sdome.fits and puts SDOME in proc table
        * TWIFLAT -  writes out \*_stwif.fits and puts STWIF in proc table

    """

    def __init__(self, action, context):
        BaseImg.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if we can build a stacked frame based on the processing table
        """
        # get list of input flats
        self.logger.info("Checking precondition for StackFlats")
        self.stacked_list = self.context.proctab.search_proctab(
            frame=self.action.args.ccddata,
            target_type=self.action.args.stack_type,
            target_group=self.action.args.groupid)
        if len(self.stacked_list) > 0:
            self.logger.info(f"already have {len(self.stacked_list)},"
                             f" stacked flats, expecting 0")
            return False
        else:
            self.combine_list = self.context.proctab.search_proctab(
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
        sig_up = 2.0
        suffix = self.action.args.stack_type.lower()

        self.logger.info("Stacking flats using method %s" % method)

        combine_list = list(self.combine_list['filename'])
        # get flat stack output name        
        stname = strip_fname(combine_list[0]) + '_' + suffix + '.fits'
        # stname = get_unique_STATEID_master_name(self.action.args.ccddata, suffix)
        stack = []
        stackf = []
        fmeds = []
        for flat in combine_list:
            # get flat intensity (int) image file name in redux directory
            stackf.append(strip_fname(flat) + '_intd.fits')
            flatfn = os.path.join(self.config.instrument.cwd,
                                  self.config.instrument.output_directory,
                                  stackf[0])
                                  # stackf[-1])
            # using [0] gets just the image data
            f = kcwi_fits_reader(flatfn)[0]
            # Set mask to None to prevent ccdproc.combine from masking
            f.mask = None
            # get scale
            fmed = float(np.median(f))
            self.logger.info("%s - median: %.2f" % (flat, fmed))
            fmeds.append(fmed)
            stack.append(f)

        stacked = ccdproc.combine(stack, method=method, sigma_clip=True,
                                  scale=inv_median,
                                  sigma_clip_low_thresh=None,
                                  sigma_clip_high_thresh=sig_up,
                                  sigma_clip_func=np.ma.median,
                                  sigma_clip_dev_func=mad_std)

        # re-scale stacked image
        stacked.data = stacked.data * np.median(fmeds)

        # Get the bad pixel mask out of one of the flats (is the same for all)
        # and add it to the stacked flat as the stack's mask
        last_flat_name = strip_fname(combine_list[-1]) + '_intd.fits'
        last_flat_path = os.path.join(self.config.instrument.cwd,
                                      self.config.instrument.output_directory,
                                      last_flat_name)
        last_flat = kcwi_fits_reader(last_flat_path)[0]
        stacked.mask = last_flat.mask
        
        stacked.header['IMTYPE'] = self.action.args.stack_type
        stacked.header['NSTACK'] = (len(combine_list),
                                    'number of images stacked')
        stacked.header['STCKMETH'] = (method, 'method used for stacking')
        stacked.header['STCKSIGU'] = (sig_up,
                                      'Upper sigma rejection for stacking')
        for ii, fname in enumerate(stackf):
            stacked.header['STACKF%d' % (ii + 1)] = (fname, "stack input file")

        log_string = StackFlats.__module__
        stacked.header['HISTORY'] = log_string

        # output stacked flat
        kcwi_fits_writer(stacked, output_file=stname,
                         output_dir=self.config.instrument.output_directory)

        self.context.proctab.update_proctab(frame=stacked, suffix=suffix,
                                            newtype=self.action.args.stack_type,
                                            filename=combine_list[0]) ### HERE
        self.context.proctab.write_proctab(tfil=self.config.instrument.procfile)

        self.logger.info(log_string)
        return self.action.args

    # END: class StackFlats()
