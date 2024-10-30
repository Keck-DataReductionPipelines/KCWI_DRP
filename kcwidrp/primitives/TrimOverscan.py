from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer

import numpy as np


class TrimOverscan(BasePrimitive):
    """
    Trim overscan region from image.

    Uses the data section (DSECn header keyword) to determine how to trim the
    image to exclude the overscan region.

    Removes raw section keywords, ASECn, BSECn, CSECn, and DSECn after trimming
    and replaces them with the ATSECn keyword giving the image section in the
    trimmed image for each amplifier.

    Uses the following configuration parameter:

        * saveintims: if set to ``True`` write out the trimmed image as \*_trim.fits.  Default is ``False``.

    If the input image is a bias frames, writes out a \*_intb.fits file,
    otherwise, just updates the image in the returned arguments.

    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):

        # parameters
        # image sections for each amp
        bsec, dsec, tsec, direc, amps, aoff = self.action.args.map_ccd
        # handle biases
        is_bias = ('BIAS' in self.action.args.ccddata.header['IMTYPE'])
        # header keyword to update
        key = 'OSCANTRM'
        keycom = 'Overscan trimmed?'
        # get output image dimensions
        max_sec = max(tsec)
        # do we have flags?
        have_flags = (self.action.args.ccddata.flags is not None)
        # create new blank image
        new = np.zeros((max_sec[1]+1, max_sec[3]+1), dtype=np.float32)
        if have_flags:
            new_flags = np.zeros((max_sec[1]+1, max_sec[3]+1), dtype=np.uint8)
        else:
            new_flags = None
        # loop over amps
        for ia in amps:
            # bias correct amp number for indexing python arrays
            iac = ia - aoff
            # input range indices
            yi0 = dsec[iac][0]
            yi1 = dsec[iac][1] + 1
            xi0 = dsec[iac][2]
            xi1 = dsec[iac][3] + 1
            # output range indices
            yo0 = tsec[iac][0]
            yo1 = tsec[iac][1] + 1
            xo0 = tsec[iac][2]
            xo1 = tsec[iac][3] + 1
            # transfer to new image
            new[yo0:yo1, xo0:xo1] = self.action.args.ccddata.data[yi0:yi1,
                                                                  xi0:xi1]
            if have_flags:
                new_flags[yo0:yo1,
                          xo0:xo1] = self.action.args.ccddata.flags[yi0:yi1,
                                                                    xi0:xi1]
            # update amp section
            sec = "[%d:" % (xo0+1)
            sec += "%d," % xo1
            sec += "%d:" % (yo0+1)
            sec += "%d]" % yo1
            self.action.args.ccddata.header[
                'ATSEC%d' % ia] = (sec, "Amp section in trimmed image")
            # remove obsolete sections
            self.action.args.ccddata.header.pop('ASEC%d' % ia)
            self.action.args.ccddata.header.pop('BSEC%d' % ia)
            self.action.args.ccddata.header.pop('DSEC%d' % ia)
            self.action.args.ccddata.header.pop('CSEC%d' % ia)
            # only in RED images
            t_key = 'TSEC%d' % ia
            if t_key in self.action.args.ccddata.header:
                self.action.args.ccddata.header.pop(t_key)
        # update with new image
        self.action.args.ccddata.data = new
        self.action.args.ccddata.flags = new_flags
        self.action.args.ccddata.header['NAXIS1'] = max_sec[3] + 1
        self.action.args.ccddata.header['NAXIS2'] = max_sec[1] + 1
        self.action.args.ccddata.header[key] = (True, keycom)

        log_string = TrimOverscan.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        if self.config.instrument.saveintims:
            kcwi_fits_writer(
                self.action.args.ccddata, table=self.action.args.table,
                output_file=self.action.args.name,
                output_dir=self.config.instrument.output_directory,
                suffix="trim")
        if is_bias:
            kcwi_fits_writer(
                self.action.args.ccddata, table=self.action.args.table,
                output_file=self.action.args.name,
                output_dir=self.config.instrument.output_directory,
                suffix="intb")
            # self.context.proctab.update_proctab(
            #     frame=self.action.args.ccddata, suffix="intb", newtype='BIAS',
            #     filename=self.action.args.ccddata.header['OFNAME'])
            self.context.proctab.update_proctab(
                frame=self.action.args.ccddata, suffix="intb", newtype='BIAS',
                filename=self.action.args.name)
            self.context.proctab.write_proctab(
                tfil=self.config.instrument.procfile)
        return self.action.args
    # END: class TrimOverscan()
