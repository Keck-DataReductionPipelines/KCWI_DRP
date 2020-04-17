from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments
from keckdrpframework.primitives.base_img import BaseImg
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer

import os
import numpy as np
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
            frame=self.action.args.ccddata,
            target_type=self.action.args.stack_type,
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
        suffix = self.action.args.new_type.lower()
        insuff = self.action.args.stack_type.lower()

        stack_list = list(self.stack_list['OFNAME'])

        if len(stack_list) <= 0:
            self.logger.warning("No flats found!")
            return self.action.args

        # get root for maps
        tab = self.context.proctab.n_proctab(
            frame=self.action.args.ccddata, target_type='ARCLAMP',
            target_group=self.action.args.groupid)
        if len(tab) <= 0:
            self.logger.error("Geometry not solved!")
            return self.action.args

        mroot = tab['OFNAME'][0].split('.fits')[0]

        # Wavelength map image
        wmf = mroot + '_wavemap.fits'
        self.logger.info("Reading image: %s" % wmf)
        wavemap = kcwi_fits_reader(
            os.path.join(os.path.dirname(self.action.args.name), 'redux',
                         wmf))[0]

        # Slice map image
        slf = mroot + '_slicemap.fits'
        self.logger.info("Reading image: %s" % slf)
        slicemap = kcwi_fits_reader(
            os.path.join(os.path.dirname(self.action.args.name), 'redux',
                         slf))[0]

        # Position map image
        pof = mroot + '_posmap.fits'
        self.logger.info("Reading image: %s" % pof)
        posmap = kcwi_fits_reader(
            os.path.join(os.path.dirname(self.action.args.name), 'redux',
                         pof))[0]

        # Read in stacked flat image
        stname = stack_list[0].split('.')[0] + '_' + insuff + '.fits'

        self.logger.info("Reading image: %s" % stname)
        stacked = kcwi_fits_reader(
            os.path.join(os.path.dirname(self.action.args.name), 'redux',
                         stname))[0]

        # get type of flat
        internal = ('SFLAT' in stacked.header['IMTYPE'])
        twiflat = ('STWIF' in stacked.header['IMTYPE'])
        domeflat = ('SDOME' in stacked.header['IMTYPE'])

        if internal:
            self.logger.info("Internal Flat")
        elif twiflat:
            self.logger.info("Twilight Flat")
        elif domeflat:
            self.logger.info("Dome Flat")
        else:
            self.logger.error("Flat of Unknown Type!")
            return self.action.args

        # get image size
        nx = stacked.header['NAXIS1']
        ny = stacked.header['NAXIS2']

        # get binning
        xbin = self.action.args.xbinsize
        ybin = self.action.args.ybinsize

        # Parameters for fitting

        # vignetted slice position range
        fitl = int(4/xbin)
        fitr = int(24/xbin)

        # un-vignetted slice position range
        flatl = int(34/xbin)
        flatr = int(72/xbin)

        # flat fitting slice position range
        ffleft = int(10/xbin)
        ffright = int(70/xbin)

        buffer = 5.0/float(xbin)

        # reference slice
        refslice = 9
        fflice = refslice
        ffslice2 = refslice
        sm = 25
        allidx = np.arange(int(140/xbin))

        # correct vignetting if we are using internal flats
        if internal:
            # get good region for fitting
            waves = wavemap.data.compress((wavemap.data > 0.).flat)
            waves = [waves.min(), waves.max()]
            dw = (waves[1] - waves[0]) / 30.0
            wavemin = (waves[0]+waves[1]) / 2.0 - dw
            wavemax = (waves[0]+waves[1]) / 2.0 + dw

            # get reference slice data

        # get master flat output name
        mfname = stack_list[0].split('.fits')[0] + '_' + suffix + '.fits'

        log_string = MakeMasterFlat.__module__ + "." + \
                     MakeMasterFlat.__qualname__
        stacked.header['IMTYPE'] = self.action.args.new_type
        stacked.header['HISTORY'] = log_string
        stacked.header['WAVMAPF'] = wmf
        stacked.header['SLIMAPF'] = slf
        stacked.header['POSMAPF'] = pof

        # output master flat
        kcwi_fits_writer(stacked, output_file=mfname)
        self.context.proctab.update_proctab(frame=stacked, suffix=suffix,
                                            newtype=self.action.args.new_type)
        self.context.proctab.write_proctab()

        self.logger.info(log_string)
        return self.action.args

    # END: class MakeMasterFlat()
