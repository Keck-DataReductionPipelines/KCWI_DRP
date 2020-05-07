from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer
from kcwidrp.core.kcwi_correct_extin import kcwi_correct_extin

import os
import numpy as np
from scipy.interpolate import interp1d

from astropy.io import fits as pf


class FluxCalibrate(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if we can calibrate flux based on the processing table
        :return:
        """
        self.logger.info("Checking precondition for FluxCalibrate")
        target_type = 'INVSENS'
        tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                             target_type=target_type,
                                             nearest=True)
        self.logger.info("pre condition got %d invsens files, expected >= 1"
                         % len(tab))
        if len(tab) <= 0:
            return False
        else:
            return True

    def _perform(self):
        # Header keyword to update
        key = 'STDCOR'
        keycom = 'std corrected?'
        target_type = 'INVSENS'

        self.logger.info("Calibrating object flux")
        tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                             target_type=target_type,
                                             nearest=True)
        self.logger.info("%d invsens files found" % len(tab))

        if len(tab) > 0:

            # read in master calibration (inverse sensitivity)
            msname = tab['OFNAME'][0].split('.')[0] + '_' + \
                target_type.lower() + ".fits"
            self.logger.info("Reading invsens: %s" % msname)
            hdul = pf.open(os.path.join(os.path.dirname(self.action.args.name),
                                        'redux', msname))
            mcal = hdul[0].data
            mchdr = hdul[0].header
            hdul.close()
            # get dimensions
            mcsz = mcal.shape
            # get master std waves
            mcw0 = mchdr['CRVAL1']
            mcdw = mchdr['CDELT1']
            mcwav = mcw0 + np.arange(mcsz[0]) * mcdw
            # get master std image number
            msimgno = mchdr['FRAMENO']
            # get input object image dimensions
            sz = self.action.args.ccddata.data.shape
            # get object waves
            w0 = self.action.args.ccddata.header['CRVAL3']
            dw = self.action.args.ccddata.header['CD3_3']
            wav = w0 + np.arange(sz[0]) * dw
            # get exposure time
            expt = self.action.args.ccddata.header['XPOSURE']
            # resample onto object waves, if needed
            if w0 != mcw0 or dw != mcdw or wav[-1] != mcwav[-1] or \
                    sz[0] != mcsz[0]:
                self.logger.warning("wavelength scales not identical, "
                                    "resampling standard")
                print(w0, mcw0, dw, mcdw, wav[-1], mcwav[-1], sz[0], mcsz[0])
                mcint = interp1d(mcwav, mcal, kind='cubic',
                                 fill_value='extrapolate')
                mscal = mcint(wav) * 1.e16 / expt
            else:
                mscal = mcal * 1.e16 / expt

            # extinction correct data
            kcwi_correct_extin(self.action.args.ccddata.data,
                               self.action.args.ccddata.header,
                               logger=self.logger)
            # do calibration
            for isl in range(sz[2]):
                for ix in range(sz[1]):
                    self.action.args.ccddata.data[:, ix, isl] *= mscal
                    self.action.args.ccddata.uncertainty.array[:, ix, isl] *= \
                        mscal ** 2

            # update header keywords
            self.action.args.ccddata.header[key] = (True, keycom)
            self.action.args.ccddata.header['MSFILE'] = (msname,
                                                         "Master std filename")
            self.action.args.ccddata.header['MSIMNO'] = (
                msimgno, 'master std image number')
        else:

            self.action.args.ccddata.header[key] = (False, keycom)

        log_string = FluxCalibrate.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string

        # write out icubes image
        kcwi_fits_writer(self.action.args.ccddata,
                         table=self.action.args.table,
                         output_file=self.action.args.name, suffix="icubes")
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix="icubes")
        self.context.proctab.write_proctab()

        self.logger.info(log_string)

        return self.action.args
    # END: class FluxCalibrate()

