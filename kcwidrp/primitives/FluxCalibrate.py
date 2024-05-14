from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer, \
    kcwi_fits_reader, get_master_name, strip_fname
from kcwidrp.core.kcwi_correct_extin import kcwi_correct_extin

import os
import numpy as np
from scipy.interpolate import interp1d

from astropy.io import fits as pf
from astropy import units as u
from astropy.nddata import CCDData


class FluxCalibrate(BasePrimitive):
    """
    Perform flux calibration.

    Uses inverse sensitivity curve derived from MakeInvsens to flux calibrate
    input observation.

    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if we can calibrate flux based on the processing table
        """
        self.logger.info("Checking precondition for FluxCalibrate")
        target_type = 'INVSENS'
        tab = self.context.proctab.search_proctab(
            frame=self.action.args.ccddata, target_type=target_type,
            nearest=True)
        self.logger.info("pre condition got %d invsens files, expected >= 1"
                         % len(tab))
        if len(tab) <= 0:
            self.action.args.invsname = None
            return False
        else:
            self.action.args.invsname = get_master_name(tab, target_type)
            return True

    def _perform(self):
        # Header keyword to update
        key = 'STDCOR'
        keycom = 'std corrected?'
        target_type = 'INVSENS'
        obj = None
        sky = None

        self.logger.info("Calibrating object flux")
        tab = self.context.proctab.search_proctab(
            frame=self.action.args.ccddata, target_type=target_type,
            nearest=True)
        self.logger.info("%d invsens files found" % len(tab))

        if self.action.args.invsname is not None:

            # read in master calibration (inverse sensitivity)
            invsname = self.action.args.invsname
            self.logger.info("Reading invsens: %s" % invsname)
            hdul = pf.open(os.path.join(self.config.instrument.cwd,
                                        'redux', invsname))
            mcal = hdul[0].data[1, :]
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
            if expt <= 0:
                self.logger.warning("XPOSURE at 0.0, trying TTIME")
                expt = self.action.args.ccddata.header['TTIME']
                if expt <= 0:
                    self.logger.warning("No valid exposure time found, "
                                        "using 1s")
                    expt = 1.0
            # resample onto object waves, if needed
            if w0 != mcw0 or dw != mcdw or wav[-1] != mcwav[-1] or \
                    sz[0] != mcsz[0]:
                self.logger.warning("wavelength scales not identical, "
                                    "resampling standard")
                self.logger.info(w0, mcw0, dw, mcdw, wav[-1], mcwav[-1], sz[0], mcsz[0])
                mcint = interp1d(mcwav, mcal, kind='cubic',
                                 fill_value='extrapolate')
                mscal = mcint(wav) * 1.e16 / expt
            else:
                mscal = mcal * 1.e16 / expt

            # extinction correct calibration
            kcwi_correct_extin(mscal, self.action.args.ccddata.header,
                               logger=self.logger)
            # do calibration
            for isl in range(sz[2]):
                for ix in range(sz[1]):
                    self.action.args.ccddata.data[:, ix, isl] *= mscal
                    if self.action.args.ccddata.noskysub is not None:
                        self.action.args.ccddata.noskysub[:, ix, isl] *= mscal
                    self.action.args.ccddata.uncertainty.array[:, ix, isl] *= \
                        mscal

            # check for obj, sky cubes
            if self.action.args.nasmask and self.action.args.numopen > 1:
                ofn = self.action.args.name
                # obj cube
                objfn = strip_fname(ofn) + '_ocubed.fits'
                full_path = os.path.join(
                    self.config.instrument.cwd,
                    self.config.instrument.output_directory, objfn)
                if os.path.exists(full_path):
                    obj = kcwi_fits_reader(full_path)[0]
                    # do calibration
                    for isl in range(sz[2]):
                        for ix in range(sz[1]):
                            obj.data[:, ix, isl] *= mscal
                # sky cube
                skyfn = strip_fname(ofn) + '_scubed.fits'
                full_path = os.path.join(
                    self.config.instrument.cwd,
                    self.config.instrument.output_directory, skyfn)
                if os.path.exists(full_path):
                    sky = kcwi_fits_reader(full_path)[0]
                    # do calibration
                    for isl in range(sz[2]):
                        for ix in range(sz[1]):
                            sky.data[:, ix, isl] *= mscal

            # units
            flam16_u = 1.e16 * u.erg / (u.angstrom * u.cm ** 2 * u.s)
            self.action.args.ccddata.unit = flam16_u
            self.action.args.ccddata.uncertainty.unit = flam16_u
            if obj is not None:
                obj.unit = flam16_u
            if sky is not None:
                sky.unit = flam16_u
            # update header keywords
            self.action.args.ccddata.header[key] = (True, keycom)
            self.action.args.ccddata.header['MSFILE'] = (invsname,
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
                         output_file=self.action.args.name,
                         output_dir=self.config.instrument.output_directory,
                         suffix="icubes")
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix="icubes",
                                            filename=self.action.args.name)
        self.context.proctab.write_proctab(tfil=self.config.instrument.procfile)

        # check for sky, obj cube
        if obj is not None:
            out_obj = CCDData(obj, meta=self.action.args.ccddata.header,
                              unit=self.action.args.ccddata.unit)
            kcwi_fits_writer(
                out_obj, output_file=self.action.args.name,
                output_dir=self.config.instrument.output_directory,
                suffix="ocubes")

        if sky is not None:
            out_sky = CCDData(sky, meta=self.action.args.ccddata.header,
                              unit=self.action.args.ccddata.unit)
            kcwi_fits_writer(
                out_sky, output_file=self.action.args.name,
                output_dir=self.config.instrument.output_directory,
                suffix="scubes")

        self.logger.info(log_string)

        return self.action.args
    # END: class FluxCalibrate()
