from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.kcwi_get_std import kcwi_get_std
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer, \
    kcwi_fits_reader

import numpy as np
from scipy.ndimage import shift
from astropy.nddata import CCDData
import math
import ref_index
import os


def atm_disper(w0, w1, airmass, temperature=10.0, pressure_pa=61100.0,
               humidity=50.0, co2=400.0):
    """Calculate atmospheric dispersion at w1 relative to w0

    Args:
        w0 (float): reference wavelength (Angstroms)
        w1 (float): offset wavelength (Angstroms)
        airmass (float): unitless airmass
        temperature (float): atmospheric temperature (C)
        pressure_pa (float): atmospheric pressure (Pa)
        humidity (float): relative humidity (%)
        co2 (float): Carbon-Dioxide (mu-mole/mole)

    """

    # Calculate
    z = math.acos(1.0/airmass)

    n0 = ref_index.ciddor(wave=w0/10., t=temperature, p=pressure_pa,
                          rh=humidity, co2=co2)
    n1 = ref_index.ciddor(wave=w1/10., t=temperature, p=pressure_pa,
                          rh=humidity, co2=co2)

    return 206265.0 * (n0 - n1) * math.tan(z)


class CorrectDar(BasePrimitive):
    """Correct for Differential Atmospheric Refraction"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """Checks if DAR correction is appropriate"""
        self.logger.info("Checking precondition for CorrectDar")
        # Check image
        if 'GEOMCOR' not in self.action.args.ccddata.header:
            self.logger.error("Can only correct DAR on geometrically corrected "
                              "images")
            self.action.args.ccddata.header['DARCOR'] = (False,
                                                         'DAR corrected?')
            return False
        else:
            if not self.action.args.ccddata.header['GEOMCOR']:
                self.logger.error(
                    "Can only correct DAR on geometrically corrected "
                    "images")
                self.action.args.ccddata.header['DARCOR'] = (False,
                                                             'DAR corrected?')
                return False
            else:
                return True

    def _perform(self):
        """Correct for differential atmospheric refraction"""
        self.logger.info("Correcting for DAR")
        # image size
        image_size = self.action.args.ccddata.data.shape

        # get wavelengths
        w0 = self.action.args.ccddata.header['CRVAL3']
        dw = self.action.args.ccddata.header['CD3_3']
        waves = w0 + np.arange(image_size[0]) * dw
        wgoo0 = self.action.args.ccddata.header['WAVGOOD0']
        wgoo1 = self.action.args.ccddata.header['WAVGOOD1']
        wref = self.action.args.ccddata.header['WAVMID']
        self.logger.info("Ref WL = %.1f, good WL range = (%.1f - %.1f" %
                         (wref, wgoo0, wgoo1))

        # spatial scales in arcsec/item
        y_scale = self.action.args.ccddata.header['PXSCL'] * 3600.
        x_scale = self.action.args.ccddata.header['SLSCL'] * 3600.

        # padding depends on grating
        if 'H' in self.action.args.grating:
            padding_as = 2.0
        elif 'M' in self.action.args.grating:
            padding_as = 3.0
        else:
            padding_as = 4.0

        padding_x = int(padding_as / x_scale)
        padding_y = int(padding_as / y_scale)

        # update WCS
        crpix1 = self.action.args.ccddata.header['CRPIX1']
        crpix2 = self.action.args.ccddata.header['CRPIX2']
        self.action.args.ccddata.header['CRPIX1'] = crpix1 + float(padding_x)
        self.action.args.ccddata.header['CRPIX2'] = crpix2 + float(padding_y)

        # airmass
        airmass = self.action.args.ccddata.header['AIRMASS']
        self.logger.info("Airmass: %.3f" % airmass)

        # IFU orientation
        ifu_pa = self.action.args.ccddata.header['IFUPA']

        # Parallactic angle
        parallactic_angle = self.action.args.ccddata.header['PARANG']

        # Projection angle in radians
        projection_angle_deg = ifu_pa - parallactic_angle
        projection_angle = math.radians(projection_angle_deg)
        self.logger.info("DAR Angles: ifu_pa, parang, projang (deg): "
                         "%.2f, %.2f, %.2f" % (ifu_pa, parallactic_angle,
                                               projection_angle_deg))

        # dispersion over goo wl range in arcsec
        dispersion_max_as = atm_disper(wgoo1, wgoo0, airmass)

        # projected onto IFU
        xdmax_as = dispersion_max_as * math.sin(projection_angle)
        ydmax_as = dispersion_max_as * math.cos(projection_angle)
        self.logger.info("DAR over GOOD WL range: total, x, y (asec): "
                         "%.2f, %.2f, %.2f" % (dispersion_max_as, xdmax_as,
                                               ydmax_as))

        # now in pixels
        xdmax_px = xdmax_as / x_scale
        ydmax_px = ydmax_as / y_scale
        dmax_px = math.sqrt(xdmax_px**2 + ydmax_px**2)
        self.logger.info("DAR over GOOD WL range: total, x, y (pix): "
                         "%.2f, %.2f, %.2f" % (dmax_px, xdmax_px, ydmax_px))

        # prepare output cubes
        output_image = np.zeros((image_size[0], image_size[1]+2*padding_y,
                                 image_size[2]+2*padding_x), dtype=np.float64)
        output_stddev = output_image.copy()
        output_mask = np.zeros((image_size[0], image_size[1]+2*padding_y,
                                image_size[2]+2*padding_x), dtype=np.uint8)
        output_flags = np.zeros((image_size[0], image_size[1] + 2 * padding_y,
                                 image_size[2] + 2 * padding_x), dtype=np.uint8)
        # DAR padded pixel flag
        output_flags += 128

        output_image[:, padding_y:(padding_y+image_size[1]),
                     padding_x:(padding_x+image_size[2])] = \
            self.action.args.ccddata.data

        output_stddev[:, padding_y:(padding_y+image_size[1]),
                      padding_x:(padding_x+image_size[2])] = \
            self.action.args.ccddata.uncertainty.array

        output_mask[:, padding_y:(padding_y+image_size[1]),
                    padding_x:(padding_x+image_size[2])] = \
            self.action.args.ccddata.mask

        output_flags[:, padding_y:(padding_y+image_size[1]),
                     padding_x:(padding_x+image_size[2])] = \
            self.action.args.ccddata.flags

        # check for obj, sky cubes
        output_obj = None
        output_sky = None
        if self.action.args.nasmask and self.action.args.numopen > 1:
            ofn = self.action.args.ccddata.header['OFNAME']

            objfn = ofn.split('.')[0] + '_ocube.fits'
            full_path = os.path.join(
                os.path.dirname(self.action.args.name),
                self.config.instrument.output_directory, objfn)
            if os.path.exists(full_path):
                obj = kcwi_fits_reader(full_path)[0]
                output_obj = np.zeros(
                    (image_size[0], image_size[1] + 2 * padding_y,
                     image_size[2] + 2 * padding_x), dtype=np.float64)
                output_obj[:, padding_y:(padding_y + image_size[1]),
                           padding_x:(padding_x + image_size[2])] = obj.data

            skyfn = ofn.split('.')[0] + '_scube.fits'
            full_path = os.path.join(
                os.path.dirname(self.action.args.name),
                self.config.instrument.output_directory, skyfn)
            if os.path.exists(full_path):
                sky = kcwi_fits_reader(full_path)[0]
                output_sky = np.zeros(
                    (image_size[0], image_size[1] + 2 * padding_y,
                     image_size[2] + 2 * padding_x), dtype=np.float64)
                output_sky[:, padding_y:(padding_y + image_size[1]),
                           padding_x:(padding_x + image_size[2])] = sky.data

        # check if we have a standard star observation
        output_del = None
        stdfile, _ = kcwi_get_std(self.action.args.ccddata.header['OBJECT'],
                                  self.logger)
        if stdfile is not None:
            afn = self.action.args.ccddata.header['ARCFL']

            delfn = afn.split('.')[0] + '_dcube.fits'
            full_path = os.path.join(
                os.path.dirname(self.action.args.name),
                self.config.instrument.output_directory, delfn)
            if os.path.exists(full_path):
                dew = kcwi_fits_reader(full_path)[0]
                output_del = np.zeros(
                    (image_size[0], image_size[1] + 2 * padding_y,
                     image_size[2] + 2 * padding_x), dtype=np.float64)
                output_del[:, padding_y:(padding_y + image_size[1]),
                           padding_x:(padding_x + image_size[2])] = dew.data

        # Perform correction
        for j, wl in enumerate(waves):
            dispersion_correction = atm_disper(wref, wl, airmass)
            x_shift = dispersion_correction * \
                math.sin(projection_angle) / x_scale
            y_shift = dispersion_correction * \
                math.cos(projection_angle) / y_scale
            output_image[j, :, :] = shift(output_image[j, :, :], (y_shift,
                                                                  x_shift))
            output_stddev[j, :, :] = shift(output_stddev[j, :, :], (y_shift,
                                                                    x_shift))
            output_mask[j, :, :] = shift(output_mask[j, :, :], (y_shift,
                                                                x_shift))
            output_flags[j, :, :] = shift(output_flags[j, :, :], (y_shift,
                                                                  x_shift))
        # for obj, sky if they exist
        if output_obj is not None:
            for j, wl in enumerate(waves):
                dispersion_correction = atm_disper(wref, wl, airmass)
                x_shift = dispersion_correction * \
                    math.sin(projection_angle) / x_scale
                y_shift = dispersion_correction * \
                    math.cos(projection_angle) / y_scale
                output_obj[j, :, :] = shift(output_obj[j, :, :], (y_shift,
                                                                  x_shift))

        if output_sky is not None:
            for j, wl in enumerate(waves):
                dispersion_correction = atm_disper(wref, wl, airmass)
                x_shift = dispersion_correction * \
                    math.sin(projection_angle) / x_scale
                y_shift = dispersion_correction * \
                    math.cos(projection_angle) / y_scale
                output_sky[j, :, :] = shift(output_sky[j, :, :], (y_shift,
                                                                  x_shift))

        # for delta wavelength cube, if it exists
        if output_del is not None:
            for j, wl in enumerate(waves):
                dispersion_correction = atm_disper(wref, wl, airmass)
                x_shift = dispersion_correction * \
                    math.sin(projection_angle) / x_scale
                y_shift = dispersion_correction * \
                    math.cos(projection_angle) / y_scale
                output_del[j, :, :] = shift(output_del[j, :, :], (y_shift,
                                                                  x_shift))

        self.action.args.ccddata.data = output_image
        self.action.args.ccddata.uncertainty.array = output_stddev
        self.action.args.ccddata.mask = output_mask
        self.action.args.ccddata.flags = output_flags

        log_string = CorrectDar.__module__

        # update header
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.action.args.ccddata.header['DARCOR'] = (True, 'DAR corrected?')
        self.action.args.ccddata.header['DARANG'] = (projection_angle_deg,
                                                     'DAR projection angle '
                                                     '(deg)')
        self.action.args.ccddata.header['DARPADX'] = (padding_x,
                                                      'DAR X padding (pix)')
        self.action.args.ccddata.header['DARPADY'] = (padding_y,
                                                      'DAR Y padding (pix)')
        self.action.args.ccddata.header['DAREFWL'] = (wref,
                                                      'DAR reference wl (Ang)')
        # write out corrected image
        kcwi_fits_writer(self.action.args.ccddata,
                         table=self.action.args.table,
                         output_file=self.action.args.name,
                         output_dir=self.config.instrument.output_directory,
                         suffix="icubed")
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix="icubed")
        self.context.proctab.write_proctab()

        # check for sky, obj cube
        if output_obj is not None:
            out_obj = CCDData(output_obj,
                              meta=self.action.args.ccddata.header,
                              unit=self.action.args.ccddata.unit)
            kcwi_fits_writer(
                out_obj, output_file=self.action.args.name,
                output_dir=self.config.instrument.output_directory,
                suffix="ocubed")
        if output_sky is not None:
            out_sky = CCDData(output_sky,
                              meta=self.action.args.ccddata.header,
                              unit=self.action.args.ccddata.unit)
            kcwi_fits_writer(
                out_sky, output_file=self.action.args.name,
                output_dir=self.config.instrument.output_directory,
                suffix="scubed")

        # check for delta wave cube
        if output_del is not None:
            out_del = CCDData(output_del,
                              meta=self.action.args.ccddata.header,
                              unit=self.action.args.ccddata.unit)
            kcwi_fits_writer(
                out_del, output_file=self.action.args.name,
                output_dir=self.config.instrument.output_directory,
                suffix="dcubed")

        self.logger.info(log_string)

        return self.action.args
    # END: class CorrectDar()
