from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer

import numpy as np
from scipy.ndimage import shift
import math
import ref_index


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

    def _perform(self):
        self.logger.info("Correcting for DAR")

        log_string = CorrectDar.__module__

        # Check image
        if 'GEOMCOR' not in self.action.args.ccddata.header:
            self.logger.error("Can only correct DAR on geometrically corrected "
                              "images")
            self.action.args.ccddata.header['DARCOR'] = (False,
                                                         'DAR corrected?')
            return self.action.args
        else:
            if not self.action.args.ccddata.header['GEOMCOR']:
                self.logger.error(
                    "Can only correct DAR on geometrically corrected "
                    "images")
                self.action.args.ccddata.header['DARCOR'] = (False,
                                                             'DAR corrected?')
                return self.action.args

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
        output_variance = output_image.copy()
        output_mask = np.zeros((image_size[0], image_size[1]+2*padding_y,
                                image_size[2]+2*padding_x), dtype=np.uint8)

        output_image[:, padding_y:(padding_y+image_size[1]),
                     padding_x:(padding_x+image_size[2])] = \
            self.action.args.ccddata.data

        output_variance[:, padding_y:(padding_y+image_size[1]),
                        padding_x:(padding_x+image_size[2])] = \
            self.action.args.ccddata.uncertainty.array

        output_mask[:, padding_y:(padding_y+image_size[1]),
                    padding_x:(padding_x+image_size[2])] = \
            self.action.args.ccddata.mask

        # Perform correction
        for j, wl in enumerate(waves):
            dispersion_correction = atm_disper(wref, wl, airmass)
            x_shift = dispersion_correction * \
                math.sin(projection_angle) / x_scale
            y_shift = dispersion_correction * \
                math.cos(projection_angle) / y_scale
            output_image[j, :, :] = shift(output_image[j, :, :], (y_shift,
                                                                  x_shift))
            output_variance[j, :, :] = shift(output_variance[j, :, :],
                                             (y_shift, x_shift))
            output_mask[j, :, :] = shift(output_mask[j, :, :], (y_shift,
                                                                x_shift))

        self.action.args.ccddata.data = output_image
        self.action.args.ccddata.uncertainty.array = output_variance
        self.action.args.ccddata.mask = output_mask

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
                         output_file=self.action.args.name, suffix="icubed")
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix="icubed")
        self.context.proctab.write_proctab()

        self.logger.info(log_string)

        return self.action.args
    # END: class CorrectDar()
