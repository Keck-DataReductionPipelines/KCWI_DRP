from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer

import numpy as np
import scipy as sp
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

        logstr = CorrectDar.__module__ + "." + CorrectDar.__qualname__

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
        sz = self.action.args.ccddata.data.shape

        # get wavelengths
        w0 = self.action.args.ccddata.header['CRVAL3']
        dw = self.action.args.ccddata.header['CD3_3']
        waves = w0 + np.arange(sz[0]) * dw
        wgoo0 = self.action.args.ccddata.header['WAVGOOD0']
        wgoo1 = self.action.args.ccddata.header['WAVGOOD1']
        wref = self.action.args.ccddata.header['WAVMID']
        self.logger.info("Ref WL = %.1f, good WL range = (%.1f - %.1f" %
                         (wref, wgoo0, wgoo1))

        # spatial scales in arcsec/item
        yscl = self.action.args.ccddata.header['PXSCL'] * 3600.
        xscl = self.action.args.ccddata.header['SLSCL'] * 3600.

        # padding depends on grating
        if 'H' in self.action.args.grating:
            pad_as = 2.0
        elif 'M' in self.action.args.grating:
            pad_as = 3.0
        else:
            pad_as = 4.0

        pad_x = int(pad_as / xscl)
        pad_y = int(pad_as / yscl)

        # update WCS
        crpix1 = self.action.args.ccddata.header['CRPIX1']
        crpix2 = self.action.args.ccddata.header['CRPIX2']
        self.action.args.ccddata.header['CRPIX1'] = crpix1 + float(pad_x)
        self.action.args.ccddata.header['CRPIX2'] = crpix2 + float(pad_y)

        # airmass
        air = self.action.args.ccddata.header['AIRMASS']
        self.logger.info("Airmass: %.3f" % air)

        # IFU orientation
        ifupa = self.action.args.ccddata.header['IFUPA']

        # Parallactic angle
        parang = self.action.args.ccddata.header['PARANG']

        # Projection angle in radians
        projang_deg = ifupa - parang
        projang = math.radians(projang_deg)
        self.logger.info("DAR Angles: ifu_pa, parang, projang (deg): "
                         "%.2f, %.2f, %.2f" % (ifupa, parang, projang_deg))

        # dispersion over goo wl range in arcsec
        dmax_as = atm_disper(wgoo1, wgoo0, air)

        # projected onto IFU
        xdmax_as = dmax_as * math.sin(projang)
        ydmax_as = dmax_as * math.cos(projang)
        self.logger.info("DAR over GOOD WL range: total, x, y (asec): "
                         "%.2f, %.2f, %.2f" % (dmax_as, xdmax_as, ydmax_as))

        # now in pixels
        xdmax_px = xdmax_as / xscl
        ydmax_px = ydmax_as / yscl
        dmax_px = math.sqrt(xdmax_px**2 + ydmax_px**2)
        self.logger.info("DAR over GOOD WL range: total, x, y (pix): "
                         "%.2f, %.2f, %.2f" % (dmax_px, xdmax_px, ydmax_px))

        # prepare output cubes
        img_out = np.zeros((sz[0], sz[1]+2*pad_y, sz[2]+2*pad_x),
                           dtype=np.float64)
        var_out = img_out.copy()
        msk_out = np.zeros((sz[0], sz[1]+2*pad_y, sz[2]+2*pad_x),
                           dtype=np.uint8)

        img_out[:, pad_y:(pad_y+sz[1]), pad_x:(pad_x+sz[2])] = \
            self.action.args.ccddata.data

        var_out[:, pad_y:(pad_y+sz[1]), pad_x:(pad_x+sz[2])] = \
            self.action.args.ccddata.uncertainty.array

        msk_out[:, pad_y:(pad_y+sz[1]), pad_x:(pad_x+sz[2])] = \
            self.action.args.ccddata.mask

        # Perform correction
        self.logger.info(img_out.shape)
        for j, wl in enumerate(waves):
            dcor = atm_disper(wref, wl, air)
            xsh = dcor * math.sin(projang) / xscl
            ysh = dcor * math.cos(projang) / yscl
            img_out[j, :, :] = sp.ndimage.shift(img_out[j, :, :], (ysh, xsh))
            var_out[j, :, :] = sp.ndimage.shift(var_out[j, :, :], (ysh, xsh))
            msk_out[j, :, :] = sp.ndimage.shift(msk_out[j, :, :], (ysh, xsh))

        self.action.args.ccddata.data = img_out
        self.action.args.ccddata.uncertainty.array = var_out
        self.action.args.ccddata.mask = msk_out

        # update header
        self.action.args.ccddata.header['HISTORY'] = logstr
        self.action.args.ccddata.header['DARCOR'] = (True, 'DAR corrected?')
        self.action.args.ccddata.header['DARANG'] = (projang_deg,
                                                     'DAR projection angle '
                                                     '(deg)')
        self.action.args.ccddata.header['DARPADX'] = (pad_x, 'DAR X padding '
                                                             '(pix)')
        self.action.args.ccddata.header['DARPADY'] = (pad_y, 'DAR Y padding '
                                                             '(pix)')
        self.action.args.ccddata.header['DAREFWL'] = (wref, 'DAR reference wl '
                                                            '(Ang)')
        # write out corrected image
        kcwi_fits_writer(self.action.args.ccddata,
                         table=self.action.args.table,
                         output_file=self.action.args.name, suffix="icubed")
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix="icubed")
        self.context.proctab.write_proctab()

        self.logger.info(logstr)

        return self.action.args
    # END: class CorrectDar()


