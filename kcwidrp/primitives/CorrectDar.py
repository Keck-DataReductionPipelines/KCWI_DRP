from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.kcwi_get_std import kcwi_get_std
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer, \
    kcwi_fits_reader, strip_fname

import numpy as np
from scipy.ndimage import shift
from astropy.nddata import CCDData
import math
import ref_index
import os
import erfa


def atm_disper(w0, w1, airmass, temperature=10.0, pressure_pa=61100.0,
               humidity=50.0, co2=400.0):
    """

    Calculate atmospheric dispersion at w1 relative to w0

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


def atm_disper_slalib(w0, w1, airmass, temperature=10.0, pressure_pa=61100.0,
                    humidity=50.0, co2=400.0):
    z = math.acos(1.0/airmass)
    HM = 4160
    TDK = temperature + 273.15
    PMB = pressure_pa/100.
    RH = humidity/100.
    PHI = math.radians(19.82833333333333)   #Keck latitude
    TLR = 0.0065
    EPS = 1.0e-8
    R0 = 206265.*refro(z, HM, TDK, PMB, RH, w0/1e4, PHI, TLR, EPS)
    R1 = 206265.*refro(z, HM, TDK, PMB, RH, w1/1e4, PHI, TLR, EPS)
    return R0-R1


def atm_disper_ERFA(w0, w1, airmass, temperature=10.0, pressure_pa=61100.0,
                    humidity=50.0, co2=400.0):
    z = math.acos(1.0/airmass)
    pmb = pressure_pa/100.
    rh = humidity/100.
    refA0, refB0 = erfa.refco(pmb, temperature, rh, w0/1e4)
    R0 = 206265*(refA0*np.tan(z) + refB0*np.tan(z)**3)
    refA1, refB1 = erfa.refco(pmb, temperature, rh, w1/1e4)
    R1 = 206265*(refA1*np.tan(z) + refB1*np.tan(z)**3)
    return R0-R1


class CorrectDar(BasePrimitive):
    """
    Correct for Differential Atmospheric Refraction

    Accounts for rotator orientation, zenith angle, and parallactic angle to
    correct input data cube into a padded, DAR corrected output cube.
    Calculates the DAR correction for each wavelength slice and adjusts position
    in cube using scipy.nddata.shift to implement correction.

    Sets flags in padded region of output cube to a value of 128.

    Also corrects delta wavelength cube if it is a standard star observation.

    """

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

        # Determine which DAR correction model to mauseke
        DAR_correction_model = self.config.instrument.dar_correction_model
        options = ["Filippenko1982", "slalib", "ERFA"]

        # If the config file has an invalid option, return
        if not bool([el for el in options if el in DAR_correction_model]):
            self.logger.error(f'Bad config option for dar_correction_model\
                , options are {options}')
            return self.action.args

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
        self.logger.info("Ref WL = %.1f, good WL range = (%.1f - %.1f)" %
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
        if DAR_correction_model == "Filippenko1982":
            dispersion_max_as = atm_disper(wgoo1, wgoo0, airmass)
        elif DAR_correction_model == "slalib":
            dispersion_max_as = atm_disper_slalib(wgoo1, wgoo0, airmass)
        elif DAR_correction_model == "ERFA":
            dispersion_max_as = atm_disper_ERFA(wgoo1, wgoo0, airmass)

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
        output_flags = np.zeros((image_size[0], image_size[1] + 2*padding_y,
                                 image_size[2] + 2 * padding_x), dtype=np.uint8)
        if self.action.args.ccddata.noskysub is not None:
            output_noskysub = np.zeros((image_size[0],
                                        image_size[1] + 2*padding_y,
                                        image_size[2] + 2 * padding_x),
                                       dtype=np.float64)
        else:
            output_noskysub = None
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

        if output_noskysub is not None:
            output_noskysub[:, padding_y:(padding_y + image_size[1]),
                            padding_x:(padding_x + image_size[2])] = \
                self.action.args.ccddata.noskysub

        # check for obj, sky cubes
        output_obj = None
        output_sky = None
        if self.action.args.nasmask and self.action.args.numopen > 1:
            ofn = self.action.args.name

            objfn = strip_fname(ofn) + '_ocube.fits'
            full_path = os.path.join(
                self.config.instrument.cwd,
                self.config.instrument.output_directory, objfn)
            if os.path.exists(full_path):
                obj = kcwi_fits_reader(full_path)[0]
                output_obj = np.zeros(
                    (image_size[0], image_size[1] + 2 * padding_y,
                     image_size[2] + 2 * padding_x), dtype=np.float64)
                output_obj[:, padding_y:(padding_y + image_size[1]),
                           padding_x:(padding_x + image_size[2])] = obj.data

            skyfn = strip_fname(ofn) + '_scube.fits'
            full_path = os.path.join(
                self.config.instrument.cwd,
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

            delfn = strip_fname(afn) + '_dcube.fits'
            full_path = os.path.join(
                self.config.instrument.cwd,
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
            output_image[j, :, :] = shift(output_image[j, :, :],
                                          (y_shift, x_shift))
            output_stddev[j, :, :] = shift(output_stddev[j, :, :],
                                           (y_shift, x_shift))
            output_mask[j, :, :] = shift(output_mask[j, :, :],
                                         (y_shift, x_shift))
            output_flags[j, :, :] = shift(output_flags[j, :, :],
                                          (y_shift, x_shift))
            if output_noskysub is not None:
                output_noskysub[j, :, :] = shift(output_noskysub[j, :, :],
                                                 (y_shift, x_shift))
        # for obj, sky if they exist
        if output_obj is not None:
            for j, wl in enumerate(waves):
                dispersion_correction = atm_disper(wref, wl, airmass)
                x_shift = dispersion_correction * \
                    math.sin(projection_angle) / x_scale
                y_shift = dispersion_correction * \
                    math.cos(projection_angle) / y_scale
                output_obj[j, :, :] = shift(output_obj[j, :, :],
                                            (y_shift, x_shift))

        if output_sky is not None:
            for j, wl in enumerate(waves):
                dispersion_correction = atm_disper(wref, wl, airmass)
                x_shift = dispersion_correction * \
                    math.sin(projection_angle) / x_scale
                y_shift = dispersion_correction * \
                    math.cos(projection_angle) / y_scale
                output_sky[j, :, :] = shift(output_sky[j, :, :],
                                            (y_shift, x_shift))

        # for delta wavelength cube, if it exists
        if output_del is not None:
            for j, wl in enumerate(waves):
                dispersion_correction = atm_disper(wref, wl, airmass)
                x_shift = dispersion_correction * \
                    math.sin(projection_angle) / x_scale
                y_shift = dispersion_correction * \
                    math.cos(projection_angle) / y_scale
                output_del[j, :, :] = shift(output_del[j, :, :],
                                            (y_shift, x_shift))

        self.action.args.ccddata.data = output_image
        self.action.args.ccddata.uncertainty.array = output_stddev
        self.action.args.ccddata.mask = output_mask
        self.action.args.ccddata.flags = output_flags
        self.action.args.ccddata.noskysub = output_noskysub

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
                                            suffix="icubed",
                                            filename=self.action.args.name)
        self.context.proctab.write_proctab(tfil=self.config.instrument.procfile)

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







'''The code for the refro, atmt, atms, and da1ps functions was pulled
from the Keck AO Imaging (KAI) data reduction pipeline at
https://github.com/Keck-DataReductionPipelines/KAI

The refco function appears to be a port of the same ERFA function here:
http://www.iausofa.org/2023_1011_F/Astrometry.html

Code below this point is covered by the following license:

Copyright 2021 J.R. Lu, A. Gautam, T. Do

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import math


def refro(ZOBS, HM, TDK, PMB, RH, WL, PHI, TLR, EPS):
    """
    Atmospheric refraction for radio and optical/IR wavelengths.
  
    Given:
      ZOBS    d  observed zenith distance of the source (radian)
      HM      d  height of the observer above sea level (metre)
      TDK     d  ambient temperature at the observer (deg K)
      PMB     d  pressure at the observer (millibar)
      RH      d  relative humidity at the observer (range 0-1)
      WL      d  effective wavelength of the source (micrometre)
      PHI     d  latitude of the observer (radian, astronomical)
      TLR     d  temperature lapse rate in the troposphere (degK/metre)
      EPS     d  precision required to terminate iteration (radian)
  
    Returned:
      REF     d  refraction: in vacuo ZD minus observed ZD (radian)
  
    Notes:
  
    1  A suggested value for the TLR argument is 0.0065.  The
       refraction is significantly affected by TLR, and if studies
       of the local atmosphere have been carried out a better TLR
       value may be available.
  
    2  A suggested value for the EPS argument is 1D-8.  The result is
       usually at least two orders of magnitude more computationally
       precise than the supplied EPS value.
  
    3  The routine computes the refraction for zenith distances up
       to and a little beyond 90 deg using the method of Hohenkerk
       and Sinclair (NAO Technical Notes 59 and 63, subsequently adopted
       in the Explanatory Supplement, 1992 edition - see section 3.281).
  
    4  The code is a development of the optical/IR refraction subroutine
       AREF of C.Hohenkerk (HMNAO, September 1984), with extensions to
       support the radio case.  Apart from merely cosmetic changes, the
       following modifications to the original HMNAO optical/IR refraction
       code have been made:
  
       .  The angle arguments have been changed to radians.
  
       .  Any value of ZOBS is allowed (see note 6, below).
  
       .  Other argument values have been limited to safe values.
  
       .  Murray's values for the gas constants have been used
          (Vectorial Astrometry, Adam Hilger, 1983).
  
       .  The numerical integration phase has been rearranged for
          extra clarity.
  
       .  A better model for Ps(T) has been adopted (taken from
          Gill, Atmosphere-Ocean Dynamics, Academic Press, 1982).
  
       .  More accurate expressions for Pwo have been adopted
          (again from Gill 1982).
  
       .  Provision for radio wavelengths has been added using
          expressions devised by A.T.Sinclair, RGO (private
          communication 1989), based on the Essen & Froome
          refractivity formula adopted in Resolution 1 of the
          13th International Geodesy Association General Assembly
          (Bulletin Geodesique 70 p390, 1963).
  
       .  Various small changes have been made to gain speed.
  
       None of the changes significantly affects the optical/IR results
       with respect to the algorithm given in the 1992 Explanatory
       Supplement.  For example, at 70 deg zenith distance the present
       routine agrees with the ES algorithm to better than 0.05 arcsec
       for any reasonable combination of parameters.  However, the
       improved water-vapour expressions do make a significant difference
       in the radio band, at 70 deg zenith distance reaching almost
       4 arcsec for a hot, humid, low-altitude site during a period of
       low pressure.
  
    5  The radio refraction is chosen by specifying WL > 100 micrometres.
       Because the algorithm takes no account of the ionosphere, the
       accuracy deteriorates at low frequencies, below about 30 MHz.
  
    6  Before use, the value of ZOBS is expressed in the range +/- pi.
       If this ranged ZOBS is -ve, the result REF is computed from its
       absolute value before being made -ve to match.  In addition, if
       it has an absolute value greater than 93 deg, a fixed REF value
       equal to the result for ZOBS = 93 deg is returned, appropriately
       signed.
  
    7  As in the original Hohenkerk and Sinclair algorithm, fixed values
       of the water vapour polytrope exponent, the height of the
       tropopause, and the height at which refraction is negligible are
       used.
  
    8  The radio refraction has been tested against work done by
       Iain Coulson, JACH, (private communication 1995) for the
       James Clerk Maxwell Telescope, Mauna Kea.  For typical conditions,
       agreement at the 0.1 arcsec level is achieved for moderate ZD,
       worsening to perhaps 0.5-1.0 arcsec at ZD 80 deg.  At hot and
       humid sea-level sites the accuracy will not be as good.
  
    9  It should be noted that the relative humidity RH is formally
       defined in terms of "mixing ratio" rather than pressures or
       densities as is often stated.  It is the mass of water per unit
       mass of dry air divided by that for saturated air at the same
       temperature and pressure (see Gill 1982).
  
    Called:  slDA1P, slATMT, slATMS
  
    P.T.Wallace   Starlink   3 June 1997
  
    Copyright (C) 1997 Rutherford Appleton Laboratory
    Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
    """
    # Fixed parameters
    
    # 93 degrees in radians
    D93 = 1.623156204
    # Universal gas constant
    GCR = 8314.32
    # Molecular weight of dry air
    DMD = 28.9644
    # Molecular weight of water vapour
    DMW = 18.0152
    # Mean Earth radius (metre)
    S = 6378120.0
    # Exponent of temperature dependence of water vapour pressure
    DELTA = 18.36
    # Height of tropopause (metre)
    HT = 11000.0
    # Upper limit for refractive effects (metre)
    HS = 80000.0
    
    # The refraction integrand
    def refi(r, dn, rdndr):
        return rdndr / (dn + rdndr)

    # Transform ZOBS into the normal range.
    ZOBS1 = da1p(ZOBS)
    ZOBS2 = min(abs(ZOBS1), 1.0e93)

    # Keep other arguments within safe bounds.
    HMOK = min(max(HM, -1.0e3), 10.0e3)
    TDKOK = min(max(TDK, 100.0), 500.0)
    PMBOK = min(max(PMB, 0.0), 10000.0)
    RHOK = min(max(RH, 0.0), 1.0)
    WLOK = max(WL, 0.1)
    ALPHA = min(max(abs(TLR), 0.001), 0.01)

    # Tolerance for iteration.
    TOL = min(max(abs(EPS), 1.0e-12), 0.1) / 2.0

    # Decide whether optical/IR or radio case - switch at 100 microns.
    OPTIC = WLOK <= 100.0

    # Set up model atmosphere parameters defined at the observer.
    WLSQ = WLOK * WLOK
    GB = 9.784 * (1.0 - 0.0026 * math.cos(PHI + PHI) - 0.00000028 * HMOK)
    if OPTIC:
        A = (287.604 + (1.6288 + 0.0136 / WLSQ) / WLSQ) * 273.15e-6 / 1013.25
    else:
        A = 77.624e-6
    
    GAMAL = (GB * DMD) / GCR
    GAMMA = GAMAL / ALPHA
    GAMM2 = GAMMA - 2.0
    DELM2 = DELTA - 2.0
    TDC = TDKOK - 273.15
    PSAT = 10.0**((0.7859 + 0.03477 * TDC) / (1.0 + 0.00412 * TDC))
    PSAT *= (1.0 + PMBOK * (4.5e-6 + 6e-10 * TDC * TDC))
    if (PMBOK > 0.0):
        PWO = RHOK * PSAT / (1.0 - (1.0 - RHOK) * PSAT / PMBOK)
    else:
        PWO = 0.0

    W = PWO * (1.0 - DMW / DMD) * GAMMA / (DELTA - GAMMA)
    C1 = A * (PMBOK + W) / TDKOK
    if OPTIC:
        C2 = (A * W + 11.2684e-6 * PWO) / TDKOK
    else:
        C2 = (A * W + 12.92e-6 * PWO) / TDKOK

    C3 = (GAMMA - 1.0) * ALPHA * C1 / TDKOK
    C4 = (DELTA - 1.0) * ALPHA * C2 / TDKOK
    if OPTIC:
        C5 = 0.0
        C6 = 0.0
    else:
        C5 = 371897e-6 * PWO / TDKOK
        C6 = C5 * DELM2 * ALPHA / (TDKOK * TDKOK)

    # Conditions at the observer.
    R0 = S + HMOK
    TEMPO, DN0, RDNDR0 = atmt(R0, TDKOK, ALPHA, GAMM2, DELM2,
                              C1, C2, C3, C4, C5, C6, R0)
                              
    SK0 = DN0 * R0 * math.sin(ZOBS2)
    F0 = refi(R0, DN0, RDNDR0)

    # Conditions in the troposphere at the tropopause.
    RT = S + HT
    TT, DNT, RDNDRT = atmt(R0, TDKOK, ALPHA, GAMM2, DELM2,
                           C1, C2, C3, C4, C5, C6, RT)
    SINE = SK0 / (RT * DNT)
    ZT = math.atan2(SINE, math.sqrt(max(1.0 - SINE * SINE, 0.0)))
    FT = refi(RT, DNT, RDNDRT)

    # Conditions in the stratosphere at the tropopause.
    DNTS, RDNDRP = atms(RT, TT, DNT, GAMAL, RT)
    SINE = SK0 / (RT * DNTS)
    ZTS = math.atan2(SINE, math.sqrt(max(1.0 - SINE * SINE,0.0)))
    FTS = refi(RT, DNTS, RDNDRP)

    # Conditions at the stratosphere limit.
    RS = S + HS
    DNS, RDNDRS = atms(RT, TT, DNT, GAMAL, RS)
    SINE = SK0 / (RS * DNS)
    ZS = math.atan2(SINE, math.sqrt(max(1.0 - SINE * SINE, 0.0)))
    FS = refi(RS, DNS, RDNDRS)

    # Integrate the refraction integral in two parts;  first in the
    # troposphere (K=1), then in the stratosphere (K=2).

    # Initialize previous refraction to ensure at least two iterations.
    REFOLD = 1.0e6

    # Start off with 8 strips for the troposphere integration, and then
    # use the final troposphere value for the stratosphere integration,
    # which tends to need more strips.
    IS = 8

    # Troposphere then stratosphere.
    for K in [1,2]:
        # Start Z, Z range, and start and end values.
        if K == 1:
            Z0 = ZOBS2
            ZRANGE = ZT - Z0
            FB = F0
            FF = FT
        else:
            Z0 = ZTS
            ZRANGE = ZS - Z0
            FB = FTS
            FF = FS

        # Sums of odd and even values.
        FO = 0.0
        FE = 0.0

        # First time through the loop we have to do every point.
        N = 1
       
        # Start of iteration loop (terminates at specified precision).
        LOOP = True
        while LOOP:
            # Strip width.
            H = ZRANGE / float(IS)

            # Initialize distance from Earth centre for quadrature pass.
            if K == 1:
                R = R0
            else:
                R = RT

            # One pass (no need to compute evens after first time).
            for I in range(1, IS, N):
                # Sine of observed zenith distance.
                SZ = math.sin(Z0 + H * float(I))

                # Find R (to the nearest metre, maximum four iterations).
                if SZ > 1e-20:
                    W = SK0 / SZ
                    RG = R
                    DR = 1e6
                    J = 0
                    
                    while (abs(DR) > 1.0) and (J < 4):
                        J = J + 1
                        if K == 1:
                            TG, DN, RDNDR = atmt(R0, TDKOK, ALPHA, GAMM2, DELM2,
                                                 C1, C2, C3, C4, C5, C6, RG)
                        else:
                            DN, RDNDR = atms(RT, TT, DNT, GAMAL, RG)

                        DR = (RG * DN - W) / (DN + RDNDR)
                        RG = RG - DR

                    R = RG

                # Find the refractive index and integrand at R.
                if K == 1:
                    T, DN, RDNDR = atmt(R0, TDKOK, ALPHA, GAMM2, DELM2,
                                        C1, C2, C3, C4, C5, C6, R)
                else:
                    DN,RDNDR = atms(RT, TT, DNT, GAMAL, R)

                F = refi(R, DN, RDNDR)

                # Accumulate odd and (first time only) even values.
                if (N == 1) and ((I % 2) == 0):
                    FE = FE + F
                else:
                    FO = FO + F

            # Evaluate the integrand using Simpson's Rule.
            REFP = H * (FB + 4.0 * FO + 2.0 * FE + FF) / 3.0

            # Has the required precision been achieved?
            if (abs(REFP - REFOLD) > TOL):

                # No: prepare for next iteration.

                # Save current value for convergence test.
                REFOLD = REFP

                # Double the number of strips.
                IS = IS + IS

                # Sum of all current values = sum of next pass's even values.
                FE = FE + FO

                # Prepare for new odd values.
                FO = 0.0

                # Skip even values next time.
                N = 2
                
            else:
                # Yes: save troposphere component and terminate the loop.
                if (K == 1):
                    REFT = REFP

                LOOP = False
            # END IF
        # END FOR
    # END WHILE

    # Result.
    REF = REFT + REFP
    if (ZOBS1 < 0.0):
        REF = -REF

    return REF


def atmt(R0, T0, ALPHA, GAMM2, DELM2, C1, C2, C3, C4, C5, C6, R):
    """
    Internal routine used by REFRO
  
    Refractive index and derivative with respect to height for the
    troposphere.
  
    Given:
      R0      d    height of observer from centre of the Earth (metre)
      T0      d    temperature at the observer (deg K)
      ALPHA   d    alpha          )
      GAMM2   d    gamma minus 2  ) see HMNAO paper
      DELM2   d    delta minus 2  )
      C1      d    useful term  )
      C2      d    useful term  )
      C3      d    useful term  ) see source
      C4      d    useful term  ) of slRFRO
      C5      d    useful term  )
      C6      d    useful term  )
      R       d    current distance from the centre of the Earth (metre)
  
    Returned:
      T       d    temperature at R (deg K)
      DN      d    refractive index at R
      RDNDR   d    R    rate the refractive index is changing at R
  
    Note that in the optical case C5 and C6 are zero.
  
    P.T.Wallace   Starlink   30 May 1997
  
    Copyright (C) 1997 Rutherford Appleton Laboratory
    Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
    """
    T = max(min(T0 - ALPHA * (R - R0), 320.0), 100.0)
    TT0 = T / T0
    TT0GM2 = TT0**GAMM2
    TT0DM2 = TT0**DELM2
    DN = 1.0 + (C1 * TT0GM2 - (C2 - C5 / T) * TT0DM2) * TT0
    RDNDR = R * (-C3 * TT0GM2 + (C4 - C6 / TT0) * TT0DM2)

    return T, DN, RDNDR


def atms(RT, TT, DNT, GAMAL, R):
    """
  
    Internal routine used by REFRO
  
    Refractive index and derivative with respect to height for the
    stratosphere.
  
    Given:
      RT      d    height of tropopause from centre of the Earth (metre)
      TT      d    temperature at the tropopause (deg K)
      DNT     d    refractive index at the tropopause
      GAMAL   d    constant of the atmospheric model = G  MD/R
      R       d    current distance from the centre of the Earth (metre)
  
    Returned:
      DN      d    refractive index at R
      RDNDR   d    R    rate the refractive index is changing at R
  
    P.T.Wallace   Starlink   14 July 1995
  
    Copyright (C) 1995 Rutherford Appleton Laboratory
    Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
    """
    B = GAMAL / TT
    W = (DNT - 1.0) * math.exp(-B * (R - RT))
    DN = 1.0 + W
    RDNDR = -R * B * W

    return DN, RDNDR
        

def da1p(ANGLE):
    """
    Normalize angle into range +/- pi  (double precision)
  
    Given:
       ANGLE     dp      the angle in radians
  
    The result (double precision) is ANGLE expressed in the range +/- pi.
  
    P.T.Wallace   Starlink   23 November 1995
  
    Copyright (C) 1995 Rutherford Appleton Laboratory
    Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
    """
    DPI = 3.141592653589793238462643
    D2PI = 6.283185307179586476925287

    slDA1P = ANGLE % D2PI
    if (abs(slDA1P) >= DPI):
        slDA1P = slDA1P - math.copysign(D2PI, ANGLE)

    return slDA1P
