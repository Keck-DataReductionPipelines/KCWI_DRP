import os

import numpy as np
from scipy.interpolate import interp1d
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation

from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer, \
                                                    kcwi_fits_reader


class WavelengthCorrections(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        return True

    def _perform(self):

        # Determine which radial velocity correction to make
        correction_mode = self.config.instrument.radial_velocity_correction
        options = ["none", "barycentric", "heliocentric"]

        # If the config file has an invalid option, return
        if not bool([el for el in options if el in correction_mode]):
            self.logger.error('Bad config option for radial_velocity_correction\
                , options are ["none", "heliocentric", "barycentric"]')
            return self.action.args
        
        suffix = 'icube' # Can be ammended to handle ocube files
        obj = self.locate_object_file(suffix)

        if "none" in correction_mode:
            self.logger.info("Skipping radial velocity correction")
        
        else:
            self.logger.info(f"Performing {correction_mode} correction")
            obj = self.heliocentric(obj, correction_mode)
        
        if self.config.instrument.air_to_vacuum:
            self.logger.info("Performing Air to Vacuum Conversion")
            obj = self.air2vac(obj)
        
        log_string = WavelengthCorrections.__module__
        obj.header['HISTORY'] = log_string
        
        kcwi_fits_writer(obj,
                        table=self.action.args.table,
                        output_file=self.action.args.name,
                        output_dir=self.config.instrument.output_directory,
                        suffix=f'{suffix}w')
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix=f'_{suffix}w')
        self.context.proctab.write_proctab()

        # Unsure here: Is this right? it seems to make DAR happy
        self.action.args.ccddata = obj
        
        return self.action.args

    def air2vac(self, obj, mask=False):
        """Covert wavelengths in a cube from standard air to vacuum.

        Args:
            fits_in (astropy HDU / HDUList): Input HDU/HDUList with 3D data.
            mask (bool): Set if the cube is a mask cube.

        Returns:
            HDU / HDUList*: Trimmed FITS object with updated header.
            *Return type matches type of fits_in argument. 
        """

        cube = np.nan_to_num(obj.data, nan=0, posinf=0, neginf=0)

        if obj.header['CTYPE3'] == 'WAVE':
            self.logger.warn("FITS already in vacuum wavelength.")
            return

        wave_air = self.get_wav_axis(obj.header) * u.AA
        wave_vac = self.a2v_conversion(wave_air)

        # resample to uniform grid
        cube_new = np.zeros_like(cube)
        for i in range(cube.shape[2]):
            for j in range(cube.shape[1]):

                spec0 = cube[:, j, i]
                if not mask:
                    f_cubic = interp1d(
                        wave_vac,
                        spec0,
                        kind='cubic',
                        fill_value='extrapolate'
                        )
                    spec_new = f_cubic(wave_air)
                else:
                    f_pre = interp1d(
                        wave_vac,
                        spec0,
                        kind='previous',
                        bounds_error=False,
                        fill_value=128
                        )
                    spec_pre = f_pre(wave_air)

                    f_nex = interp1d(
                        wave_vac,
                        spec0,
                        kind='next',
                        bounds_error=False,
                        fill_value=128
                        )
                    spec_nex = f_nex(wave_air)

                    spec_new = np.zeros_like(spec0)
                    for k in range(spec0.shape[0]):
                        spec_new[k] = max(spec_pre[k], spec_nex[k])

                cube_new[:, j, i] = spec_new

        obj.header['CTYPE3'] = 'WAVE'
        obj.data = cube_new
        return obj

    def a2v_conversion(self, wave):
        """ Convert air-based wavelengths to vacuum

        Adapted from wave.py in: https://github.com/pypeit/PypeIt/
        Formula from https://ui.adsabs.harvard.edu/abs/1996ApOpt..35.1566C/

        Parameters
        ----------
        wave: Quantity array
            Wavelengths
        Returns
        -------
        wave: Quantity array
            Wavelength array corrected to vacuum wavelengths
        """
        # Convert to AA
        wave = wave.to(u.AA)
        wavelength = wave.value

        # Standard conversion format
        sigma_sq = (1.e4/wavelength)**2. #wavenumber squared
        factor = 1 + (5.792105e-2/(238.0185-sigma_sq)) + (1.67918e-3/(57.362-sigma_sq))
        factor = factor*(wavelength>=2000.) + 1.*(wavelength<2000.) #only modify above 2000A

        # Convert
        wavelength = wavelength*factor
        # Units
        new_wave = wavelength*u.AA
        new_wave.to(wave.unit)

        return new_wave

    def heliocentric(self, obj, correction_mode, mask=False, resample=True, vcorr=None):
        """Apply heliocentric correction to the cubes. 
        *Note that this only works for KCWI data because the location of the Keck 
        Observatory is hard-coded in the function.*

        Adapted from https://github.com/dbosul/cwitools.git

        Args:
            fits_in (astropy HDU / HDUList): Input HDU/HDUList with 3D data.
            mask (bool): Set if the cube is a mask cube. This only works for
                resampled cubes.
            return_vcorr (bool): If set, return the correction velocity (in km/s)
                as well.
            resample (bool): Resample the cube to the original wavelength grid?
            vcorr (float): Use a different correction velocity.
            barycentric (bool): Use barycentric correction instead of heliocentric.

        Returns:
            HDU / HDUList*: Trimmed FITS object with updated header.
            vcorr (float): (if vcorr is True) Correction velocity in km/s.
            *Return type matches type of fits_in argument.
            
        Examples: 
            
            To apply heliocentric correction,
            
            >>> hdu_new = heliocentric(hdu_old)
            
            However, this resamples the wavelengths back to the original grid. To
            use the new grid without resampling the data,
            
            >>> hdu_new = heliocentric(hdu_old, resample=False)
        """

        barycentric = ("barycentric" in correction_mode)

        cube = np.nan_to_num(obj.data,
                            nan=0, posinf=0, neginf=0)

        v_old = 0.
        if 'VCORR' in obj.header:
            v_old = obj.header['VCORR']
            self.logger.info("Rolling back the existing correction with:")
            self.logger.info("Vcorr = %.2f km/s." % (v_old))

        if vcorr is None:
            targ = SkyCoord(
                obj.header['TARGRA'],
                obj.header['TARGDEC'],
                unit='deg',
                obstime=obj.header['DATE-BEG']
            )

            lat = self.config.instrument.latitude
            lon = self.config.instrument.longitude
            alt = self.config.instrument.altitude

            keck = EarthLocation.from_geodetic(lat=lat, lon=lon, height=alt)
            if barycentric:
                vcorr = targ.radial_velocity_correction(
                    kind='barycentric', location=keck)
            else:
                vcorr = targ.radial_velocity_correction(
                    kind='heliocentric', location=keck)
            vcorr = vcorr.to('km/s').value

        self.logger.info("Helio/Barycentric correction:")
        self.logger.info("Vcorr = %.2f km/s." % (vcorr))

        v_tot = vcorr-v_old

        if not resample:
            obj.header['CRVAL3'] *= (1 + v_tot / 2.99792458e5)
            obj.header['CD3_3'] *= (1 + v_tot / 2.99792458e5)
            obj.header['VCORR'] = vcorr
            return obj

        wav_old = self.get_wav_axis(obj.header)
        wav_hel = wav_old * (1 + v_tot / 2.99792458e5)

        # resample to uniform grid
        self.logger.info("Resampling to uniform grid")
        cube_new = np.zeros_like(cube)
        for i in range(cube.shape[2]):
            for j in range(cube.shape[1]):

                spc0 = cube[:, j, i]
                if not mask:
                    f_cubic = interp1d(wav_hel, spc0, kind='cubic',
                                    fill_value='extrapolate')
                    spec_new = f_cubic(wav_old)

                else:
                    f_pre = interp1d(wav_hel, spc0, kind='previous',
                                    bounds_error=False, fill_value=128)
                    spec_pre = f_pre(wav_old)
                    f_nex = interp1d(wav_hel, spc0, kind='next',
                                    bounds_error=False, fill_value=128)
                    spec_nex = f_nex(wav_old)
                    spec_new = np.zeros_like(spc0)
                    for k in range(spc0.shape[0]):
                        spec_new[k] = max(spec_pre[k], spec_nex[k])
                cube_new[:, j, i] = spec_new
        
        obj.header['VCORR'] = vcorr
        obj.data = cube_new
        return obj

    def get_wav_axis(self, header):
        """Returns a NumPy array representing the wavelength axis of a cube.

        Adapted from https://github.com/dbosul/cwitools.git

        Args:
            header (astropy.io.fits.Header): header that contains wavelength
                or velocity axis that is specified in 'CTYPE' keywords in any 
                dimension.

        Returns:
            numpy.ndarray: Wavelength axis for this data.

        """

        #Select the appropriate axis.
        naxis = header['NAXIS']
        flag = False
        for i in range(naxis):
            #Keyword entry
            card = "CTYPE{0}".format(i+1)
            if not card in header:
                self.logger.warning.error("Header must contain 'CTYPE' keywords.")
            
            #Possible wave types.
            if header[card] in ['AWAV', 'WAVE', 'VELO']:
                axis = i+1
                flag = True
                break

        #No wavelength axis
        if flag is False:
            self.logger.error("Header must contain a wavelength/velocity axis.")

        #Get keywords defining wavelength axis
        nwav = header["NAXIS{0}".format(axis)]
        wav0 = header["CRVAL{0}".format(axis)]
        dwav = header["CD{0}_{0}".format(axis)]
        pix0 = header["CRPIX{0}".format(axis)]

        #Calculate and return
        return np.array([wav0 + (i - pix0) * dwav for i in range(nwav)])
    
    def locate_object_file(self, suffix):
        ofn = self.action.args.ccddata.header['OFNAME']
        objfn = ofn.split('.')[0] + f'_{suffix}.fits'
        full_path = os.path.join(
            os.path.dirname(self.action.args.name),
            self.config.instrument.output_directory, objfn)
        if os.path.exists(full_path):
            return kcwi_fits_reader(full_path)[0]
        else:
            self.logger.error(f'Unable to read file {objfn}')
            return None

    