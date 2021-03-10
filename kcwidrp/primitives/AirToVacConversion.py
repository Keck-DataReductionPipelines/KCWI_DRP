import os

import numpy as np
from scipy.interpolate import interp1d
from astropy import units as u

from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer, \
                                                    kcwi_fits_reader


class AirToVacConversion(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        return True

    def _perform(self):

        if self.config.instrument.air_to_vacuum:
            self.logger.info("Performing Air to Vacuum Conversion")
            ofn = self.action.args.ccddata.header['OFNAME']
            objfn = ofn.split('.')[0] + '_icube.fits'
            full_path = os.path.join(
                os.path.dirname(self.action.args.name),
                self.config.instrument.output_directory, objfn)
            if os.path.exists(full_path):
                obj = kcwi_fits_reader(full_path)[0]
                obj = self.air2vac(obj)
                kcwi_fits_writer(obj,
                                table=self.action.args.table,
                                output_file=objfn,
                                output_dir=self.config.instrument.output_directory,
                                suffix="vac")
                self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                                    suffix="_vac")
                self.context.proctab.write_proctab()
            else:
                self.logger.warn(f"No file found at location {full_path}")
            
            return self.action.args
        
        self.logger.info("Already in air wavelengths.")
        return self.action.args

    def air2vac(self, hdu, mask=False):
        """Covert wavelengths in a cube from standard air to vacuum.

        Args:
            fits_in (astropy HDU / HDUList): Input HDU/HDUList with 3D data.
            mask (bool): Set if the cube is a mask cube.

        Returns:
            HDU / HDUList*: Trimmed FITS object with updated header.
            *Return type matches type of fits_in argument. 
        """

        cube = np.nan_to_num(hdu.data, nan=0, posinf=0, neginf=0)

        if hdu.header['CTYPE3'] == 'WAVE':
            self.logger.warn("FITS already in vacuum wavelength.")
            return

        wave_air = self.get_wav_axis(hdu.header) * u.nm
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

        hdu.header['CTYPE3'] = 'WAVE'
        hdu.data = cube_new
        return hdu

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

    def get_wav_axis(self, header):
        """Returns a NumPy array representing the wavelength axis of a cube.

        Adapted from https://github.com/dbosul/cwitools.git

        Args:
            header (astropy.io.fits.Header): header that contains wavelength
                or velocity axis that is specified in 'CTYPE' keywords in any 
                dimension.

        Returns:
            numpy.ndarray: Wavelength axis for this data.

        Examples:

            If you wanted to plot your spectrum vs. wavelength in matplotlib:

            >>> import matplotlib.pyplot as plt
            >>> from cwitools import cubes
            >>> from astropy.io import fits
            >>> spec,header = fits.getdata("myspectrum.fits",header=True)
            >>> wav_axis = cubes.get_wav_axis(header)
            >>> fig,ax = plt.subplots(1,1)
            >>> ax.plot(wav_axis,spec)
            >>> fig.show()

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
        if flag == False:
            self.logger.error("Header must contain a wavelength/velocity axis.")

        #Get keywords defining wavelength axis
        nwav = header["NAXIS{0}".format(axis)]
        wav0 = header["CRVAL{0}".format(axis)]
        dwav = header["CD{0}_{0}".format(axis)]
        pix0 = header["CRPIX{0}".format(axis)]

        #Calculate and return
        return np.array([wav0 + (i - pix0) * dwav for i in range(nwav)])
