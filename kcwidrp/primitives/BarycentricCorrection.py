import os

import numpy as np
from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord, EarthLocation

from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer, \
                                                    kcwi_fits_reader


class BarycentricCorrection(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        return True

    def _perform(self):

        self.logger.info("Performing Barycentric Correction")

        ofn = self.action.args.ccddata.header['OFNAME']
        objfn = ofn.split('.')[0] + '_icubed.fits'
        full_path = os.path.join(
            os.path.dirname(self.action.args.name),
            self.config.instrument.output_directory, objfn)
        if os.path.exists(full_path):
            obj = kcwi_fits_reader(full_path)[0]
            obj = self.heliocentric(obj)
            kcwi_fits_writer(obj,
                            table=self.action.args.table,
                            output_file=objfn,
                            output_dir=self.config.instrument.output_directory,
                            suffix="h")
            self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                                suffix="h")
            self.context.proctab.write_proctab()
        else:
            self.logger.warn(f"No file found at location {full_path}")
        
        return self.action.args

    def heliocentric(self, obj, mask=False, resample=True, vcorr=None):
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

        barycentric = self.config.instrument.barycentric
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
        if flag is False:
            self.logger.error("Header must contain a wavelength/velocity axis.")

        #Get keywords defining wavelength axis
        nwav = header["NAXIS{0}".format(axis)]
        wav0 = header["CRVAL{0}".format(axis)]
        dwav = header["CD{0}_{0}".format(axis)]
        pix0 = header["CRPIX{0}".format(axis)]

        #Calculate and return
        return np.array([wav0 + (i - pix0) * dwav for i in range(nwav)])
