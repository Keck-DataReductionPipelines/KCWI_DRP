import os

import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import PchipInterpolator
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation

from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer, \
                                                    kcwi_fits_reader, \
                                                    strip_fname

class WavelengthCorrections(BasePrimitive):
    """
    Perform wavelength corrections

    Based on the type specified in the kcwidrp/configs/kcwi.cfg file in the
    ``radial_velocity_correction`` parameter, perform the requested correction.
    The current options are "heliocentric", "barycentric", and "none".  The
    default is "heliocentric".

    Also, if the ``air_to_vacuum`` parameter is set to ``True`` (default),
    apply air-to-vacuum correction.

    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if we can correct wavelengths based on the processing table
        """
        self.logger.info("Checking precondition for WavelengthCorrections")
        precondition = False
        suffix = 'icube'
        obj = self.locate_object_file(suffix)
        if obj:
            precondition = True
        else:
            self.logger.error("Precondition for WavelengthCorrections failed!")
        return precondition

    def _perform(self):

        # Determine which radial velocity correction to make
        correction_mode = self.config.instrument.radial_velocity_correction
        options = ["none", "barycentric", "heliocentric"]

        # If the config file has an invalid option, return
        if not bool([el for el in options if el in correction_mode]):
            self.logger.error('Bad config option for radial_velocity_correction\
                , options are ["none", "heliocentric", "barycentric"]')
            return self.action.args

        suffix = 'icube'  # Can be ammended to handle ocube files
        obj = self.locate_object_file(suffix)

        if "none" in correction_mode:
            perform_helio = False
            self.logger.info("Skipping radial velocity correction")
            obj.header['VCORR'] = (0.0, 'km/s')
            obj.header['VCORRTYP'] = (correction_mode, 'Vcorr type')
        else:
            perform_helio = True

        perform_a2v = self.config.instrument.air_to_vacuum

        # does helio corr, then a2v, then resamples
        obj = self.helio_a2v_conversion(obj, perform_helio, perform_a2v, correction_mode)

        log_string = WavelengthCorrections.__module__
        obj.header['HISTORY'] = log_string

        kcwi_fits_writer(obj,
                         table=self.action.args.table,
                         output_file=self.action.args.name,
                         output_dir=self.config.instrument.output_directory,
                         suffix=f'{suffix}w')
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix=f'{suffix}w',
                                            filename=self.action.args.name)
        self.context.proctab.write_proctab(tfil=self.config.instrument.procfile)

        # Unsure here: Is this right? it seems to make DAR happy
        self.action.args.ccddata = obj

        return self.action.args

    def helio_a2v_conversion(self, obj, perform_helio, perform_a2v, correction_mode, mask=False, vcorr=None, resample=True):
        """
        Apply heliocentric and/or air to vacuum conversion to wavelengths in a cube.

        Adapted from https://github.com/dbosul/cwitools.git

        Args:
            obj (astropy HDU / HDUList): Input HDU/HDUList with 3D data.
            perform_helio (bool): Set to apply heliocentric conversion.
            perform_a2v (bool): Set to apply air to vac conversion.
            correction_mode (str): "none", "barycentric", or "heliocentric"
            vcorr (float): Use a different correction velocity.
            mask (bool): Set if the cube is a mask cube.

        Returns:
            HDU / HDUList: Trimmed FITS object with updated header.
            float: (if vcorr is True) Correction velocity in km/s.
            Return type matches type of fits_in argument.
        """

        cube = np.nan_to_num(obj.data, nan=0, posinf=0, neginf=0)

        wav_old = self.get_wav_axis(obj.header)

        wav_air = wav_old * u.AA

        if perform_helio:

            self.logger.info(f"Performing {correction_mode} correction")

            v_old = 0.
            if 'VCORR' in obj.header:
                v_old = obj.header['VCORR']
                self.logger.info("Rolling back the existing correction with:")
                self.logger.info("Vcorr = %.2f km/s." % v_old)

            # get correction velocity
            vcorr = self.heliocentric(obj, correction_mode, vcorr)
            v_tot = vcorr - v_old
            if not resample:
                obj.header['CRVAL3'] = str(obj.header['CRVAL3'] * (1 + v_tot / 2.99792458e5))
                obj.header['CD3_3'] = str(obj.header['CD3_3'] * (1 + v_tot / 2.99792458e5))

            # update headers
            obj.header['VCORR'] = (vcorr, 'km/s')
            obj.header['VCORRTYP'] = (correction_mode, 'Vcorr type')

            # calculate new wavelengths with helio corrections
            wav_new = wav_old * (1 + v_tot / 2.99792458e5)
            cwave_hel = wav_old[int(wav_old.shape[0] / 2)] * (1 + v_tot / 2.99792458e5)
            self.logger.info("Vcorr for CWAVE (%.3f) gives %.3f" % 
                             (wav_old[int(wav_old.shape[0] / 2)], cwave_hel))

            # if also doing air to vac, propogate the helio correction and convert units
            if perform_a2v:
                wav_air = wav_new * u.AA

        if perform_a2v:

            self.logger.info("Performing Air to Vacuum Conversion")

            # get vacuum wavelengths
            wav_new = self.air_to_vac(wav_air)
            cwave_air = wav_air[int(wav_air.shape[0] / 2)]
            cwave_vac = wav_new[int(wav_new.shape[0] / 2)]

            # update headers
            self.logger.info("New Air to Vacuum for (%.3f) gives %.3f" % (cwave_air.value, cwave_vac.value))
            obj.header['CTYPE3'] = ('WAVE', 'Vacuum Wavelengths')

        # resample to uniform grid
        if resample:
            self.logger.info("Resampling to uniform grid")
            cube_new = np.zeros_like(cube)
            for i in range(cube.shape[2]):
                for j in range(cube.shape[1]):

                    spec0 = cube[:, j, i]
                    if not mask:

                        f_cubic = interp1d(wav_new, spec0, kind='cubic',
                                         fill_value='extrapolate')

                        spec_new = f_cubic(wav_old)

                        f_cubic = PchipInterpolator(wav_new, spec0)
                        spec_new = f_cubic(wav_old)
                    else:
                        f_pre = interp1d(wav_new, spec0, kind='previous',
                                            bounds_error=False, fill_value=128)
                        spec_pre = f_pre(wav_old)

                        f_nex = interp1d(wav_new, spec0, kind='next',
                                            bounds_error=False, fill_value=128)
                        spec_nex = f_nex(wav_old)

                        spec_new = np.zeros_like(spec0)
                        for k in range(spec0.shape[0]):
                            spec_new[k] = max(spec_pre[k], spec_nex[k])
                    cube_new[:, j, i] = spec_new

            obj.data = cube_new
        return obj

    def air_to_vac(self, wave):
        """
        Convert air-based wavelengths to vacuum wavelengths

        Adapted from wave.py in: https://github.com/pypeit/PypeIt/.
        Formulas from https://ui.adsabs.harvard.edu/abs/1996ApOpt..35.1566C/.

        Args:
            wave (ndarray): Wavelengths.

        Returns:
            ndarray: Wavelength array corrected to vacuum wavelengths.

        """
        # Convert to AA
        wave = wave.to(u.AA)
        wavelength = wave.value

        # wave number (inverse wavelength in micrometers)
        sigma = (1.e4/wavelength)
        sigma_sq = sigma**2

        # humidity: convert from % to decimal
        h = self.action.args.dome_hum / 100
        # temperature: convert from [C] to [K]
        T = self.action.args.dome_temp + 273.15
        # pressure: convert from mbar to Pa
        p = self.action.args.pres * 100
        # ppm of CO2
        x_c = 450

        # Constants for equation (1) from Appendix A
        k0 = 238.0185
        k1 = 5792105
        k2 = 57.362
        k3 = 167917

        # Index of refraction for standard air. From equation (1)
        n_as = ((k1 / (k0 - sigma_sq) + k3 / (k2 - sigma_sq)) / 10**8) + 1

        # Index of refraction for standard air with x_c ppm CO2. From equation (2)
        n_axs = ((n_as - 1) * (1 + 0.534 * 10**-6 * (x_c - 450))) + 1

        # Constants for equation (3) from Appendix A
        cf = 1.022 # correction factor
        w0 = 295.235
        w1 = 2.6422
        w2 = -0.032380
        w3 = 0.004028

        # Index of refraction for standard water vapor. From equation (3)
        n_ws = ((cf * (w0 + (w1 * sigma**2) + (w2 * sigma**4) + (w3 * sigma**6))) / 10e8) + 1

        def get_molar_frac_water(T, p, h):
            """
            Gets the molar fraction of water vapor in moist air [kg/m^3].
            Part of equation (4) in reference paper.

            Args:
                T (float): Temperature in [K] Kelvin.
                p (float): Pressure in [Pa] Pascal.
                h (float): Percent relative humidity.

            Returns:
                float: Molar fraction of water vapor in moist air.
            """

            # temperature [C]
            t = T - 273.15

            # constants for f from Appendix A
            alpha = 1.00062
            beta = 3.14 * 10**-8
            gamma  = 5.6 * 10**-7

            # enchancement factor. part of Equation (4)
            f = alpha + beta * p + gamma * t**2

            # constants for svp from Appendix A
            A = 1.2378847 * 10**-5
            B = -1.9121316 * 10**-2
            C = 33.93711047
            D = -6.3431645 * 10**3 

            # saturation vapor pressure (svp) of water vapor in air at temperature, T [K]. Part of Equation (4)
            svp = np.exp((A * T**2 + B * T + C + D/T))

            # molar fraction of water vapor in moist air. Part of Equation (4)
            x_w = f * h * (svp/p)

            return x_w 
        
        def get_compress(T, p, x_w):
            """
            Gets the compressability of moist air.
            From Equation (12) in Appendix A.

            Args:
                T (float): Temperature in [K] Kelvin.
                p (float): Pressure in [Pa] Pascals.
                x_w (float): Molar fraction of water vapor.

            Returns:
                float: Compressability of moist air.
            """

            # t = temperature [C]
            t = T - 273.15

            # constants for Z from Appendix A
            a0 = 1.58123 * 10**-6
            a1 = -2.9331 * 10**-8
            a2 = 1.1043 * 10**-10
            b0 = 5.707 * 10**-6
            b1 = -2.051 * 10**-8
            c0 = 1.9898 * 10**-4
            c1 = -2.376 * 10**-6
            d = 1.83 * 10**-11
            e = -0.765 * 10**-8

            # Compresibility of moist air. From equation (12) in Appendix A
            Z = 1 - (p/T) * (a0 + a1*t + a2 * t**2 + (b0 + b1*t) * x_w + (c0 + c1*t) * x_w**2) + (p/T)**2 * (d + e * x_w**2)
            return Z

        def get_density(T, p, h, x_c, x_w=None):
            """
            Gets the density of moist air in [kg/m^3].
            From equation (4) in reference paper.

            Args:
                T (float): Temperature in [K] Kelvin.
                p (float): Pressure in [Pa] Pascal.
                h (float): Percent relative humidity.
                x_c (int): Ppm of CO2.
                x_w (float, optional): Set a default water vaport molar fraction. Defaults to None.

            Returns:
                float: Density of moist air.
            """

            # molar mass of water vapor [kg/mol]
            M_w = 0.018015
            # molar mass of dry air with x_c ppm of CO2 [kg/mol]
            M_a = 10**-3 * (28.9635 + 12.011 * 10**-6 * (x_c - 400))
            # gas constant [J/(mol*K)]
            R = 8.314510

            # molar fraction of water vapor in moist air
            if x_w is None:
                x_w = get_molar_frac_water(T, p, h)

            # compresibility of moist air
            Z = get_compress(T, p, x_w)

            # density of moist air. From equation (4)
            rho = (p * M_a / Z * R * T) * (1 - x_w * (1 - M_w / M_a))

            return rho

        def get_rind(T, p, h, x_c):
            """
            Gets the refractive index of moist air.
            From Equation (5) in the reference paper.
            Uses equations for dry air (found in step 8 of Appendix B) 
            and moist air (found in step 9 of Appendix B) components.

            Args:
                T (float): Temperature in [K] Kelvin.
                p (float): Pressure in [Pa] Pascal.
                h (float): Percent relative humidity.
                x_c (int): Ppm of CO2.

            Returns:
                float: Refractive index of moist air.
            """

            # molar mass of water vapor [kg/mol]
            M_w = 0.018015
            # molar mass of dry air with x_c ppm of CO2 [kg/mol]
            M_a = 10**-3 * (28.9635 + 12.011 * 10**-6 * (x_c - 400))
            # gas constant [J/(mol*K)]
            R = 8.314510

            # density of standard dry air (15 deg C, 101325 Pa, x_w = 0, 450 ppm CO2)
            rho_axs = get_density(288.15, 101315, 0, 450, 0)
            # density of standard water vapor (20 deg C, 1333 Pa, x_w = 1)
            rho_ws = get_density(293.15, 1333, 1, 0, 1)
            # molar fraction of water vapor in moist air
            x_w = get_molar_frac_water(T, p, h)
            # compressability of moist air
            Z = get_compress(T, p, x_w)

            # density of dry air component of moist air with actual conditions. From step 8 in Appendix B
            rho_a = p * M_a * (1 - x_w) / Z * R * T 
            # density of water vapor component of moist air with actual conditons. From step 9. in Appendix B
            rho_w = p * M_w * x_w / Z * R * T

            # refractive index of moist air. From Equation (5)
            n_prop = 1 + ((rho_a/rho_axs) * (n_axs - 1) + (rho_w/rho_ws) * (n_ws - 1))

            return n_prop

        # refractive index
        n_prop = get_rind(T, p, h, x_c)

        new_rind = n_prop[int(n_prop.shape[0] / 2)]
        rwav = wavelength[int(wavelength.shape[0] / 2)]
        self.logger.info("Refractive index = %.10f at %.3f Ang" % (new_rind, rwav))
        
        # only modify above 2000A
        n_prop = n_prop*(wavelength >= 2000.) + 1.*(wavelength < 2000.)

        # Convert
        new_wavelength = wavelength*n_prop

        # Units
        new_wave = new_wavelength*u.AA
        new_wave.to(wave.unit)

        return new_wave

    def heliocentric(self, obj, correction_mode, vcorr=None):
        """
        Apply heliocentric correction to old wavelengths.

        Adapted from https://github.com/dbosul/cwitools.git

        Args:
            obj (astropy HDU / HDUList): Input HDU/HDUList with 3D data.
            correction_mode (str): "none", "barycentric", or "heliocentric".
            vcorr (float): Use a different correction velocity.

        Returns:
            float: Correction velocity in km/s.
        """

        # check if barycentric correction
        barycentric = ("barycentric" in correction_mode)

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
        self.logger.info("Vcorr = %.2f km/s." % vcorr)

        return vcorr

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

        # Select the appropriate axis.
        naxis = header['NAXIS']
        axis = None
        for i in range(naxis):
            # Keyword entry
            card = "CTYPE{0}".format(i + 1)
            if card not in header:
                self.logger.warning.error(
                    "Header must contain 'CTYPE' keywords.")

            # Possible wave types.
            if header[card] in ['AWAV', 'WAVE', 'VELO']:
                axis = i + 1
                break

        # No wavelength axis
        if axis is None:
            self.logger.error("Header must contain a wavelength/velocity axis.")
            retval = None

        else:
            # Get keywords defining wavelength axis
            nwav = header["NAXIS{0}".format(axis)]
            wav0 = header["CRVAL{0}".format(axis)]
            dwav = header["CD{0}_{0}".format(axis)]
            pix0 = header["CRPIX{0}".format(axis)]

            # Calculate and return
            retval = np.array([wav0 + (i - pix0) * dwav for i in range(nwav)])

        return retval

    def locate_object_file(self, suffix):
        """
        Return FITS HDU list if current file with requested suffix can be found.

        Will return ``None`` if file cannot be found.

        Args:
            suffix (str): The file suffix that you want to operate on.

        Returns:
            FITS HDU List from kcwi_fits_reader routine, or ``None`` if file
            with suffix not found.

        """
        ofn = self.action.args.name
        objfn = strip_fname(ofn) + f'_{suffix}.fits'
        full_path = os.path.join(
            self.config.instrument.cwd,
            self.config.instrument.output_directory, objfn)
        if os.path.exists(full_path):
            return kcwi_fits_reader(full_path)[0]
        else:
            self.logger.error(f'Unable to read file {objfn}')
            return None
