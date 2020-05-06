from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_correct_extin import kcwi_correct_extin
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer

from bokeh.plotting import figure
from scipy.signal import find_peaks
from scipy.interpolate import interp1d
from astropy.io import fits as pf

import os
import pkg_resources
import numpy as np


class MakeInvsens(BasePrimitive):

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """Checks if object is a standard star"""
        self.logger.info("Checking precondition for MakeInvsens")

        stdfile = None
        if 'object' in self.action.args.imtype.lower():
            obname = self.action.args.ccddata.header['OBJECT'].lower()
            path = 'data/stds/%s.fits' % obname
            package = __name__.split('.')[0]
            full_path = pkg_resources.resource_filename(package, path)
            if os.path.exists(full_path):
                self.logger.info("Found std file: %s" % full_path)
                stdfile = full_path
            else:
                self.logger.info("Not found in data/stds: %s" % full_path)
        else:
            self.logger.warning("Not object type: %s" %
                                self.action.args.imtype)
        self.action.args.stdfile = stdfile
        if self.action.args.stdfile is not None:
            return True
        else:
            return False

    def _perform(self):
        self.logger.info("Making inverse sensitivity curve")

        suffix = 'invsens'

        # get size
        sz = self.action.args.ccddata.data.shape
        print(sz)
        # default pixel ranges
        z = np.arange(sz[0])
        z0 = 175
        z1 = sz[0] - 175
        # get exposure time
        expt = self.action.args.ccddata.header['XPOSURE']
        if expt == 0.:
            self.logger.warning("No exposure time found, setting to 1s")
            expt = 1.
        else:
            self.logger.info("Using exposure time of %.1f" % expt)
        # get wavelength scale
        w0 = self.action.args.ccddata.header['CRVAL3']
        dw = self.action.args.ccddata.header['CD3_3']
        crpixw = self.action.args.ccddata.header['CRPIX3']
        # get all good wavelength range
        wgoo0 = self.action.args.ccddata.header['WAVGOOD0']
        if wgoo0 < 3650:
            wgoo0 = 3650.
        wgoo1 = self.action.args.ccddata.header['WAVGOOD1']
        # get all inclusive wavelength range
        # wall0 = self.action.args.ccddata.header['WAVALL0']
        # wall1 = self.action.args.ccddata.header['WAVALL1']
        # get DAR padding in y
        pad_y = self.action.args.ccddata.header['DARPADY']
        # get sky subtraction status
        skycor = self.action.args.ccddata.header['SKYCOR']
        # get telescope and atm. correction
        tel = self.action.args.ccddata.header['TELESCOP']
        if 'Keck' in tel:
            area = 760000.0
        else:
            area = -1.0
        # tlab = tel
        # compute good y pixel ranges
        if w0 > 0. and dw > 0. and wgoo0 > 0. and wgoo1 > 0.:
            z0 = int((wgoo0 - w0) / dw) + 10
            z1 = int((wgoo1 - w0) / dw) - 10
        # gz = [i for i, v in enumerate(z) if z0 <= v <= z1]
        # wavelength scale
        w = w0 + z * dw
        # good spatial range
        gy0 = pad_y if pad_y > 1 else 1
        gy1 = sz[1] - (pad_y if pad_y > 2 else 2)
        # log results
        self.logger.info("Invsen. Pars: Y0, Y1, Z0, Z1, Wav0, Wav1: "
                         "%d, %d, %d, %d, %.3f, %.3f" %
                         (gy0, gy1, z0, z1, w[z0], w[z1]))
        # central wavelength
        # cwv = self.action.args.cwave
        # find standard in slices
        # sum over wavelength
        tot = np.sum(self.action.args.ccddata.data[z0:z1, gy0:gy1, :], 0)
        # yy = np.arange(gy1-gy0) + gy0
        mxsl = -1.
        mxsg = 0.
        # for each slice
        for i in range(sz[2]):
            tstd = float(np.nanstd(tot[:, i]))
            if tstd > mxsg:
                mxsg = tstd
                mxsl = i

        # relevant slices
        sl0 = (mxsl - 3) if mxsl >= 3 else 0
        sl1 = (mxsl + 3) if (mxsl + 3) <= sz[2]-1 else sz[2]-1
        # get y position of std
        cy, _ = find_peaks(tot[:, mxsl], height=np.nanmean(tot[:, mxsl]))
        cy = int(cy[0]) + gy0
        # log results
        self.logger.info("Std lices: max, sl0, sl1, spatial cntrd: "
                         "%d, %d, %d, %.2f" % (mxsl, sl0, sl1, cy))
        # get dwave spectrum
        dwspec = np.zeros(sz[0]) + dw
        # copy of input cube
        scub = self.action.args.ccddata.data.copy()
        # sky window width in pixels
        skywin = int(self.config.instrument.psfwid / self.action.args.xbinsize)
        # do sky subtraction, if needed
        if not skycor:
            self.logger.warning("Sky should have been subtraced already")
        # apply extinction correction
        ucub = scub.copy()  # uncorrected cube
        kcwi_correct_extin(scub, self.action.args.ccddata.header,
                           logger=self.logger)

        # get slice spectra y limits
        sy0 = (cy - skywin) if (cy - skywin) > 0 else 0
        sy1 = (cy + skywin) if (cy + skywin) < (sz[1]-1) else (sz[1]-1)
        # sum over y range
        slspec = np.sum(scub[:, sy0:sy1, :], 1)
        ulspec = np.sum(ucub[:, sy0:sy1, :], 1)
        # sum over slices
        obsspec = np.sum(slspec[:, sl0:sl1], 1)
        ubsspec = np.sum(ulspec[:, sl0:sl1], 1)
        # convert to e-/second
        obsspec /= expt
        ubsspec /= expt
        # read in standard star spectrum
        hdul = pf.open(self.action.args.stdfile)
        swl = hdul[1].data['WAVELENGTH']
        sflx = hdul[1].data['FLUX']
        sfw = hdul[1].data['FWHM']
        hdul.close()
        # get region of interest
        sroi = [i for i, v in enumerate(swl) if w[0] <= v <= w[-1]]
        nsroi = len(sroi)
        if nsroi <= 0:
            self.logger.error("No standard wavelengths in common")
            return self.action.args
        # expand range after checking for edges
        if sroi[0] > 0:
            sroi.insert(0, sroi[0]-1)
            nsroi += 1
        if sroi[-1] < (len(swl)-1):
            sroi.append(sroi[-1]+1)
            nsroi += 1
        # very sparsely sampled w.r.t. object
        if nsroi <= 1:
            self.logger.error("Not enough standard points")
            return self.action.args
        self.logger.info("Number of standard point = %d" % nsroi)
        # how many points?
        # do_line = nsroi > 20

        swl = swl[sroi]
        sflx = sflx[sroi]
        sfw = sfw[sroi]
        fwhm = np.max(sfw)
        self.logger.info("Reference spectrum FWHM used = %.1f (A)" % fwhm)
        # resample standard onto our wavelength grid
        rsint = interp1d(swl, sflx, kind='cubic')
        rsflx = rsint(w)
        # get effective inverse sensitivity
        invsen = rsflx / obsspec
        # convert to photons/s/cm^2/(wl bin = dw)
        rspho = 5.03411250e7 * rsflx * w * dw
        # get effective area
        earea = ubsspec / rspho
        # correct to native bins
        earea *= dw/dwspec
        # Balmer lines
        blines = [6563., 4861., 4341., 4102., 3970., 3889., 3835.]
        # default values (for BM)
        bwid = 0.017    # fractional width to mask
        ford = 9        # fit order
        sigf = 3.0      # rejection sigma
        if 'BL' in self.action.args.grating:
            bwid = 0.008
            ford = 7
            sigf = 4.0
        elif 'BH' in self.action.args.grating:
            bwid = 0.008
            ford = 9
            sigf = 3.0
        # fit inverse sensitivity and effective area
        t = [i for i, v in enumerate(w) if wgoo0 <= v <= wgoo1]
        nt = len(t)
        if nt > 0:
            # TODO: put in interactive setting of these values
            wlm0 = wgoo0
            wlm1 = wgoo1
        else:
            self.logger.error("No good wavelengths to fit")
            return self.action.args
        # set up fitting vectors, flux, waves, measure errors
        sf = invsen[t]
        wf = w[t]
        mw = np.ones(sf.shape, dtype=float)
        use = np.ones(nt)
        # loop over Balmer lines
        for bl in blines:
            roi = [i for i, v in enumerate(wf)
                   if (bl - bl*bwid) <= v < (bl + bl*bwid)]
            nroi = len(roi)
            if nroi > 0:
                use[roi] = 0
                self.logger.info("Masking Balmer line at %.1f +- %.4f (A)"
                                 % (bl, bwid))
        # ignore bad points by setting large errors
        mf = []
        ww = []

        for i in range(len(use)):
            if use[i] != 1:
                mw[i] = 1.e-9
                mf.append(sf[i])
                ww.append(wf[i])
        # initial polynomial fit of inverse sensitivity
        wf0 = np.min(wf)
        res = np.polyfit(wf-wf0, sf, deg=ford, w=mw)
        finvsen = np.polyval(res, w-wf0)
        calspec = obsspec * finvsen

        # TODO: put in interactive tweaking of invsens fitting

        yran = [np.min(sf), np.max(sf)]
        p = figure(
            title=self.action.args.plotlabel + ' Invsens',
            x_axis_label='Wave (A)',
            y_axis_label='Invserse Sensitivity (Flux/e-/s)',
            plot_width=self.config.instrument.plot_width,
            plot_height=self.config.instrument.plot_height)
        p.line(wf, sf, line_color='black', legend_label='Data')
        p.line(w, finvsen, line_color='green', legend_label='Fit')
        p.scatter(ww, mf, marker='x')
        p.line([wgoo0, wgoo0], yran, line_color='orange')
        p.line([wgoo1, wgoo1], yran, line_color='orange')
        p.y_range.start = yran[0]
        p.y_range.end = yran[1]
        bokeh_plot(p, self.context.bokeh_session)
        input("Next? <cr>: ")

        yran = [np.min(calspec), np.max(calspec)]
        p = figure(
            title=self.action.args.plotlabel + ' Calib',
            x_axis_label='Wave (A)',
            y_axis_label='Flux (ergs/s/cm^s/A)',
            plot_width=self.config.instrument.plot_width,
            plot_height=self.config.instrument.plot_height)
        p.line(w, calspec, line_color='black', legend_label='Data')
        p.line(w, rsflx, line_color='red', legend_label='Ref')
        p.line([wgoo0, wgoo0], yran, line_color='orange')
        p.line([wgoo1, wgoo1], yran, line_color='orange')
        p.y_range.start = yran[0]
        p.y_range.end = yran[1]
        bokeh_plot(p, self.context.bokeh_session)
        input("Next? <cr>: ")
        log_string = MakeInvsens.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        # write out effictive inverse sensitivity

        # update invsens header
        hdr = self.action.args.ccddata.header.copy()
        hdr['HISTORY'] = log_string
        hdr['INVSENS'] = (True, 'effective inv. sens. spectrum?')
        hdr['INVSW0'] = (w[z0], 'low wavelength for eff. inv. sens.')
        hdr['INVSW1'] = (w[z1], 'high wavelength for eff. inv. sens.')
        hdr['INVSZ0'] = (z0, 'low wave pixel for eff. inv. sens.')
        hdr['INVSZ1'] = (z1, 'high wave pixel for eff. inv. sens.')
        hdr['INVSY0'] = (gy0, 'low spatial pixel for eff. inv. sens.')
        hdr['INVSY1'] = (gy1, 'high spatial pixel for eff. inv. sens.')
        hdr['INVSLMX'] = (mxsl, 'brightest std star slice')
        hdr['INVSL0'] = (sl0, 'lowest std star slice summed')
        hdr['INVSL1'] = (sl1, 'highest std star slice summed')
        hdr['INVSLY'] = (cy, 'spatial pixel position of std within slice')
        # this apparently doesn't work (ends up being electron)'
        hdr['BUNIT'] = ('erg/cm^2/A/e-', 'brightness units')
        hdr['EXPTIME'] = (1., 'effective exposure time (seconds)')
        hdr['XPOSURE'] = (1., 'effective exposure time (seconds)')

        # remove old WCS
        del hdr['RADESYS']
        del hdr['EQUINOX']
        del hdr['LONPOLE']
        del hdr['LATPOLE']
        del hdr['NAXIS2']
        if 'NAXIS3' in hdr:
            del hdr['NAXIS3']
        del hdr['CTYPE1']
        del hdr['CTYPE2']
        del hdr['CTYPE3']
        del hdr['CUNIT1']
        del hdr['CUNIT2']
        del hdr['CUNIT3']
        del hdr['CNAME1']
        del hdr['CNAME2']
        del hdr['CNAME3']
        del hdr['CRVAL1']
        del hdr['CRVAL2']
        del hdr['CRVAL3']
        del hdr['CRPIX1']
        del hdr['CRPIX2']
        del hdr['CRPIX3']
        del hdr['CD1_1']
        del hdr['CD1_2']
        del hdr['CD2_1']
        del hdr['CD2_2']
        del hdr['CD3_3']

        # set wavelength axis WCS values
        hdr['WCSDIM'] = 1
        hdr['CTYPE1'] = ('AWAV', 'Air Wavelengths')
        hdr['CUNIT1'] = ('Angstrom', 'Wavelength units')
        hdr['CNAME1'] = ('KCWI INVSENS Wavelength', 'Wavelength name')
        hdr['CRVAL1'] = (w0, 'Wavelength zeropoint')
        hdr['CRPIX1'] = (crpixw, 'Wavelength reference pixel')
        hdr['CDELT1'] = (dw, 'Wavelength Angstroms per pixel')

        ofn = self.action.args.ccddata.header['OFNAME']
        invsname = ofn.split('.fits')[0] + '_' + suffix + '.fits'

        st_hdr = self.action.args.ccddata.header.copy()
        st_img = self.action.args.ccddata.data

        self.action.args.ccddata.header = hdr
        self.action.args.ccddata.data = finvsen

        # output invsense
        kcwi_fits_writer(self.action.args.ccddata, output_file=invsname)
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix=suffix,
                                            newtype='INVSENS')
        self.context.proctab.write_proctab()

        # restore original image
        self.action.args.ccddata.data = st_img
        self.action.args.ccddata.header = st_hdr

        return self.action.args
    # END: class MakeInvsens()
