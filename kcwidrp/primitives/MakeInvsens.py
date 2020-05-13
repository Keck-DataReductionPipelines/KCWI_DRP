from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_correct_extin import kcwi_correct_extin
from kcwidrp.core.kcwi_plotting import set_plot_lims, save_plot
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer

from bokeh.plotting import figure, ColumnDataSource
from scipy.signal import find_peaks
from scipy.interpolate import interp1d
from astropy.io import fits as pf
from astropy.nddata import CCDData
from astropy import units as u

import time
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
        stdname = None
        if 'object' in self.action.args.imtype.lower():
            obname = self.action.args.ccddata.header['OBJECT'].lower()
            path = 'data/stds/%s.fits' % obname
            package = __name__.split('.')[0]
            full_path = pkg_resources.resource_filename(package, path)
            if os.path.exists(full_path):
                self.logger.info("Found std file: %s" % full_path)
                stdfile = full_path
                stdname = obname
            else:
                self.logger.info("Not found in data/stds: %s" % full_path)
        else:
            self.logger.warning("Not object type: %s" %
                                self.action.args.imtype)
        self.action.args.stdfile = stdfile
        self.action.args.stdname = stdname
        # check pre condition
        if self.action.args.stdfile is not None:
            # does file already exist?
            ofn = self.action.args.ccddata.header['OFNAME']
            msname = ofn.split('.fits')[0] + '_invsens.fits'
            rdir = self.config.instrument.output_directory
            invsensf = os.path.join(rdir, msname)
            if os.path.exists(invsensf):
                self.logger.info("Master cal already exists: %s" % invsensf)
                return False
            else:
                self.logger.info("Master cal will be generated.")
                return True
        else:
            return False

    def _perform(self):
        self.logger.info("Making inverse sensitivity curve")

        suffix = 'invsens'
        stdname = self.action.args.stdname

        # get size
        sz = self.action.args.ccddata.data.shape
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
        wall0 = self.action.args.ccddata.header['WAVALL0']
        wall1 = self.action.args.ccddata.header['WAVALL1']
        # get DAR padding in y
        pad_y = self.action.args.ccddata.header['DARPADY']
        # get sky subtraction status
        skycor = self.action.args.ccddata.header['SKYCOR']
        # get telescope and atm. correction
        if 'TELESCOP' in self.action.args.ccddata.header:
            tel = self.action.args.ccddata.header['TELESCOP']
        else:
            tel = 'KeckI'
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
        cwv = self.action.args.cwave
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
        bwid = 0.008    # fractional width to mask
        ford = 9        # fit order
        if 'BL' in self.action.args.grating:
            bwid = 0.004
            ford = 7
        elif 'BH' in self.action.args.grating:
            bwid = 0.012
            ford = 9
        bwids = [bl * bwid for bl in blines]
        # fit inverse sensitivity and effective area
        # get initial region of interest
        wl_good = [i for i, v in enumerate(w) if wgoo0 <= v <= wgoo1]
        nwl_good = len(wl_good)
        if nwl_good <= 0:
            self.logger.error("No good wavelengths to fit")
            return self.action.args
        wlm0 = wgoo0
        wlm1 = wgoo1
        # interactively set wavelength limits
        if self.config.instrument.plot_level >= 2:
            yran = [np.min(obsspec), np.max(obsspec)]
            source = ColumnDataSource(data=dict(x=w, y=obsspec))
            done = False
            while not done:
                p = figure(
                    tooltips=[("x", "@x{0,0.0}"), ("y", "@y{0,0.0}")],
                    title=self.action.args.plotlabel + ' %s Obs Spec' % stdname,
                    x_axis_label='Wave (A)',
                    y_axis_label='Intensity (e-)',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                p.line('x', 'y', line_color='black', source=source)
                p.line([wgoo0, wgoo0], yran, line_color='green',
                       legend_label='WAVGOOD')
                p.line([wgoo1, wgoo1], yran, line_color='green')
                p.line([wlm0, wlm0], yran, line_color='blue',
                       legend_label='LIMITS')
                p.line([wlm1, wlm1], yran, line_color='blue')
                p.line([cwv, cwv], yran, line_color='red', legend_label='CWAV')
                set_plot_lims(p, xlim=[wall0, wall1], ylim=yran)
                bokeh_plot(p, self.context.bokeh_session)
                print("WL limits: %.1f - %.1f" % (wlm0, wlm1))
                qstr = input("New? <float> <float>, <cr> - done: ")
                if len(qstr) <= 0:
                    done = True
                else:
                    try:
                        wlm0 = float(qstr.split()[0])
                        wlm1 = float(qstr.split()[1])
                        if wlm1 < wlm0 or wlm0 < wall0 or wlm1 > wall1:
                            wlm0 = wgoo0
                            wlm1 = wgoo1
                            print("range/order error, try again")
                    except (IndexError, ValueError):
                        wlm0 = wgoo0
                        wlm1 = wgoo1
                        print("format error, try again")
            # update region of interest
            wl_good = [i for i, v in enumerate(w) if wlm0 <= v <= wlm1]
            nwl_good = len(wl_good)
        # END: interactively set wavelength limits
        # Now interactively identify lines
        if self.config.instrument.plot_level >= 2:
            yran = [np.min(obsspec), np.max(obsspec)]
            # source = ColumnDataSource(data=dict(x=w, y=obsspec))
            done = False
            while not done:
                p = figure(
                    tooltips=[("x", "@x{0.0}"), ("y", "@y{0.0}")],
                    title=self.action.args.plotlabel + ' %s Obs Spec' % stdname,
                    x_axis_label='Wave (A)',
                    y_axis_label='Intensity (e-)',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                p.line(w, obsspec, line_color='black')
                p.line([wgoo0, wgoo0], yran, line_color='green',
                       legend_label='WAVGOOD')
                p.line([wgoo1, wgoo1], yran, line_color='green')
                p.line([wlm0, wlm0], yran, line_color='blue',
                       legend_label='LIMITS')
                p.line([wlm1, wlm1], yran, line_color='blue')
                p.line([cwv, cwv], yran, line_color='red', legend_label='CWAV')
                for il, bl in enumerate(blines):
                    if wall0 < bl < wall1:
                        p.line([bl, bl], yran, line_color='orange')
                set_plot_lims(p, xlim=[wall0, wall1], ylim=yran)
                bokeh_plot(p, self.context.bokeh_session)
                qstr = input("New lines? <float> [<float>] ... (A), "
                             "<cr> - done: ")
                if len(qstr) <= 0:
                    done = True
                else:
                    for lstr in qstr.split():
                        try:
                            new_line = float(lstr)
                        except ValueError:
                            print("bad line: %s" % lstr)
                            continue
                        if wlm0 < new_line < wlm1:
                            blines.append(new_line)
                            bwids.append(bwid * new_line)
                        else:
                            print("line outside range: %s" % lstr)
        # END: interactively identify lines
        # set up fitting vectors, flux, waves, measure errors
        sf = invsen[wl_good]   # dependent variable
        af = earea[wl_good]    # effective area
        wf = w[wl_good]        # independent variable
        bf = obsspec[wl_good]  # input electrons
        rf = rsflx[wl_good]    # reference flux
        mw = np.ones(nwl_good, dtype=float)     # weights
        use = np.ones(nwl_good, dtype=int)      # toggles for usage
        # loop over Balmer lines
        for il, bl in enumerate(blines):
            roi = [i for i, v in enumerate(wf)
                   if (bl - bwids[il]) <= v <= (bl + bwids[il])]
            nroi = len(roi)
            if nroi > 0:
                use[roi] = 0
                self.logger.info("Masking line at %.1f +- %.4f (A)"
                                 % (bl, bwids[il]))
        # ignore bad points by setting large errors
        mf = []
        ef = []
        ww = []
        used = []
        not_used = []
        for i in range(len(use)):
            if use[i] == 1:
                used.append(i)
            else:
                mw[i] = 1.e-9
                mf.append(sf[i])
                ef.append(100.*af[i]/area)
                ww.append(wf[i])
                not_used.append(i)

        # initial polynomial fit of inverse sensitivity
        wf0 = np.min(wf)
        res = np.polyfit(wf-wf0, sf, deg=ford, w=mw)
        finvsen = np.polyval(res, w-wf0)
        sinvsen = np.polyval(res, wf-wf0)
        calspec = obsspec * finvsen
        scalspec = bf * sinvsen
        # initial polynomial fit of effective area
        res = np.polyfit(wf-wf0, af, ford, w=mw)
        fearea = np.polyval(res, w-wf0)
        # calculate residuals
        resid = 100.0 * (rf - scalspec) / rf
        if len(not_used) > 0:
            rbad = resid[not_used]
        else:
            rbad = None
        rsd_mean = float(np.nanmean(resid[used]))
        rsd_stdv = float(np.nanstd(resid[used]))
        self.logger.info("Calibration residuals = %f +- %f %%" %
                         (rsd_mean, rsd_stdv))
        # plots
        peff = None
        pivs = None
        pcal = None
        prsd = None
        # interactively adjust fit
        if self.config.instrument.plot_level >= 1:
            done = False
            while not done:
                yran = [np.min(100.*af/area), np.max(100.*af/area)]
                peff = figure(
                    title=self.action.args.plotlabel +
                    ' %s Efficiency' % stdname,
                    x_axis_label='Wave (A)',
                    y_axis_label='Effective Efficiency (%)',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                peff.line(wf, 100.*af/area, line_color='black',
                          legend_label='Data')
                peff.line(w, 100.*fearea/area, line_color='green',
                          legend_label='Fit')
                peff.scatter(ww, ef, marker='x', legend_label='Rej')
                peff.line([wlm0, wlm0], yran, line_color='orange')
                peff.line([wlm1, wlm1], yran, line_color='orange')
                set_plot_lims(peff, xlim=[wall0, wall1], ylim=yran)
                bokeh_plot(peff, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(2. * self.config.instrument.plot_pause)

                yran = [np.min(sf), np.max(sf)]
                pivs = figure(
                    title=self.action.args.plotlabel +
                    ' %s Inverse sensitivity' % stdname,
                    x_axis_label='Wave (A)',
                    y_axis_label='Invserse Sensitivity (Flux/e-/s)',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                pivs.line(wf, sf, line_color='black', legend_label='Data')
                pivs.line(w, finvsen, line_color='green', legend_label='Fit')
                pivs.scatter(ww, mf, marker='x', legend_label='Rej')
                pivs.line([wlm0, wlm0], yran, line_color='orange')
                pivs.line([wlm1, wlm1], yran, line_color='orange')
                set_plot_lims(pivs, xlim=[wall0, wall1], ylim=yran)
                bokeh_plot(pivs, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(2. * self.config.instrument.plot_pause)

                yran = [np.min(calspec), np.max(calspec)]
                pcal = figure(
                    title=self.action.args.plotlabel +
                    ' %s Calibrated' % stdname,
                    x_axis_label='Wave (A)',
                    y_axis_label='Flux (ergs/s/cm^s/A)',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                pcal.line(w, calspec, line_color='black', legend_label='Obs')
                pcal.line(w, rsflx, line_color='red', legend_label='Ref')
                pcal.line([wlm0, wlm0], yran, line_color='orange')
                pcal.line([wlm1, wlm1], yran, line_color='orange')
                set_plot_lims(pcal, xlim=[wall0, wall1], ylim=yran)
                bokeh_plot(pcal, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(2. * self.config.instrument.plot_pause)

                yran = [np.min(resid), np.max(resid)]
                prsd = figure(
                    title=self.action.args.plotlabel +
                    ' %s Residuals = %.1f +- %.1f (%%)' % (stdname, rsd_mean,
                                                           rsd_stdv),
                    x_axis_label='Wave (A)',
                    y_axis_label='Ref - Cal / Ref (%)',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                prsd.line(wf, resid, line_color='black',
                          legend_label='Ref - Cal (%)')
                if len(not_used) > 0:
                    prsd.scatter(ww, rbad, marker='x', legend_label='Rej')
                prsd.line([wlm0, wlm0], yran, line_color='orange')
                prsd.line([wlm1, wlm1], yran, line_color='orange')
                prsd.line([wall0, wall1], [0, 0], line_color='blue')
                prsd.line([wall0, wall1],
                          [rsd_mean+rsd_stdv, rsd_mean+rsd_stdv],
                          line_color='green')
                prsd.line([wall0, wall1],
                          [rsd_mean-rsd_stdv, rsd_mean-rsd_stdv],
                          line_color='green')
                set_plot_lims(prsd, xlim=[wall0, wall1], ylim=yran)
                bokeh_plot(prsd, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    qstr = input("Current fit order = %d, "
                                 "New fit order? <int>, <cr> - done: " % ford)
                    if len(qstr) <= 0:
                        done = True
                    else:
                        try:
                            ford = int(qstr)
                            # update fit of inverse sensitivity
                            res = np.polyfit(wf - wf0, sf, deg=ford, w=mw)
                            finvsen = np.polyval(res, w - wf0)
                            sinvsen = np.polyval(res, wf - wf0)
                            calspec = obsspec * finvsen
                            scalspec = bf * sinvsen
                            # update polynomial fit of effective area
                            res = np.polyfit(wf - wf0, af, ford, w=mw)
                            fearea = np.polyval(res, w - wf0)
                            # re-calculate residuals
                            resid = 100.0 * (rf - scalspec) / rf
                            if len(not_used) > 0:
                                rbad = resid[not_used]
                            else:
                                rbad = None
                            rsd_mean = float(np.nanmean(resid[used]))
                            rsd_stdv = float(np.nanstd(resid[used]))
                            self.logger.info(
                                "Calibration residuals = %f +- %f %%" %
                                (rsd_mean, rsd_stdv))
                        except ValueError:
                            print("Bad fit order, try again")
                else:
                    done = True
                    time.sleep(2. * self.config.instrument.plot_pause)
            # output plots
            imstr = "%05d" % self.action.args.ccddata.header['FRAMENO']
            cwv = "%d" % int(self.action.args.cwave)
            grat = self.action.args.grating.strip()
            ifu = self.action.args.ifuname.strip()
            pfname = os.path.join(self.config.instrument.output_directory,
                                  stdname + '_' + grat + '_' + cwv + '_' +
                                  ifu + '_' + imstr)
            # Save plots
            save_plot(peff, filename=pfname + '_eff.png')      # Efficiency
            save_plot(pivs, filename=pfname + '_invsens.png')  # Inv. Sens.
            save_plot(prsd, filename=pfname + '_resid.png')    # Residuals
            save_plot(pcal, filename=pfname + '_cal.png')      # Calibrated

        log_string = MakeInvsens.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        # write out effective inverse sensitivity

        # update inverse sensitivity header
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
        eaname = ofn.split('.fits')[0] + '_ea.fits'

        # set units
        invsens_u = u.erg / (u.angstrom * u.cm ** 2 * u.s * u.electron)
        # output inverse sensitivity
        out_invsens = CCDData(np.asarray([invsen, finvsen, obsspec]),
                              meta=hdr, unit=invsens_u)
        kcwi_fits_writer(out_invsens, output_file=invsname,
                         output_dir=self.config.instrument.output_directory)
        self.context.proctab.update_proctab(frame=out_invsens, suffix=suffix,
                                            newtype='INVSENS')
        # output effective area
        ea_u = u.cm ** 2 / u.angstrom
        out_ea = CCDData(np.asarray([earea, fearea]), meta=hdr, unit=ea_u)
        kcwi_fits_writer(out_ea, output_file=eaname,
                         output_dir=self.config.instrument.output_directory)
        self.context.proctab.update_proctab(frame=out_ea, suffix='ea',
                                            newtype='EAREA')

        return self.action.args
    # END: class MakeInvsens()
