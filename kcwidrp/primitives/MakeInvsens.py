from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_correct_extin import kcwi_correct_extin
from kcwidrp.core.kcwi_get_std import kcwi_get_std
from kcwidrp.core.kcwi_plotting import set_plot_lims, save_plot
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer, \
    kcwi_fits_reader, strip_fname

from bokeh.plotting import figure, ColumnDataSource
from scipy.signal import find_peaks
from scipy.interpolate import interp1d
from astropy.io import fits as pf
from astropy.nddata import CCDData
from astropy import units as u

import time
import os
import numpy as np


class MakeInvsens(BasePrimitive):
    """
    Generate inverse sensitivity curve from a standard star observation.

    Uses object name to determine if the observation is a standard star.  Then
    checks that the image has been processed through DAR correction.  Reads
    in the reference spectrum and compares it with the observed spectrum to
    generate the inverse sensitivity curve that can be used to flux calibrate
    science observations.

    Outputs FITS inverse sensitivity spectrum along with diagnostic plots that
    include residuals between reference spectrum and calibrated observed
    spectrum and effective area and efficiency plots as a function of
    wavelength.

    """
    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """Checks if object is a standard star"""
        self.logger.info("Checking precondition for MakeInvsens")

        stdfile = None
        stdname = None
        if 'object' in self.action.args.imtype.lower():
            self.logger.info("Checking OBJECT keyword")
            stdfile, stdname = kcwi_get_std(
                self.action.args.ccddata.header['OBJECT'], self.logger)
            if not stdfile:
                self.logger.info("Checking TARGNAME keyword")
                stdfile, stdname = kcwi_get_std(
                    self.action.args.ccddata.header['TARGNAME'], self.logger)
        else:
            self.logger.warning("Not object type: %s" %
                                self.action.args.imtype)
        self.action.args.stdfile = stdfile
        self.action.args.stdname = stdname
        # have we been processed correctly?
        if 'DARCOR' in self.action.args.ccddata.header:
            if self.action.args.ccddata.header['DARCOR']:
                # does the standard file exist?
                if self.action.args.stdfile is not None:
                    # have we been stacked?
                    if 'NSTACK' in self.action.args.ccddata.header:
                        nstack = self.action.args.ccddata.header['NSTACK']
                    else:
                        nstack = 0
                    # does file already exist?
                    ofn = self.action.args.name
                    msname = strip_fname(ofn) + '_invsens.fits'
                    rdir = self.config.instrument.output_directory
                    invsensf = os.path.join(self.config.instrument.cwd,
                                            rdir,
                                            msname)
                    if os.path.exists(invsensf):
                        if nstack <= 0:
                            self.logger.warning("Master cal already exists: %s"
                                                % invsensf)
                            return False
                        else:
                            os.unlink(invsensf)
                            self.logger.info("Master cal will re-generated "
                                             "from stacked image")
                            return True
                    else:
                        self.logger.info("Master cal will be generated.")
                        return True
                else:
                    self.logger.warning("Not a KCWI standard star observation.")
                    return False
            else:
                self.logger.warning("DAR not corrected, cannot generate "
                                    "inverse sensitivity")
        else:
            self.logger.warning("Not processed enough to generate "
                                "inverse sensitivity")
            return False

    def _perform(self):
        self.logger.info("Making inverse sensitivity curve")

        suffix = 'invsens'
        stdname = self.action.args.stdname

        do_plots = self.config.instrument.plot_level >= 3

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
        # get all-inclusive wavelength range
        wall0 = self.action.args.ccddata.header['WAVALL0']
        wall1 = self.action.args.ccddata.header['WAVALL1']
        # get DAR padding in y
        pad_y = self.action.args.ccddata.header['DARPADY']
        # get sky subtraction status
        if 'SKYCOR' in self.action.args.ccddata.header:
            skycor = self.action.args.ccddata.header['SKYCOR']
        else:
            skycor = False
        # get wavelength correction status
        if 'CTYPE3' in self.action.args.ccddata.header:
            wavecor = ('WAVE' in self.action.args.ccddata.header['CTYPE3'])
        else:
            wavecor = False
        # get telescope and atm. correction
        if 'TELESCOP' in self.action.args.ccddata.header:
            tel = self.action.args.ccddata.header['TELESCOP']
        else:
            tel = 'KeckI'
        if 'Keck' in tel:
            area = 760000.0
        else:
            area = -1.0
        # compute good y pixel ranges
        if w0 > 0. and dw > 0. and wgoo0 > 0. and wgoo1 > 0.:
            z0 = int((wgoo0 - w0) / dw) + 10
            z1 = int((wgoo1 - w0) / dw) - 10
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
            # plot slice data, if requested
            if do_plots:
                py = tot[:, i]
                px = np.arange(len(py))
                p = figure(
                    title=self.action.args.stdlabel +
                    ' Std Slice %d (DAR padded)' % i,
                    x_axis_label="Position along slice",
                    y_axis_label="Flux summed over WLs",
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                p.scatter(px, py, marker='x')
                bokeh_plot(p, self.context.bokeh_session)
                print("sl, std: %02d, %.3f" % (i, tstd))
                qstr = input("Next? <cr>, q - quit: ")
                if 'Q' in qstr.upper():
                    do_plots = False

        # relevant slices
        if self.action.args.ifunum == 1:    # Large slicer
            slset = 3
        elif self.action.args.ifunum == 2:  # Medium slicer
            slset = 5
        elif self.action.args.ifunum == 3:  # Small slicer
            slset = 12
        else:
            slset = 3
            self.logger.error("Assuming Large slicer: id undefined! - %d" %
                              self.action.args.ifunum)
        # Get set of slices needed for calculation
        sl0 = (mxsl - slset) if mxsl >= slset else 0
        sl1 = (mxsl + slset) if (mxsl + slset) <= sz[2]-1 else sz[2]-1
        # get y position of std
        cy, _ = find_peaks(tot[:, mxsl], height=np.nanmean(tot[:, mxsl]))
        cy = int(cy[0]) + gy0
        # log results
        self.logger.info("Std slices (DAR padded): "
                         "max, sl0, sl1, spatial cntrd: "
                         "%d, %d, %d, %.2f" % (mxsl, sl0, sl1, cy))
        # get dwave spectrum
        ofn = self.action.args.name
        delfn = strip_fname(ofn) + '_dcubed.fits'
        full_path = os.path.join(
            self.config.instrument.cwd,
            self.config.instrument.output_directory, delfn)
        if os.path.exists(full_path):
            dew = kcwi_fits_reader(full_path)[0]
            dwspec = dew.data[:, cy, mxsl]
            zeros = np.where(dwspec == 0)
            if len(zeros) > 0:
                dwspec[zeros] = dw
        else:
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
        # check for zeros
        zeros = np.where(obsspec == 0.)
        if len(zeros) > 0:
            obsmean = np.nanmean(obsspec)
            obsspec[zeros] = obsmean
        # read in standard star reference spectrum
        stdfile = self.action.args.stdfile
        hdul = pf.open(stdfile)
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
        # how many points?
        self.logger.info("Number of standard points = %d" % nsroi)
        # very sparsely sampled w.r.t. object
        if nsroi <= 5:
            self.logger.error("Not enough standard points")
            return self.action.args

        swl = swl[sroi]
        sflx = sflx[sroi]
        sfw = sfw[sroi]
        fwhm = np.max(sfw)
        self.logger.info("Reference spectrum FWHM used = %.1f (A)" % fwhm)
        # resample standard onto our wavelength grid
        rsint = interp1d(swl, sflx, kind='cubic', fill_value='extrapolate')
        rsflx = rsint(w)
        # get effective inverse sensitivity
        invsen = rsflx / obsspec
        # convert to photons/s/cm^2/(wl bin = dw)
        rspho = 5.03411250e7 * rsflx * w * dw
        # get effective area
        earea = ubsspec / rspho
        # correct to native bins
        earea *= dw/dwspec
        # Line masks
        lmasks = []
        # Header for line mask file
        lmhdr = []
        # default values (for BM)
        ford = 9        # fit order
        if 'BL' in self.action.args.grating:
            ford = 7
        elif 'BH' in self.action.args.grating:
            ford = 9
        # Adjust for dichroic fraction
        try:
            dichroic_fraction = self.action.args.ccddata.header['DICHFRAC']
        except KeyError:
            dichroic_fraction = 1.
        ford = int(ford * dichroic_fraction)
        if ford < 3:
            ford = 3
        self.logger.info("Fitting Invsens with polynomial order %d" % ford)
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
        if self.config.instrument.plot_level > 1:
            print("CHECKING WAVELENGTH LIMITS")
            print("Current WL limits: %.1f - %.1f Angstroms "
                  "(blue vertical lines)" % (wlm0, wlm1))
            print("A <cr> will accept current values, "
                  "or enter new values to avoid edge problems (if present).")
            print("Hover the cursor over spectrum to find wavelengths.")
            # yran = [np.min(obsspec), np.max(obsspec)]
            yran = [np.min(obsspec[wl_good]), np.max(obsspec[wl_good])]
            source = ColumnDataSource(data=dict(x=w, y=obsspec))
            done = False
            while not done:
                p = figure(
                    tooltips=[("x", "@x{0,0.0}"), ("y", "@y{0,0.0}")],
                    title=self.action.args.stdlabel + ' Obs Spec',
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
                qstr = input("New values? <float> <float> (or <cr> - done): ")
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

        # Look for mask data file
        # first in local directory
        local_lmfile = os.path.basename(stdfile).split('.fits')[0] + '.lmsk'
        if not os.path.exists(local_lmfile):
            # then in stds directory
            # lmfile = stdfile.split('.fits')[0] + '.lmsk'
            lmfile = strip_fname(stdfile) + '.lmsk'
            if not os.path.exists(lmfile):
                lmfile = None
        else:
            lmfile = local_lmfile
        # if we found a mask file, read it in
        if lmfile is None:
            self.logger.info("No line mask file found")
        else:
            self.logger.info("Using line mask file %s" % lmfile)
            with open(lmfile) as lmf:
                lmasks_str = lmf.readlines()
            # parse file into mask list
            for lmws in lmasks_str:
                lmws = lmws.strip()
                # Collect header lines
                if lmws.startswith('#'):
                    lmhdr.append(lmws)
                    continue
                # Skip blank lines
                if len(lmws.strip()) < 1:
                    continue
                try:
                    data = lmws.split('#')[0]
                    # Parse comment on line
                    if '#' in lmws:
                        comm = lmws.split('#')[1].lstrip()
                    else:
                        comm = ''
                    # Parse line mask range
                    lm = [float(v) for v in data.split()]
                except ValueError:
                    print("bad line: %s" % lmws)
                    continue
                # Only collect good lines
                if len(lm) == 2 and lm[0] < lm[1]:
                    lmdict = {'w0': lm[0], 'w1': lm[1], 'com': comm}
                    lmasks.append(lmdict)
                else:
                    print("bad line: %s" % lmws)
        # Now interactively identify lines if requested
        if self.config.instrument.plot_level >= 2:
            yran = [np.min(obsspec[wl_good]), np.max(obsspec[wl_good])]
            # source = ColumnDataSource(data=dict(x=w, y=obsspec))
            print("MASKING SHARP FEATURES: ABSORPTION LINES/COSMIC RAYS")
            print("To mask, enter starting and stopping wavelengths "
                  "followed by an optional comment (string) for each line "
                  "you want to mask.")
            print("Current masks are shown as vertical dashed yellow lines.")
            print("Hover the cursor over spectrum to find wavelengths.")
            print("Mask values and comments are saved to this file: %s." %
                  local_lmfile)
            print("A <cr> with no line limits will accept current masks.")
            done = False
            while not done:
                p = figure(
                    tooltips=[("x", "@x{0.0}"), ("y", "@y{0.0}")],
                    title=self.action.args.stdlabel + ' Obs Spec',
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
                for ml in lmasks:
                    if wall0 < ml['w0'] < wall1 and wall0 < ml['w1'] < wall1:
                        p.line([ml['w0'], ml['w0']], yran, line_color='orange',
                               line_dash='dashed')
                        p.line([ml['w1'], ml['w1']], yran, line_color='orange',
                               line_dash='dashed')
                set_plot_lims(p, xlim=[wall0, wall1], ylim=yran)
                bokeh_plot(p, self.context.bokeh_session)
                qstr = input("Mask line: wavelength start stop (Ang) comment "
                             "(str) <float> <float> <str> (<cr> - done): ")
                if len(qstr) <= 0:
                    done = True
                else:
                    try:
                        newl = {'w0': float(qstr.split()[0]),
                                'w1': float(qstr.split()[1]),
                                'com': " ".join(qstr.split()[2:]).strip()}
                        # newl = [float(val) for val in qstr.split()]
                    except ValueError:
                        print("bad line: %s" % qstr)
                        continue
                    if wlm0 < newl['w0'] < wlm1 and wlm0 < newl['w1'] < wlm1:
                        lmasks.append(newl)
                    else:
                        print("line mask outside wl range: %s" % qstr)
        # END: interactively identify lines
        # Write out a local copy of line masks so user can edit
        if len(lmasks) > 0:
            local_lmfile = os.path.basename(
                stdfile).split('.fits')[0] + '.lmsk'
            with open(local_lmfile, 'w') as lmf:
                for lmh in lmhdr:
                    lmf.write(lmh)
                for lm in lmasks:
                    if len(lm['com']) == 0:
                        lmf.write("%.2f %.2f\n" % (lm['w0'], lm['w1']))
                    else:
                        lmf.write("%.2f %.2f # %s\n" % (lm['w0'], lm['w1'],
                                                        lm['com']))
        # set up fitting vectors, flux, waves, measure errors
        sf = invsen[wl_good]   # dependent variable
        af = earea[wl_good]    # effective area
        wf = w[wl_good]        # independent variable
        bf = obsspec[wl_good]  # input electrons
        rf = rsflx[wl_good]    # reference flux
        mw = np.ones(nwl_good, dtype=float)     # weights
        use = np.ones(nwl_good, dtype=int)      # toggles for usage
        # loop over line masks and apply them
        for ml in lmasks:
            roi = [i for i, v in enumerate(wf) if ml['w0'] <= v <= ml['w1']]
            nroi = len(roi)
            if nroi > 0:
                use[roi] = 0
                self.logger.info("Masking line at %.2f to %.2f (A)"
                                 % (ml['w0'], ml['w1']))
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
        resid = 100.0 * (scalspec - rf) / rf
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
                effmax = np.nanmax(100.*fearea/area)
                effmean = np.nanmean(100.*fearea/area)
                peff = figure(
                    title=self.action.args.stdlabel + ' Efficiency',
                    x_axis_label='Wave (A)',
                    y_axis_label='Effective Efficiency (%)',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                peff.line(wf, 100.*af/area, line_color='black',
                          legend_label='Data')
                peff.line(w, 100.*fearea/area, line_color='red',
                          legend_label='Fit')
                peff.scatter(ww, ef, marker='x', legend_label='Rej')
                peff.line([wlm0, wlm0], yran, line_color='green',
                          legend_label='WL lims')
                peff.line([wlm1, wlm1], yran, line_color='green')
                peff.line([wall0, wall1], [effmax, effmax],
                          line_color='black', line_dash='dashed')
                peff.line([wall0, wall1], [effmean, effmean],
                          line_color='black', line_dash='dashdot')
                set_plot_lims(peff, xlim=[wall0, wall1], ylim=yran)
                bokeh_plot(peff, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(2. * self.config.instrument.plot_pause)

                yran = [np.min(sf), np.max(sf)]
                pivs = figure(
                    title=self.action.args.stdlabel + ' Inverse sensitivity',
                    x_axis_label='Wave (A)',
                    y_axis_label='Invserse Sensitivity (Flux/e-/s)',
                    y_axis_type='log',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                pivs.line(wf, sf, line_color='black', legend_label='Data')
                pivs.line(w, finvsen, line_color='red', legend_label='Fit')
                pivs.scatter(ww, mf, marker='x', legend_label='Rej')
                pivs.line([wlm0, wlm0], yran, line_color='green',
                          legend_label='WL lims')
                pivs.line([wlm1, wlm1], yran, line_color='green')
                set_plot_lims(pivs, xlim=[wall0, wall1], ylim=yran)
                bokeh_plot(pivs, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(2. * self.config.instrument.plot_pause)

                yran = [np.min(calspec[wl_good]), np.max(calspec[wl_good])]
                pcal = figure(
                    title=self.action.args.stdlabel + ' Calibrated',
                    x_axis_label='Wave (A)',
                    y_axis_label='Flux (ergs/s/cm^2/A)',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                pcal.line(w, calspec, line_color='black', legend_label='Obs')
                pcal.line(w, rsflx, line_color='red', legend_label='Ref')
                pcal.line([wlm0, wlm0], yran, line_color='green',
                          legend_label='WL lims')
                pcal.line([wlm1, wlm1], yran, line_color='green')
                set_plot_lims(pcal, xlim=[wall0, wall1], ylim=yran)
                bokeh_plot(pcal, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(2. * self.config.instrument.plot_pause)

                yran = [np.min(resid), np.max(resid)]
                prsd = figure(
                    title=self.action.args.stdlabel +
                    ' Residuals = %.1f +- %.1f (%%)' % (rsd_mean, rsd_stdv),
                    x_axis_label='Wave (A)',
                    y_axis_label='Obs - Ref / Ref (%)',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                prsd.line(wf, resid, line_color='black',
                          legend_label='Obs - Ref (%)')
                if len(not_used) > 0:
                    prsd.scatter(ww, rbad, marker='x', legend_label='Rej')
                prsd.line([wlm0, wlm0], yran, line_color='green',
                          legend_label='WL lims')
                prsd.line([wlm1, wlm1], yran, line_color='green')
                prsd.line([wall0, wall1], [rsd_mean, rsd_mean],
                          line_color='red')
                prsd.line([wall0, wall1],
                          [rsd_mean+rsd_stdv, rsd_mean+rsd_stdv],
                          line_color='black', line_dash='dashed')
                prsd.line([wall0, wall1],
                          [rsd_mean-rsd_stdv, rsd_mean-rsd_stdv],
                          line_color='black', line_dash='dashed')
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
                            resid = 100.0 * (scalspec - rf) / rf
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
            # log results
            effmax = float(np.nanmax(100. * fearea / area))
            effmean = float(np.nanmean(100. * fearea / area))
            self.logger.info("Peak, mean efficiency (%%): %.1f, %.1f" %
                             (effmax, effmean))
            self.logger.info("Fit order = %d" % ford)
            # output plots
            pfname = "std_%05d_%s_%s_%s_%d" % (
                self.action.args.ccddata.header['FRAMENO'],
                stdname, self.action.args.grating.strip(),
                self.action.args.ifuname.strip(), int(self.action.args.cwave))
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
        hdr['INVFW0'] = (wlm0, 'low wavelength for fits')
        hdr['INVFW1'] = (wlm1, 'high wavelength for fits')
        hdr['INVFORD'] = (ford, 'fit order')
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
        if wavecor:
            hdr['CTYPE1'] = ('WAVE', 'Vacuum Wavelengths')
        else:
            hdr['CTYPE1'] = ('AWAV', 'Air Wavelengths')
        hdr['CUNIT1'] = ('Angstrom', 'Wavelength units')
        hdr['CNAME1'] = ('KCWI INVSENS Wavelength', 'Wavelength name')
        hdr['CRVAL1'] = (w0, 'Wavelength zeropoint')
        hdr['CRPIX1'] = (crpixw, 'Wavelength reference pixel')
        hdr['CDELT1'] = (dw, 'Wavelength Angstroms per pixel')

        ofn = self.action.args.name
        invsname = strip_fname(ofn) + '_' + suffix + '.fits'
        eaname = strip_fname(ofn) + '_ea.fits'

        # set units
        invsens_u = u.erg / (u.angstrom * u.cm ** 2 * u.s * u.electron)
        # output inverse sensitivity
        out_invsens = CCDData(np.asarray([invsen, finvsen, obsspec]),
                              meta=hdr, unit=invsens_u)
        kcwi_fits_writer(out_invsens, output_file=invsname,
                         output_dir=self.config.instrument.output_directory)
        self.context.proctab.update_proctab(frame=out_invsens, suffix=suffix,
                                            newtype='INVSENS',
                                            filename=self.action.args.name)
        # output effective area
        ea_u = u.cm ** 2 / u.angstrom
        out_ea = CCDData(np.asarray([earea, fearea]), meta=hdr, unit=ea_u)
        kcwi_fits_writer(out_ea, output_file=eaname,
                         output_dir=self.config.instrument.output_directory)
        self.context.proctab.update_proctab(frame=out_ea, suffix='ea',
                                            newtype='EAREA',
                                            filename=self.action.args.name)

        return self.action.args
    # END: class MakeInvsens()
