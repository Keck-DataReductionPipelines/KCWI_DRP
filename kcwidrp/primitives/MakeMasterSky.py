from keckdrpframework.primitives.base_img import BaseImg
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer, strip_fname
from kcwidrp.primitives.GetAtlasLines import gaus
from kcwidrp.core.kcwi_get_std import kcwi_get_std
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_plotting import save_plot
from kcwidrp.core.bspline import Bspline
from bokeh.plotting import figure

import os
import time
import numpy as np
from scipy.optimize import curve_fit
from astropy.io import fits


class MakeMasterSky(BaseImg):
    """Make master sky image"""

    def __init__(self, action, context):
        BaseImg.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """
        Checks if we can create a master sky
        :return:
        """
        self.logger.info("Checking precondition for MakeMasterSky")

        suffix = 'sky'  # self.action.args.new_type.lower()
        ofn = self.action.args.name
        rdir = self.config.instrument.output_directory

        # Are we a standard star?
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

        # Is there a kcwi.sky file?
        skyfile = None
        skymask = None
        # check if kcwi.sky exists
        if os.path.exists('kcwi.sky'):
            self.logger.info("Reading kcwi.sky")
            f = open('kcwi.sky')
            skyproc = f.readlines()
            f.close()
            # is our file in the list?
            for row in skyproc:
                # skip comments
                if row.startswith('#'):
                    continue
                # Parse row:
                # <raw sci file> <raw sky file> <optional mask file>
                #  OR
                # <raw sci file> skip
                # to disable sky subtraction
                # Find match to current file
                if ofn in row.split()[0]:
                    skyfile = row.split()[1]
                    # Should we skip sky subtraction?
                    if 'skip' in skyfile:
                        self.logger.info("Skipping sky subtraction for %s" %
                                         ofn)
                        keycom = 'sky corrected?'
                        self.action.args.ccddata.header['SKYCOR'] = (False,
                                                                     keycom)
                        return False
                    self.logger.info("Found sky entry for %s: %s" % (ofn,
                                                                     skyfile))
                    # Do we have an optional sky mask file?
                    if len(row.split()) > 2:
                        skymask = row.split()[2]
                        self.logger.info("Found sky mask entry for %s: %s"
                                         % (ofn, skymask))
            # Do have a mask file?
            if skymask:
                # Does it exist?
                if os.path.exists(skymask):
                    self.logger.info("Using sky mask file: %s" % skymask)
                else:
                    self.logger.warning("Sky mask file not found: %s" % skymask)
                    skymask = None
        # Record results
        self.action.args.skyfile = skyfile
        self.action.args.skymask = skymask
        # Do we have a sky alternate?
        if skyfile:
            # Generate sky file name
            msname = skyfile.split('.fits')[0] + '_' + suffix + '.fits'
            mskyf = os.path.join(rdir, msname)
            # Does it exist?
            if os.path.exists(mskyf):
                self.logger.info("Master sky already exists: %s" % mskyf)
                return False
            else:
                self.logger.warning("Alternate master sky %s not found."
                                    % mskyf)
                return True
        else:
            self.logger.info("No alternate master sky requested.")
            return True

    def _perform(self):
        """
        Returns an Argument() with the parameters that depends on this operation
        """
        self.logger.info("Creating master sky")

        suffix = 'sky'

        # get root for maps
        tab = self.context.proctab.search_proctab(
            frame=self.action.args.ccddata, target_type='MARC',
            target_group=self.action.args.groupid)
        if len(tab) <= 0:
            self.logger.error("Geometry not solved!")
            return self.action.args

        groot = strip_fname(tab['filename'][0])

        # Wavelength map image
        wmf = groot + '_wavemap.fits'
        self.logger.info("Reading image: %s" % wmf)
        wavemap = kcwi_fits_reader(os.path.join(
            self.config.instrument.cwd, 'redux', wmf))[0]

        # Slice map image
        slf = groot + '_slicemap.fits'
        self.logger.info("Reading image: %s" % slf)
        slicemap = kcwi_fits_reader(os.path.join(
            self.config.instrument.cwd, 'redux', slf))[0]

        # Position map image
        pof = groot + '_posmap.fits'
        self.logger.info("Reading image: %s" % pof)
        posmap = kcwi_fits_reader(os.path.join(
            self.config.instrument.cwd, 'redux', pof))[0]
        posmax = np.nanmax(posmap.data)
        posbuf = int(10. / self.action.args.xbinsize)

        ny = posmap.data.shape[0]

        if self.action.args.camera == 1:    # Red
            # Get ymap for trimming junk at ends
            ymap = posmap.copy()
            for i in range(ny):
                ymap.data[i, :] = float(i)
        else:
            ymap = None

        # wavelength region
        wavegood0 = wavemap.header['WAVGOOD0']
        wavegood1 = wavemap.header['WAVGOOD1']
        waveall0 = wavemap.header['WAVALL0']
        waveall1 = wavemap.header['WAVALL1']
        wavemid = wavemap.header['WAVMID']

        # get image size
        sm_sz = self.action.args.ccddata.data.shape

        # sky masking
        # default is no masking (True = mask, False = don't mask)
        binary_mask = np.zeros(sm_sz, dtype=bool)

        # was sky masking requested?
        if self.action.args.skymask:
            if os.path.exists(self.action.args.skymask):
                self.logger.info("Reading sky mask file: %s"
                                 % self.action.args.skymask)
                hdul = fits.open(self.action.args.skymask)
                binary_mask = hdul[0].data
                # verify size match
                bm_sz = binary_mask.shape
                if bm_sz[0] != sm_sz[0] or bm_sz[1] != sm_sz[1]:
                    self.logger.warning("Sky mask size mis-match: "
                                        "masking disabled")
                    binary_mask = np.zeros(sm_sz, dtype=bool)
            else:
                self.logger.warning("Sky mask image not found: %s"
                                    % self.action.args.skymask)

        # if we are a standard, get mask for bright continuum source
        if self.action.args.stdname is not None:
            self.logger.info("Standard star observation of "
                             "%s will be auto-masked" %
                             self.action.args.stdname)
            # Use 10% of wavelength range at wavemid
            std_wav_ran = (wavemid - 0.05 * (wavegood1 - wavegood0),
                           wavemid + 0.05 * (wavegood1 - wavegood0))
            std_sl_max = -1
            std_sl_sig_max = -1.
            std_sl_max_pos_data = None
            std_sl_max_flx_data = None

            self.logger.info("Finding the std max slice")
            for si in range(24):
                sq = [i for i, v in enumerate(slicemap.data.flat) if v == si and
                      std_wav_ran[0] < wavemap.data.flat[i] < std_wav_ran[1] and
                      posbuf < posmap.data.flat[i] < (posmax - posbuf)]
                xplt = posmap.data.flat[sq]
                yplt = self.action.args.ccddata.data.flat[sq]
                sig = float(np.nanstd(yplt))
                self.logger.info("Slice %d - StDev = %.2f" % (si, sig))
                if sig > std_sl_sig_max:
                    std_sl_sig_max = sig
                    std_sl_max = si
                    std_sl_max_pos_data = xplt.copy()
                    std_sl_max_flx_data = yplt.copy()
            ipk = np.argmax(std_sl_max_flx_data)
            ppk = std_sl_max_pos_data[ipk]
            fpk = std_sl_max_flx_data[ipk]

            # gaussian fit to max slice
            res, _ = curve_fit(gaus, std_sl_max_pos_data, std_sl_max_flx_data,
                               p0=[fpk, ppk, 1.])
            self.logger.info("Std max at %.2f in slice %d with width %.2f px"
                             % (res[1], std_sl_max, res[2]))
            std_pos_mask_0 = res[1] - 5. * res[2]
            std_pos_mask_1 = res[1] + 5. * res[2]
            self.logger.info("Masking between %.2f and %.2f" %
                             (std_pos_mask_0, std_pos_mask_1))

            # Mask standard from sky calculation
            for i, v in enumerate(binary_mask.flat):
                if std_pos_mask_0 < posmap.data.flat[i] < std_pos_mask_1:
                    binary_mask.flat[i] = True

            # plot, if requested
            if self.config.instrument.plot_level >= 1:
                xx = np.arange(np.min(std_sl_max_pos_data),
                               np.max(std_sl_max_pos_data), 1)
                yy = gaus(xx, res[0], res[1], res[2])
                p = figure(
                    title=self.action.args.plotlabel +
                    ' Std max sl %d' % std_sl_max,
                    x_axis_label='Pos (x px)',
                    y_axis_label='Flux (e-)',
                    plot_width=self.config.instrument.plot_width,
                    plot_height=self.config.instrument.plot_height)
                p.circle(std_sl_max_pos_data, std_sl_max_flx_data, size=1,
                         line_alpha=0., fill_color='purple', legend_label='Data')
                p.line([ppk, ppk], [0, fpk], color='green')
                p.line([std_pos_mask_0, std_pos_mask_0], [0, fpk], color='blue')
                p.line([std_pos_mask_1, std_pos_mask_1], [0, fpk], color='blue')
                p.line(xx, yy, color='red')
                bokeh_plot(p, self.context.bokeh_session)
                if self.config.instrument.plot_level >= 2:
                    input("Next? <cr>: ")
                else:
                    time.sleep(self.config.instrument.plot_pause)

        # count masked pixels
        tmsk = len(np.nonzero(np.where(binary_mask.flat, True, False))[0])
        self.logger.info("Number of pixels masked = %d" % tmsk)

        finiteflux = np.isfinite(self.action.args.ccddata.data.flat)

        # get un-masked points mapped to exposed regions on CCD
        # handle dichroic bad region
        if self.action.args.dich:
            if self.action.args.camera == 0:    # Blue
                q = [i for i, v in enumerate(slicemap.data.flat)
                     if 0 <= v <= 23 and
                     posbuf < posmap.data.flat[i] < (posmax - posbuf) and
                     waveall0 <= wavemap.data.flat[i] <= waveall1 and
                     not (v > 20 and wavemap.data.flat[i] > 5600.) and
                     finiteflux[i] and not binary_mask.flat[i]]
            else:                               # Red
                q = [i for i, v in enumerate(slicemap.data.flat)
                     if 0 <= v <= 23 and
                     posbuf < posmap.data.flat[i] < (posmax - posbuf) and
                     waveall0 <= wavemap.data.flat[i] <= waveall1 and
                     not (v > 20 and wavemap.data.flat[i] < 5600.) and
                     finiteflux[i] and not binary_mask.flat[i] and
                     50 <= ymap.data.flat[i] <= (ny - 50)]
        else:
            if self.action.args.camera == 0:    # Blue
                q = [i for i, v in enumerate(slicemap.data.flat)
                     if 0 <= v <= 23 and
                     posbuf < posmap.data.flat[i] < (posmax - posbuf) and
                     waveall0 <= wavemap.data.flat[i] <= waveall1 and
                     finiteflux[i] and not binary_mask.flat[i]]
            else:
                q = [i for i, v in enumerate(slicemap.data.flat)
                     if 0 <= v <= 23 and
                     posbuf < posmap.data.flat[i] < (posmax - posbuf) and
                     waveall0 <= wavemap.data.flat[i] <= waveall1 and
                     finiteflux[i] and not binary_mask.flat[i] and
                     50 <= ymap.data.flat[i] <= (ny - 50)]

        # get all points mapped to exposed regions on the CCD (for output)
        qo = [i for i, v in enumerate(slicemap.data.flat)
              if 0 <= v <= 23 and posmap.data.flat[i] >= 0 and
              waveall0 <= wavemap.data.flat[i] <= waveall1 and
              finiteflux[i]]

        # extract relevant image values
        fluxes = self.action.args.ccddata.data.flat[q]

        # relevant wavelengths
        waves = wavemap.data.flat[q]
        self.logger.info("Number of fit waves = %d" % len(waves))

        # keep output wavelengths
        owaves = wavemap.data.flat[qo]
        self.logger.info("Number of output waves = %d" % len(owaves))

        # sort on wavelength
        s = np.argsort(waves)
        waves = waves[s]
        fluxes = fluxes[s]

        # knots per pixel
        knotspp = self.config.instrument.KNOTSPP
        n = int(sm_sz[0] * knotspp)

        # calculate break points for b splines
        bkpt = np.min(waves) + np.arange(n+1) * \
            (np.max(waves) - np.min(waves)) / n

        # log
        self.logger.info("Nknots = %d, min = %.2f, max = %.2f (A)" %
                         (n, np.min(bkpt), np.max(bkpt)))

        # do bspline fit
        sft0, gmask = Bspline.iterfit(waves, fluxes, fullbkpt=bkpt,
                                      upper=1, lower=1)
        gp = [i for i, v in enumerate(gmask) if v]
        yfit1, _ = sft0.value(waves)
        self.logger.info("Number of good points = %d" % len(gp))

        # check result
        if np.max(yfit1) < 0:
            self.logger.warning("B-spline failure")
            if n > 2000:
                if n == 5000:
                    n = 2000
                if n == 8000:
                    n = 5000
                # calculate breakpoints
                bkpt = np.min(waves) + np.arange(n + 1) * \
                    (np.max(waves) - np.min(waves)) / n
                # log
                self.logger.info("Nknots = %d, min = %.2f, max = %.2f (A)" %
                                 (n, np.min(bkpt), np.max(bkpt)))
                # do bspline fit
                sft0, gmask = Bspline.iterfit(waves, fluxes, fullbkpt=bkpt,
                                              upper=1, lower=1)
                yfit1, _ = sft0.value(waves)
            if np.max(yfit1) <= 0:
                self.logger.warning("B-spline final failure, sky is zero")

        # get values at original wavelengths
        yfit, _ = sft0.value(owaves)

        # for plotting
        gwaves = waves[gp]
        gfluxes = fluxes[gp]
        npts = len(gwaves)
        stride = int(npts / 8000.)
        xplt = gwaves[::stride]
        yplt = gfluxes[::stride]
        fplt, _ = sft0.value(xplt)
        yrng = [np.min(yplt), np.max(yplt)]
        self.logger.info("Stride = %d" % stride)

        # plot, if requested
        if self.config.instrument.plot_level >= 1:
            # output filename stub
            skyfnam = "sky_%05d_%s_%s_%s" % \
                     (self.action.args.ccddata.header['FRAMENO'],
                      self.action.args.illum, self.action.args.grating,
                      self.action.args.ifuname)
            p = figure(
                title=self.action.args.plotlabel + ' Master Sky',
                x_axis_label='Wave (A)',
                y_axis_label='Flux (e-)',
                plot_width=self.config.instrument.plot_width,
                plot_height=self.config.instrument.plot_height)
            p.circle(xplt, yplt, size=1, line_alpha=0., fill_color='purple',
                     legend_label='Data')
            p.line(xplt, fplt, line_color='red', legend_label='Fit')
            p.line([wavegood0, wavegood0], yrng, line_color='green')
            p.line([wavegood1, wavegood1], yrng, line_color='green')
            p.y_range.start = yrng[0]
            p.y_range.end = yrng[1]
            bokeh_plot(p, self.context.bokeh_session)
            if self.config.instrument.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.config.instrument.plot_pause)
            save_plot(p, filename=skyfnam+".png")

        # create sky image
        sky = np.zeros(self.action.args.ccddata.data.shape, dtype=float)
        sky.flat[qo] = yfit

        # store original data, header
        img = self.action.args.ccddata.data
        hdr = self.action.args.ccddata.header.copy()
        self.action.args.ccddata.data = sky

        # get master sky output name
        ofn_full = self.action.args.name
        ofn = os.path.basename(ofn_full)
        msname = strip_fname(ofn) + '_' + suffix + '.fits'

        log_string = MakeMasterSky.__module__
        self.action.args.ccddata.header['IMTYPE'] = 'SKY'
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.action.args.ccddata.header['SKYMODEL'] = (True, 'sky model image?')
        self.action.args.ccddata.header['SKYIMAGE'] = \
            (ofn, 'image used for sky model')
        if tmsk > 0:
            self.action.args.ccddata.header['SKYMSK'] = (True,
                                                         'was sky masked?')
            # self.action.args.ccddata.header['SKYMSKF'] = (skymf,
            # 'sky mask file')
        else:
            self.action.args.ccddata.header['SKYMSK'] = (False,
                                                         'was sky masked?')
        self.action.args.ccddata.header['WAVMAPF'] = wmf
        self.action.args.ccddata.header['SLIMAPF'] = slf
        self.action.args.ccddata.header['POSMAPF'] = pof

        # output master sky
        kcwi_fits_writer(self.action.args.ccddata, output_file=msname,
                         output_dir=self.config.instrument.output_directory)
        self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                            suffix=suffix,
                                            newtype="SKY",
                                            filename=self.action.args.name)
        self.context.proctab.write_proctab(tfil=self.config.instrument.procfile)

        # restore original image
        self.action.args.ccddata.data = img
        self.action.args.ccddata.header = hdr

        self.logger.info(log_string)
        return self.action.args

    # END: class MakeMasterSky()
