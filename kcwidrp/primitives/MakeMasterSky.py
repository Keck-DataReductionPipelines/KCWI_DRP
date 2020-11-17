from keckdrpframework.primitives.base_img import BaseImg
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_reader, \
    kcwi_fits_writer
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_plotting import save_plot
from kcwidrp.core.bspline import Bspline
from bokeh.plotting import figure

import os
import time
import numpy as np
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

        suffix = self.action.args.new_type.lower()
        ofn = self.action.args.ccddata.header['OFNAME']
        rdir = self.config.instrument.output_directory

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
                if ofn in row.split()[0]:
                    skyfile = row.split()[1]
                    self.logger.info("Found sky entry for %s: %s" % (ofn,
                                                                     skyfile))
                    if len(row.split()) > 2:
                        skymask = row.split()[2]
                        self.logger.info("Found sky mask entry for %s: %s"
                                         % (ofn, skymask))
            if skymask:
                if os.path.exists(skymask):
                    self.logger.info("Using sky mask file: %s" % skymask)
                else:
                    self.logger.warning("Sky mask file not found: %s" % skymask)
                    skymask = None
        self.action.args.skyfile = skyfile
        self.action.args.skymask = skymask
        if skyfile:
            msname = skyfile.split('.fits')[0] + '_' + suffix + '.fits'
            mskyf = os.path.join(rdir, msname)
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

        suffix = self.action.args.new_type.lower()

        # get root for maps
        tab = self.context.proctab.n_proctab(
            frame=self.action.args.ccddata, target_type='ARCLAMP',
            target_group=self.action.args.groupid)
        if len(tab) <= 0:
            self.logger.error("Geometry not solved!")
            return self.action.args

        groot = tab['OFNAME'][0].split('.fits')[0]

        # Wavelength map image
        wmf = groot + '_wavemap.fits'
        self.logger.info("Reading image: %s" % wmf)
        wavemap = kcwi_fits_reader(
            os.path.join(os.path.dirname(self.action.args.name), 'redux',
                         wmf))[0]

        # Slice map image
        slf = groot + '_slicemap.fits'
        self.logger.info("Reading image: %s" % slf)
        slicemap = kcwi_fits_reader(
            os.path.join(os.path.dirname(self.action.args.name), 'redux',
                         slf))[0]

        # Position map image
        pof = groot + '_posmap.fits'
        self.logger.info("Reading image: %s" % pof)
        posmap = kcwi_fits_reader(
            os.path.join(os.path.dirname(self.action.args.name), 'redux',
                         pof))[0]
        posmax = np.nanmax(posmap.data)
        posbuf = int(10. / self.action.args.xbinsize)

        # wavelength region
        wavegood0 = wavemap.header['WAVGOOD0']
        wavegood1 = wavemap.header['WAVGOOD1']
        waveall0 = wavemap.header['WAVALL0']
        waveall1 = wavemap.header['WAVALL1']

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
                     finiteflux[i] and not binary_mask.flat[i]]
        else:
            q = [i for i, v in enumerate(slicemap.data.flat)
                 if 0 <= v <= 23 and
                 posbuf < posmap.data.flat[i] < (posmax - posbuf) and
                 waveall0 <= wavemap.data.flat[i] <= waveall1 and
                 finiteflux[i] and not binary_mask.flat[i]]

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
        ofn = self.action.args.ccddata.header['OFNAME']
        msname = ofn.split('.fits')[0] + '_' + suffix + '.fits'

        log_string = MakeMasterSky.__module__
        self.action.args.ccddata.header['IMTYPE'] = self.action.args.new_type
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
                                            newtype=self.action.args.new_type)
        self.context.proctab.write_proctab()

        # restore original image
        self.action.args.ccddata.data = img
        self.action.args.ccddata.header = hdr

        self.logger.info(log_string)
        return self.action.args

    # END: class MakeMasterSky()
