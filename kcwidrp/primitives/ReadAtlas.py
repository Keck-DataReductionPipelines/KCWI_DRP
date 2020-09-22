from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot

from bokeh.plotting import figure
from bokeh.models import Range1d
import pkg_resources
import os
from astropy.io import fits as pf
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interpolate
from scipy import signal
import time


class ReadAtlas(BasePrimitive):
    """Read in atlas spectrum and derive alignment offset"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        # What lamp are we using?
        lamp = self.action.args.illum
        # Does the atlas file exist?
        path = "data/%s.fits" % lamp.lower()  # always use slash
        pkg = __name__.split('.')[0]
        atpath = pkg_resources.resource_filename(pkg, path)
        if os.path.exists(atpath):
            self.logger.info("Reading atlas spectrum in: %s" % atpath)
        else:
            self.logger.error("Atlas spectrum not found for %s" % atpath)
        # Read the atlas
        ff = pf.open(atpath)
        reflux = ff[0].data
        refdisp = ff[0].header['CDELT1']
        refwav = np.arange(0, len(reflux)) * refdisp + ff[0].header['CRVAL1']
        ff.close()
        # Convolve with appropriate Gaussian
        self.logger.info("Convolving Atlas with Gaussian having sigma of %.2f "
                         "px" % self.action.args.atsig)
        reflux = gaussian_filter1d(reflux, self.action.args.atsig)
        # Observed arc spectrum
        obsarc = self.context.arcs[self.config.instrument.REFBAR]
        # Preliminary wavelength solution
        xvals = np.arange(0, len(obsarc)) - int(len(obsarc)/2)
        obswav = xvals * self.context.prelim_disp + self.action.args.cwave
        # Get central third
        minow = int(len(obsarc)/3)
        maxow = int(2.*len(obsarc)/3)
        # Unless we are low dispersion, then get central 3 5ths
        if 'BL' in self.action.args.grating or 'RL' in self.action.args.grating:
            minow = int(len(obsarc)/5)
            maxow = int(4.*len(obsarc)/5)
        minwav = obswav[minow]
        maxwav = obswav[maxow]
        # Check for dichroic
        if self.action.args.dich:
            obs_extent = maxow - minow
            if self.action.args.camera == 0:  # Blue
                # test if correction needed
                if maxwav > 5600:
                    maxwav = [w for w in obswav if w < 5620.][-1]
                    maxow = [i for i, w in enumerate(obswav) if w < 5620.][-1]
                    minow = maxow - obs_extent
                    if minow < 0:
                        minow = 0
                    self.logger.info("Dichroic adjusted - x0, x1, w0, w1: "
                                     "%d, %d, %.2f, %.2f" % (minow, maxow,
                                                             minwav, maxwav))
            elif self.action.args.camera == 1:  # Red
                # test if correction needed
                if minwav < 5600:
                    minwav = [w for w in obswav if w > 5580.][0]
                    minow = [i for i, w in enumerate(obswav) if w > 5580.][0]
                    maxow = minow + obs_extent
                    if maxow > (len(obsarc) - 1):
                        maxow = len(obsarc) - 1
                    self.logger.info("Dichroic adjusted - x0, x1, w0, w1: "
                                     "%d, %d, %.2f, %.2f" % (minow, maxow,
                                                             minwav, maxwav))
            else:
                self.logger.warning("Camera undefined!!")
        # Get corresponding ref range
        minrw = [i for i, v in enumerate(refwav) if v >= minwav][0]
        maxrw = [i for i, v in enumerate(refwav) if v <= maxwav][-1]
        # Subsample for cross-correlation
        cc_obsarc = obsarc[minow:maxow]
        cc_obswav = obswav[minow:maxow]
        cc_reflux = reflux[minrw:maxrw]
        cc_refwav = refwav[minrw:maxrw]
        # Resample onto reference wavelength scale
        obsint = interpolate.interp1d(cc_obswav, cc_obsarc, kind='cubic',
                                      bounds_error=False,
                                      fill_value='extrapolate'
                                      )
        cc_obsarc = obsint(cc_refwav)
        # Apply cosign bell taper to both
        cc_obsarc *= signal.windows.tukey(
            len(cc_obsarc), alpha=self.config.instrument.TAPERFRAC)
        cc_reflux *= signal.windows.tukey(
            len(cc_reflux), alpha=self.config.instrument.TAPERFRAC)
        nsamp = len(cc_refwav)
        offar = np.arange(1 - nsamp, nsamp)
        # Cross-correlate
        xcorr = np.correlate(cc_obsarc, cc_reflux, mode='full')
        # Get central region
        x0c = int(len(xcorr)/3)
        x1c = int(2*(len(xcorr)/3))
        xcorr_central = xcorr[x0c:x1c]
        offar_central = offar[x0c:x1c]
        # Calculate offset
        offset_pix = offar_central[xcorr_central.argmax()]
        offset_wav = offset_pix * refdisp
        self.logger.info("Initial arc-atlas offset (px, Ang): %d, %.1f" %
                         (offset_pix, offset_wav))
        if self.config.instrument.plot_level >= 1:
            # Plot
            p = figure(title=self.action.args.plotlabel +
                       "ATLAS OFFSET = %d px" % offset_pix,
                       x_axis_label="Offset(px)", y_axis_label="X-corr",
                       plot_width=self.config.instrument.plot_width,
                       plot_height=self.config.instrument.plot_height)

            p.line(offar_central, xcorr_central, legend_label='Data')
            ylim_min = min(xcorr_central)
            ylim_max = max(xcorr_central)
            p.line([offset_pix, offset_pix], [ylim_min, ylim_max],
                   color='red', legend_label='Peak')
            bokeh_plot(p, self.context.bokeh_session)
            if self.config.instrument.plot_level >= 2:
                input("Next? <cr>: ")
            else:
                time.sleep(self.config.instrument.plot_pause)
            # Get central wavelength
            cwave = self.action.args.cwave
            # Set up offset tweaking
            q = 'test'
            while q:
                # Plot the two spectra
                p = figure(title=self.action.args.plotlabel +
                           "ATLAS OFFSET = %.1f Ang (%d px)" %
                           (offset_wav, offset_pix),
                           x_axis_label="Wave(A)", y_axis_label="Rel. Flux",
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.line(obswav[minow:maxow] - offset_wav,
                       obsarc[minow:maxow]/np.nanmax(obsarc[minow:maxow]),
                       legend_label="ref bar (%d)" %
                       self.config.instrument.REFBAR)
                p.line(refwav[minrw:maxrw],
                       reflux[minrw:maxrw]/np.nanmax(reflux[minrw:maxrw]),
                       color="red", legend_label="Atlas")
                p.x_range = Range1d(minwav, maxwav)
                ylim_min = min(
                    obsarc[minow:maxow]/np.nanmax(obsarc[minow:maxow]))
                ylim_max = max(
                    obsarc[minow:maxow]/np.nanmax(obsarc[minow:maxow]))
                p.line([cwave, cwave], [ylim_min, ylim_max], color="magenta",
                       legend_label="CWAVE", line_dash="dashdot")
                bokeh_plot(p, self.context.bokeh_session)

                if self.config.instrument.plot_level >= 2:
                    q = input("Enter: <cr> - next, new offset (int px): ")
                    if q:
                        try:
                            offset_pix = int(q)
                            offset_wav = offset_pix * refdisp
                        except ValueError:
                            print("Try again: integer pixel values accepted")
                            q = 'test'
                else:
                    time.sleep(self.config.instrument.plot_pause)
                    q = None
            self.logger.info("Final   arc-atlas offset (px, Ang): %d, %.1f" %
                             (offset_pix, offset_wav))
        # Store atlas spectrum
        self.action.args.reflux = reflux
        self.action.args.refwave = refwav
        # Store offsets
        self.action.args.offset_pix = offset_pix
        self.action.args.offset_wave = offset_wav
        # Store reference dispersion
        self.action.args.refdisp = refdisp
        # Store central limits
        self.action.args.minrow = minow
        self.action.args.maxrow = maxow
        # Store x values
        self.action.args.xvals = xvals
        self.action.args.x0 = int(len(obsarc)/2)

        log_string = ReadAtlas.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class ReadAtlas()
