from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.core.bokeh_plotting import bokeh_plot
from kcwidrp.core.kcwi_plotting import save_plot
from kcwidrp.core.kcwi_pkg_resources import get_resource_path
from kcwidrp.primitives.kcwi_file_primitives import plotlabel
from kcwidrp.core.kcwi_interactive import ask_value, pause


from kcwidrp.core.kcwi_plotting import figure
from kcwidrp.core import kcwi_plot_colors as pc
from kcwidrp.core.kcwi_plotting import oplot_cwave
from kcwidrp.core.kcwi_plotting import Range1d
import os
from astropy.io import fits as pf
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interpolate
from scipy import signal


class ReadAtlas(BasePrimitive):
    """
    Read in atlas spectrum and derive alignment offset.

    Reads in spectral atlas spectrum of observed arc lamp.  Uses instrument
    configuration to determine convoution kernel that will match the observed
    spectrum resolution and applies it.  Then picks a central part of the
    observed spectrum and tapers the ends and does a cross-correlation to
    determine the pixel offset that will align the observed and atlas spectra.

    Uses the following configuration parameters:

        * MIDFRAC: what central fraction to use for finding the atlas offset - defaults to 0.3.
        * TAPERFRAC: how much of the ends to taper off (in case of bright lines at the ends) - defaults to 0.2.
        * ATOFF: force a specific pixel offset for the atlas - defaults to 0.

    Places convolved atlas spectrum with wavelength scale in the return
    argument so subsequent processing can use it to fit the observed wavelength
    scale accurately.

    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        self.logger.info("Checking for master arc")
        if 'MARC' in self.action.args.ccddata.header['IMTYPE']:
            return True
        else:
            return False

    def _perform(self):
        # What lamp are we using?
        lamp = self.action.args.illum
        # Does the atlas file exist?
        path = "data/%s.fits" % lamp.lower()  # always use slash
        pkg = __name__.split('.')[0]
        atpath = get_resource_path(pkg, path)
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
        # Define Central Region
        mf = self.config.instrument.MIDFRAC
        # Set on command line
        if 0 < mf <= 1.0:
            self.logger.info("Using command line value for middle fraction")
        # Default values based on grating
        else:
            self.logger.info("Using default calculated middle fraction")
            # Get central third
            mf = 1./3.
            # Unless we are low disp., then get central 3 5ths
            if 'BL' in self.action.args.grating or \
                    'RL' in self.action.args.grating:
                mf = 3./5.
            # Unless we are red high disp., then get central 8 10ths
            if 'RH' in self.action.args.grating:
                mf = 8./10.
        self.logger.info("Central %.3f used for initial fit" % mf)
        minow = int(len(obsarc) * (0.5 - mf / 2.))
        maxow = int(len(obsarc) * (0.5 + mf / 2.))
        if self.context.prelim_disp > 0:
            minwav = obswav[minow]
            maxwav = obswav[maxow]
        else:
            minwav = obswav[maxow]
            maxwav = obswav[minow]
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
                                      fill_value='extrapolate')
        cc_obsarc = obsint(cc_refwav)
        # Apply cosin bell taper to both
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
        calc_offset_pix = offar_central[xcorr_central.argmax()]
        calc_offset_wav = calc_offset_pix * refdisp
        self.logger.info("Calculated arc-atlas offset (px, Ang): %d, %.1f" %
                         (calc_offset_pix, calc_offset_wav))
        atoff = self.config.instrument.ATOFF
        req_offset_pix = 0
        if atoff != 0:
            req_offset_pix = atoff
            req_offset_wav = atoff * refdisp
            self.logger.info("Command line requested offset (px, Ang):"
                             " %d, %.1f" % (req_offset_pix, req_offset_wav))
            offset_pix = req_offset_pix
            offset_wav = req_offset_wav
        else:
            offset_pix = calc_offset_pix
            offset_wav = calc_offset_wav
        plab = plotlabel(self.action.args)
        if self.config.instrument.plot_level >= 1:
            # Plot
            p = figure(title=plab + "ATLAS OFFSET = %d px" % offset_pix,
                       x_axis_label="Offset(px)", y_axis_label="X-corr",
                       plot_width=self.config.instrument.plot_width,
                       plot_height=self.config.instrument.plot_height)

            p.line(offar_central, xcorr_central, color=pc.DATA,
                   legend_label='Data')
            ylim_min = min(xcorr_central)
            ylim_max = max(xcorr_central)
            p.line([calc_offset_pix, calc_offset_pix], [ylim_min, ylim_max],
                   color=pc.FOUND, legend_label='Peak')
            if atoff != 0:
                p.line([req_offset_pix, req_offset_pix], [ylim_min, ylim_max],
                       color=pc.PREDICTED, line_dash='dashed',
                       legend_label='Reqested')
            bokeh_plot(p, self.context.bokeh_session)
            pause(self.config.instrument, self.logger,
                  what="Cross-correlation of the reference-bar arc with the lamp atlas: offset found.",
                  help=("What you see:\n"
                        "  Title \"ATLAS OFFSET = <n> px\": cross-correlation of the reference-bar\n"
                        "  arc spectrum (on its calculated wavelengths) with the atlas spectrum,\n"
                        "  against offset in atlas pixels (central third of the range).\n"
                        "  Blue: cross-correlation.\n"
                        "  Magenta vertical (\"Peak\"): offset found.\n"
                        "  Black dashed vertical: offset requested with ATOFF (only if set).\n"
                        "What to look for:\n"
                        "  Good: one sharp peak clearly above the rest, with the magenta line on\n"
                        "  it.\n"
                        "  Bad: several similar peaks or a broad hump: the offset is ambiguous.\n"
                        "  Check the overlay in the next plot and set the offset there.\n"
                        "Your options:\n"
                        "  <cr>  continue to the overlay. Display only.\n"
                        "Default:\n"
                        "  <cr> continues. At plot_level 1 the plot is shown for plot_pause\n"
                        "  seconds; at plot_level 0 it is not drawn."),
                  who="ReadAtlas")
            # Get central wavelength
            cwave = self.action.args.cwave
            # Set up offset tweaking
            q = 'test'
            while q is not None:
                # Plot the two spectra
                p = figure(title=plab + "ATLAS OFFSET = %.1f Ang (%d px)" %
                           (offset_wav, offset_pix),
                           x_axis_label="Wave(A)", y_axis_label="Rel. Flux",
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                # atlas first, so the arc being inspected is drawn on top
                p.line(refwav[minrw:maxrw],
                       reflux[minrw:maxrw]/np.nanmax(reflux[minrw:maxrw]),
                       color=pc.REFERENCE, legend_label="Atlas")
                p.line(obswav[minow:maxow] - offset_wav,
                       obsarc[minow:maxow]/np.nanmax(obsarc[minow:maxow]),
                       color=pc.DATA, legend_label="ref bar (%d)" %
                       self.config.instrument.REFBAR)
                p.x_range = Range1d(minwav, maxwav)
                ylim_min = min(
                    obsarc[minow:maxow]/np.nanmax(obsarc[minow:maxow]))
                ylim_max = max(
                    obsarc[minow:maxow]/np.nanmax(obsarc[minow:maxow]))
                oplot_cwave(p, cwave, [ylim_min, ylim_max])
                bokeh_plot(p, self.context.bokeh_session)

                if self.config.instrument.plot_level >= 2:
                    q = ask_value(
                        self.logger, "New offset <int px> (<cr> = accept)",
                        help=("What you see:\n"
                              "  Title \"ATLAS OFFSET = <A> Ang (<n> px)\".\n"
                              "  Blue: reference-bar arc spectrum (bar REFBAR), normalised, on its\n"
                              "    calculated wavelengths minus the current offset.\n"
                              "  Black: atlas spectrum of the lamp, normalised.\n"
                              "  Orange dash-dot-dot: CWAVE.\n"
                              "What to look for:\n"
                              "  Good: blue and black lines coincide line for line over the range\n"
                              "  (relative strengths differ; small misalignments at the ends are\n"
                              "  fixed later, since the dispersion is not yet refined).\n"
                              "  Bad: a constant shift between the two patterns. Enter a new offset:\n"
                              "  a larger offset moves the blue spectrum to shorter wavelengths; the\n"
                              "  title gives the offset in A and px, i.e. the atlas A per px.\n"
                              "Your options:\n"
                              "  <cr>      accept the offset shown and continue.\n"
                              "  <number>  try that offset (integer atlas px); the plot is redrawn and\n"
                              "            you are asked again.\n"
                              "  The offset sets the starting wavelengths of every bar's fit, so it\n"
                              "  changes the whole wavelength solution and all cubes using this arc.\n"
                              "Default:\n"
                              "  <cr> keeps the offset shown: the calculated one (or ATOFF if set)\n"
                              "  unless you changed it. Below plot_level 2 that offset is used."),
                        parse=int, default=None,
                        effect="atlas offset set by user to {value} px "
                               "(calculated %d px)" % calc_offset_pix,
                        who="ReadAtlas",
                        what="Reference-bar arc over the lamp atlas: accept the offset or enter a new one (atlas px).")
                    if q is not None:
                        offset_pix = q
                        offset_wav = offset_pix * refdisp
                else:
                    # sleeps plot_pause (plot_level < 2), no prompt
                    pause(self.config.instrument, self.logger)
                    q = None
            # save final plot of arc/atlas offset
            pfname = "arc_%05d_%s_%s_%s" % (
                self.action.args.ccddata.header['FRAMENO'],
                self.action.args.illum, self.action.args.grating,
                self.action.args.ifuname)
            save_plot(p, filename=pfname + '_atoff.png')
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
