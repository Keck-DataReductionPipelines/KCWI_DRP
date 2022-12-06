from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import strip_fname
from kcwidrp.core.bokeh_plotting import bokeh_plot
import os
import numpy as np
from scipy.interpolate import interpolate
import pickle
from bokeh.plotting import figure
from astropy.io import fits as pf


class SolveAIT(BasePrimitive):
    """Solve the AIT data of the IFU"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.action.args.geometry_file = None
        self.action.args.x0out = None
        self.action.args.wave0out = None
        self.action.args.wave1out = None
        self.action.args.wavegood0 = None
        self.action.args.wavegood1 = None
        self.action.args.waveall0 = None
        self.action.args.waveall1 = None
        self.action.args.wavemid = None
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Writing out AIT spectra")

        log_string = SolveAIT.__module__

        do_plot = (self.config.instrument.plot_level >= 3)

        # Get some geometry constraints
        goody0 = 0
        goody1 = max(self.action.args.xsvals)
        # N&S limits
        goodnsy0 = self.action.args.shufrows
        goodnsy1 = goodnsy0 + self.action.args.shufrows
        # Calculate wavelength ranges
        y0wvs = []
        y1wvs = []
        # N&S ranges
        y0nswvs = []
        y1nswvs = []
        # Get wavelength extremes for each bar
        for fcfs in self.action.args.fincoeff:
            y0wvs.append(float(np.polyval(fcfs, goody0)))
            y1wvs.append(float(np.polyval(fcfs, goody1)))
            y0nswvs.append(float(np.polyval(fcfs, goodnsy0)))
            y1nswvs.append(float(np.polyval(fcfs, goodnsy1)))
        # Now get ensemble extremes
        y0max = max(y0wvs)
        y0min = min(y0wvs)
        y1max = max(y1wvs)
        y1min = min(y1wvs)
        y0nsmax = max(y0nswvs)
        y0nsmin = min(y0nswvs)
        y1nsmax = max(y1nswvs)
        y1nsmin = min(y1nswvs)
        # Cube trimming wavelengths
        trimw0 = y0min
        trimw1 = y1max
        # Check for negative dispersion
        if trimw0 > trimw1:
            trimw0 = y1min
            trimw1 = y0max
        # Calculate output wavelengths
        dwout = self.action.args.dwout
        self.logger.info("Output Delta WAVE: %.3f", dwout)
        ndels = int((trimw0 - self.config.instrument.WAVEFID) / dwout)
        self.action.args.wave0out = \
            self.config.instrument.WAVEFID + float(ndels) * dwout
        ndels = int((trimw1 - self.config.instrument.WAVEFID) / dwout)
        self.action.args.wave1out = \
            self.config.instrument.WAVEFID + float(ndels) * dwout
        self.logger.info("WAVE RANGE: %.2f - %.2f" %
                         (self.action.args.wave0out, self.action.args.wave1out))
        # Calculate wavelength limits
        self.action.args.wavegood0 = min([y0max, y1max])
        self.action.args.wavegood1 = max([y0min, y1min])
        self.action.args.waveall0 = min([y0min, y1min])
        self.action.args.waveall1 = max([y0max, y1max])
        self.action.args.wavemid = np.average([self.action.args.wavegood0,
                                               self.action.args.wavegood1,
                                               self.action.args.waveall0,
                                               self.action.args.waveall1])
        self.action.args.wavensgood0 = min([y0nsmax, y1nsmax])
        self.action.args.wavensgood1 = max([y0nsmin, y1nsmin])
        self.action.args.wavensall0 = min([y0nsmin, y1nsmin])
        self.action.args.wavensall1 = max([y0nsmax, y1nsmax])
        self.action.args.wavensmid = np.average([self.action.args.wavensgood0,
                                                 self.action.args.wavensgood1,
                                                 self.action.args.wavensall0,
                                                 self.action.args.wavensall1])
        self.logger.info("WAVE  GOOD: %.2f - %.2f" %
                         (self.action.args.wavegood0,
                          self.action.args.wavegood1))
        self.logger.info("WAVE   ALL: %.2f - %.2f" %
                         (self.action.args.waveall0, self.action.args.waveall1))
        self.logger.info("WAVE   MID: %.2f" % self.action.args.wavemid)

        # Use extremes to define output size
        ysize = int((self.action.args.waveall1 - self.action.args.wave0out)
                    / dwout)
        yout = np.arange(ysize) * dwout + self.action.args.wave0out
        print(len(yout), np.nanmin(yout), np.nanmax(yout))
        xsize = self.config.instrument.NBARS * 2
        self.logger.info("Output spectra will be %d x %d px" % (xsize, ysize))
        specout = np.zeros((ysize, xsize), dtype=float)
        # Loop over AIT traces
        for isl in range(0, self.config.instrument.NBARS):
            fc = self.action.args.fincoeff[isl]
            arc = self.context.arcs[isl]
            bar = self.context.bars[isl]
            xwv = np.polyval(fc, np.arange(len(arc)))
            arcint = interpolate.interp1d(xwv, arc, kind='cubic',
                                          bounds_error=False,
                                          fill_value='extrapolate')
            barint = interpolate.interp1d(xwv, bar, kind='cubic',
                                          bounds_error=False,
                                          fill_value='extrapolate')
            arcout = arcint(yout)
            barout = barint(yout)
            specout[:, isl] = arcout[:]
            specout[:, isl + self.config.instrument.NBARS] = barout[:]
            if do_plot:
                p = figure(title=self.action.args.plotlabel + "ARC # %d" %
                           (isl + 1),
                           x_axis_label="Wavelength",
                           y_axis_label="Flux",
                           plot_width=self.config.instrument.plot_width,
                           plot_height=self.config.instrument.plot_height)
                p.line(xwv, arc, legend_label='Arc Data', color='blue')
                p.line(xwv, bar, legend_label='Bar Data', color='red')
                p.line(yout, arcout, legend_label='Arc Int', color='purple')
                p.line(yout, barout, legend_label='Bar Int', color='orange')
                bokeh_plot(p, self.context.bokeh_session)
                q = input("Next? <cr>, q to quit: ")
                if 'Q' in q.upper():
                    do_plot = False

        # Pixel scales
        pxscl = self.config.instrument.PIXSCALE * self.action.args.xbinsize
        ifunum = self.action.args.ifunum
        if ifunum == 2:
            slscl = self.config.instrument.SLICESCALE / 2.0
        elif ifunum == 3:
            slscl = self.config.instrument.SLICESCALE / 4.0
        else:
            slscl = self.config.instrument.SLICESCALE

        # Update header
        # Geometry corrected?
        self.action.args.ccddata.header['GEOMCOR'] = (
            True, 'Geometry corrected?')
        #
        # Spatial geometry
        # self.action.args.ccddata.header['BARSEP'] = (
        #    self.action.args.barsep, 'separation of bars (binned pix)')
        # self.action.args.ccddata.header['BAR0'] = (
        #    self.action.args.bar0, 'first bar pixel position')
        # Wavelength ranges
        if self.action.args.nasmask:
            self.action.args.ccddata.header['WAVALL0'] = (
                self.action.args.wavensall0, 'Low inclusive wavelength')
            self.action.args.ccddata.header['WAVALL1'] = (
                self.action.args.wavensall1, 'High inclusive wavelength')
            self.action.args.ccddata.header['WAVGOOD0'] = (
                self.action.args.wavensgood0, 'Low good wavelength')
            self.action.args.ccddata.header['WAVGOOD1'] = (
                self.action.args.wavensgood1, 'High good wavelength')
            self.action.args.ccddata.header['WAVMID'] = (
                self.action.args.wavensmid, 'middle wavelength')
        else:
            self.action.args.ccddata.header['WAVALL0'] = (
                self.action.args.waveall0, 'Low inclusive wavelength')
            self.action.args.ccddata.header['WAVALL1'] = (
                self.action.args.waveall1, 'High inclusive wavelength')
            self.action.args.ccddata.header['WAVGOOD0'] = (
                self.action.args.wavegood0, 'Low good wavelength')
            self.action.args.ccddata.header['WAVGOOD1'] = (
                self.action.args.wavegood1, 'High good wavelength')
            self.action.args.ccddata.header['WAVMID'] = (
                self.action.args.wavemid, 'middle wavelength')
        # Dichroic fraction
        try:
            dichroic_fraction = self.action.args.dich_frac
        except AttributeError:
            dichroic_fraction = 1.
        self.action.args.ccddata.header['DICHFRAC'] = (
            dichroic_fraction, 'Dichroic Fraction')
        # Wavelength fit statistics
        self.action.args.ccddata.header['AVWVSIG'] = (
            self.action.args.av_bar_sig, 'Avg. bar wave sigma (Ang)')
        self.action.args.ccddata.header['SDWVSIG'] = (
            self.action.args.st_bar_sig, 'Stdev. var wave sigma (Ang)')
        # Pixel scales
        self.action.args.ccddata.header['PXSCL'] = (
            pxscl, 'Pixel scale along slice (deg)')
        self.action.args.ccddata.header['SLSCL'] = (
            slscl, 'Pixel scale perp. to slices (deg)')
        # Geometry origins
        self.action.args.ccddata.header['CBARSNO'] = (
            self.action.args.contbar_image_number,
            'Continuum bars image number')
        self.action.args.ccddata.header['CBARSFL'] = (
            self.action.args.contbar_image, 'Continuum bars image filename')
        self.action.args.ccddata.header['ARCNO'] = (
            self.action.args.arc_number, 'Arc image number')
        self.action.args.ccddata.header['ARCFL'] = (
            self.action.args.arc_image, 'Arc image filename')
        # self.action.args.ccddata.header['GEOMFL'] = (
        #     geom_file.split('/')[-1], 'Geometry file')
        # WCS
        # self.action.args.ccddata.header['IFUPA'] = (
        #     skypa, 'IFU position angle (degrees)')
        # self.action.args.ccddata.header['IFUROFF'] = (
        #    self.config.instrument.ROTOFF, 'IFU-SKYPA offset (degrees)')
        self.action.args.ccddata.header['WCSDIM'] = (
            2, 'number of dimensions in WCS')
        self.action.args.ccddata.header['WCSNAME'] = 'KCRM_AIT'
        self.action.args.ccddata.header['CTYPE1'] = 'Pix'
        self.action.args.ccddata.header['CTYPE2'] = ('AWAV',
                                                     'Air Wavelengths')
        self.action.args.ccddata.header['CUNIT1'] = ('pixel', 'image units')
        self.action.args.ccddata.header['CUNIT2'] = ('Angstrom',
                                                     'Wavelength units')
        self.action.args.ccddata.header['CNAME1'] = ('KCWI Pix', 'pix name')
        self.action.args.ccddata.header['CNAME2'] = ('KCWI Wavelength',
                                                     'Wavelength name')
        self.action.args.ccddata.header['CRVAL1'] = (1, 'Pixel zeropoint')
        self.action.args.ccddata.header['CRVAL2'] = (self.action.args.wave0out,
                                                     'Wavelength zeropoint')
        self.action.args.ccddata.header['CDELT1'] = (1.0, 'Pixel scale')
        self.action.args.ccddata.header['CDELT2'] = (dwout, 'Wavelength scale')
        self.action.args.ccddata.header['CRPIX1'] = (1,
                                                     'Pixel reference pixel')
        self.action.args.ccddata.header['CRPIX2'] = (
            1, 'Wavelength reference pixel')

        # write out cube
        self.action.args.ccddata.header['HISTORY'] = log_string

        # Package AIT spectra
        ofname = self.action.args.name
        self.action.args.ait_file = os.path.join(
            self.config.instrument.output_directory,
            strip_fname(ofname) + '_aitspec.fits')
        if os.path.exists(self.action.args.ait_file):
            self.logger.error("AIT spec file already exists: %s" %
                              self.action.args.ait_file)
        else:

            hdu = pf.PrimaryHDU(specout, header=self.action.args.ccddata.header)
            hdul = pf.HDUList([hdu])
            hdul.writeto(self.action.args.ait_file)

            self.logger.info("AIT spectra written to: %s" %
                             self.action.args.ait_file)

        self.logger.info(log_string)

        return self.action.args
    # END: class SolveAIT()
