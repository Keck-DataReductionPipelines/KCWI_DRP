from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import strip_fname
import os
import numpy as np
from kcwidrp.core import geometric as tf
import pickle


class SolveGeom(BasePrimitive):
    """
    Solve the overall geometry of the IFU.

    Takes individual bar arc spectra and generates a fit for each slice, which
    contains five bar arc spectra each.  Given that there are only five points
    in the 'x' or spatial direction and many more in the 'y' or wavelength
    direction, an asymmetric polynomial is fit with order 2 in the x and order
    4 in the y directions.  The wavelength coverage of the observation is
    recorded in parameters with the range that includes all data being in
    the waveall0 and waveall1 parameters, and the range that includes only
    good data in wavegood0 and wavegood1 parameters.  The middle of the
    wavelength range is recorded in the wavemid parameter.

    Forward and inverse transforms, along with all the parameters for the
    geometric fit are written out as a python pickled dictionary in a
    \*_geom.pkl file.

    """

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

    def _pre_condition(self):
        self.logger.info("Checking for master arc")
        if 'MARC' in self.action.args.ccddata.header['IMTYPE']:
            return True
        else:
            return False

    def _perform(self):
        self.logger.info("Solving overall geometry")

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
        # Start setting up slice transforms
        self.action.args.x0out = \
            int(self.action.args.reference_bar_separation / 2.) + 1
        self.refoutx = np.arange(0, 5) * \
            self.action.args.reference_bar_separation + self.action.args.x0out
        # Variables for output control points
        srcw = []
        # Loop over source control points
        for ixy, xy in enumerate(self.action.args.source_control_points):
            # Calculate y wavelength
            yw = float(np.polyval(
                self.action.args.fincoeff[self.action.args.bar_id[ixy]], xy[1]))
            # Convert to output pixels
            yw = (yw - self.action.args.wave0out) / dwout
            srcw.append([xy[0], yw])
        # Use extremes to define output size
        ysize = int((self.action.args.waveall1 - self.action.args.wave0out)
                    / dwout)
        xsize = int(5. * self.action.args.reference_bar_separation) + 1
        self.logger.info("Output slices will be %d x %d px" % (xsize, ysize))
        # Now loop over slices and get relevant control points for each slice
        # Output variables
        xl0_out = []
        xl1_out = []
        tform_list = []
        invtf_list = []
        # Loop over 24 slices
        for isl in range(0, 24):
            # Get control points
            xw = []
            yw = []
            xi = []
            yi = []
            # Loop over all control points
            for ixy, xy in enumerate(srcw):
                # Only use the ones for this slice
                if self.action.args.slice_id[ixy] == isl:
                    # Index in to reference output x array
                    ib = self.action.args.bar_id[ixy] % 5
                    # Geometrically corrected control points
                    xw.append(self.refoutx[ib])
                    yw.append(xy[1])
                    # Input control points
                    xi.append(self.action.args.destination_control_points[
                                  ixy][0])
                    yi.append(self.action.args.destination_control_points[
                                  ixy][1])
            # get image limits
            xl0 = int(min(xi) - self.action.args.reference_bar_separation)
            if xl0 < 0:
                xl0 = 0
            xl1 = int(max(xi) + self.action.args.reference_bar_separation)
            if xl1 > (self.action.args.ccddata.data.shape[0] - 1):
                xl1 = self.action.args.ccddata.data.shape[0] - 1
            # Store for output
            xl0_out.append(xl0)
            xl1_out.append(xl1)
            self.logger.info("Slice %d arc image x limits: %d - %d" %
                             (isl, xl0, xl1))
            # adjust control points
            xit = [x - float(xl0) for x in xi]
            # fit transform
            dst = np.column_stack((xit, yi))
            src = np.column_stack((xw, yw))
            self.logger.info("Fitting wavelength and spatial control points")
            tform = tf.estimate_transform('asympolynomial', src, dst,
                                          order=(2, 4))
            invtf = tf.estimate_transform('asympolynomial', dst, src,
                                          order=(2, 4))
            # Store for output
            tform_list.append(tform)
            invtf_list.append(invtf)
        # Pixel scales
        pxscl = self.config.instrument.PIXSCALE * self.action.args.xbinsize
        ifunum = self.action.args.ifunum
        if ifunum == 2:
            slscl = self.config.instrument.SLICESCALE / 2.0
        elif ifunum == 3:
            slscl = self.config.instrument.SLICESCALE / 4.0
        else:
            slscl = self.config.instrument.SLICESCALE
        # Dichroic fraction
        try:
            dichroic_fraction = self.action.args.dichroic_fraction
        except AttributeError:
            dichroic_fraction = 1.

        # Package geometry data
        # name = self.action.args.ccddata.header['OFNAME']
        name = self.action.args.name
        self.action.args.geometry_file = os.path.join(
            self.config.instrument.output_directory,
            strip_fname(name) + '_geom.pkl')
        if os.path.exists(self.action.args.geometry_file):
            self.logger.error("Geometry file already exists: %s" %
                              self.action.args.geometry_file)
        else:
            geom = {
                "geom_file": self.action.args.geometry_file,            # output save file for geom struct
                "xsize": xsize, "ysize": ysize,                         # xy size of image
                "pxscl": pxscl, "slscl": slscl,                         # degrees per unbinned spatial and slice pixel
                "cbarsno": self.action.args.contbar_image_number,       # cont bars image number
                "cbarsfl": self.action.args.contbar_image,              # cont bars image name
                "arcno": self.action.args.arc_number,                   # arc image number
                "arcfl": self.action.args.arc_image,                    # arc image name
                "barsep": self.action.args.reference_bar_separation,    # reference delta x between bars (pix)
                "bar0": self.action.args.x0out,                         # output spatial zeropoint (unbinned pix)
                "waveall0": self.action.args.waveall0,                  # low wavelength that includes all data
                "waveall1": self.action.args.waveall1,                  # high wavelength that includes all data
                "wavegood0": self.action.args.wavegood0,                # low wavelength that includes good data
                "wavegood1": self.action.args.wavegood1,                # high wavelength that includes good data
                "wavemid": self.action.args.wavemid,                    # wavelength of the middle of the range
                "wavensall0": self.action.args.wavensall0,              # low wavelength that includes all masked data
                "wavensall1": self.action.args.wavensall1,              # high wavelength that includes all masked data
                "wavensgood0": self.action.args.wavensgood0,            # low wavelength that includes good masked data
                "wavensgood1": self.action.args.wavensgood1,            # high wavelength that includes good masked data
                "wavensmid": self.action.args.wavensmid,                # wavelength of the middle of the masked range
                "dich_frac": dichroic_fraction,                         # fraction of wavelength range impacted by dichroic
                "dwout": dwout,                                         # output disperson (Ang/pix)
                "wave0out": self.action.args.wave0out,                  # output wavelength zeropoint
                "wave1out": self.action.args.wave1out,                  # output ending wavelength
                "avwvsig": self.action.args.av_bar_sig,                 # average bar wavelength sigma (Ang)
                "sdwvsig": self.action.args.st_bar_sig,                 # standard deviation of bar sigmas (Ang)
                "xl0": xl0_out, "xl1": xl1_out,                         # slice spatial limits
                "tform": tform_list, "invtf": invtf_list                # transform and inverse transforms for each slice
            }
            with open(self.action.args.geometry_file, 'wb') as ofile:
                pickle.dump(geom, ofile)
            self.logger.info("Geometry written to: %s" %
                             self.action.args.geometry_file)

        log_string = SolveGeom.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        return self.action.args
    # END: class SolveGeom()
