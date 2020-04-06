from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments

import os
import numpy as np
from skimage import transform as tf
import pickle


class SolveGeom(BasePrimitive):
    """Solve the overall geometry of the IFU"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.action.args.geom_file = None
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
        self.logger.info("Solving overall geometry")

        # Get some geometry constraints
        if self.action.args.nasmask:
            goody0 = self.action.args.shufrows + 1
            goody1 = goody0 + self.action.args.shufrows
        else:
            goody0 = 0
            goody1 = max(self.action.args.xsvals)
        # Calculate wavelength ranges
        y0wvs = []
        y1wvs = []
        # Get wavelength extremes for each bar
        for fcfs in self.action.args.fincoeff:
            y0wvs.append(float(np.polyval(fcfs, goody0)))
            y1wvs.append(float(np.polyval(fcfs, goody1)))
        # Now get ensemble extremes
        y0max = max(y0wvs)
        y0min = min(y0wvs)
        y1max = max(y1wvs)
        y1min = min(y1wvs)
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
        self.logger.info("WAVE  GOOD: %.2f - %.2f" %
                         (self.action.args.wavegood0,
                          self.action.args.wavegood1))
        self.logger.info("WAVE   ALL: %.2f - %.2f" %
                         (self.action.args.waveall0, self.action.args.waveall1))
        self.logger.info("WAVE   MID: %.2f" % self.action.args.wavemid)
        # Start setting up slice transforms
        self.action.args.x0out = int(self.action.args.refdelx / 2.) + 1
        self.refoutx = np.arange(0, 5) * self.action.args.refdelx + \
            self.action.args.x0out
        # Variables for output control points
        srcw = []
        max_srcw = 0
        min_srcw = 4096 / self.action.args.ybinsize
        # Loop over source control points
        for ixy, xy in enumerate(self.action.args.src):
            # Calculate y wavelength
            yw = float(np.polyval(
                self.action.args.fincoeff[self.action.args.barid[ixy]], xy[1]))
            # Convert to output pixels
            yw = (yw - self.action.args.wave0out) / dwout
            # Calculate extreme values
            if yw > max_srcw:
                max_srcw = yw
            if yw < min_srcw:
                min_srcw = yw
            srcw.append([xy[0], yw])
        # Use extremes to define output size
        ysize = int(max_srcw + min_srcw + 20 / self.action.args.ybinsize)
        xsize = int(5. * self.action.args.refdelx) + 1
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
                if self.action.args.slid[ixy] == isl:
                    # Index in to reference output x array
                    ib = self.action.args.barid[ixy] % 5
                    # Geometrically corrected control points
                    xw.append(self.refoutx[ib])
                    yw.append(xy[1])
                    # Input control points
                    xi.append(self.action.args.dst[ixy][0])
                    yi.append(self.action.args.dst[ixy][1])
            # get image limits
            xl0 = int(min(xi) - self.action.args.refdelx)
            if xl0 < 0:
                xl0 = 0
            xl1 = int(max(xi) + self.action.args.refdelx)
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
            tform = tf.estimate_transform('polynomial', src, dst, order=3)
            invtf = tf.estimate_transform('polynomial', dst, src, order=3)
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
        # Package geometry data
        ofname = self.action.args.ccddata.header['OFNAME']
        self.action.args.geom_file = os.path.join(
            self.config.instrument.output_directory,
            ofname.split('.')[0] + '_geom.pkl')
        if os.path.exists(self.action.args.geom_file):
            self.logger.error("Geometry file already exists: %s" %
                              self.action.args.geom_file)
        else:
            geom = {
                "geom_file": self.action.args.geom_file,
                "xsize": xsize, "ysize": ysize,
                "pxscl": pxscl, "slscl": slscl,
                "cbarsno": self.action.args.cbarsno,
                "cbarsfl": self.action.args.cbarsfl,
                "arcno": self.action.args.arcno,
                "arcfl": self.action.args.arcfl,
                "barsep": self.action.args.refdelx,
                "bar0": self.action.args.x0out,
                "waveall0": self.action.args.waveall0,
                "waveall1": self.action.args.waveall1,
                "wavegood0": self.action.args.wavegood0,
                "wavegood1": self.action.args.wavegood1,
                "wavemid": self.action.args.wavemid,
                "dwout": dwout,
                "wave0out": self.action.args.wave0out,
                "wave1out": self.action.args.wave1out,
                "avwvsig": self.action.args.av_bar_sig,
                "sdwvsig": self.action.args.st_bar_sig,
                "xl0": xl0_out, "xl1": xl1_out,
                "tform": tform_list, "invtf": invtf_list
            }
            with open(self.action.args.geom_file, 'wb') as ofile:
                pickle.dump(geom, ofile)
            self.logger.info("Geometry written to: %s" %
                             self.action.args.geom_file)

        logstr = SolveGeom.__module__ + "." + SolveGeom.__qualname__
        self.action.args.ccddata.header['HISTORY'] = logstr
        self.logger.info(logstr)

        return self.action.args
    # END: class SolveGeom()


