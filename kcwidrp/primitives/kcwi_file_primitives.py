from keckdrpframework.models.arguments import Arguments
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.table import Table
# from astropy import units as u
import numpy as np
from datetime import datetime

from keckdrpframework.primitives.base_primitive import BasePrimitive
import os
import logging
import pkg_resources
import subprocess
from pathlib import Path

logger = logging.getLogger('KCWI')

red_amp_dict = {'L1': 0, 'L2': 1, 'U1': 2, 'U2': 3}
red_amp_gain = {'L1': 1.54, 'L2': 1.55, 'U1': 1.61, 'U2': 1.49}


def parse_imsec(section=None):

    xsec = section[1:-1].split(',')[0]
    ysec = section[1:-1].split(',')[1]
    xparts = xsec.split(':')
    yparts = ysec.split(':')
    p1 = int(xparts[0])
    p2 = int(xparts[1])
    p3 = int(yparts[0])
    p4 = int(yparts[1])
    # check for scale factor
    if len(xparts) == 3:
        xstride = int(xparts[2])
    else:
        xstride = 1
    if len(yparts) == 3:
        ystride = int(yparts[2])
    else:
        ystride = 1
    # is x axis in descending order?
    if p1 > p2:
        x0 = p2 - 1
        x1 = p1 - 1
        xstride = -abs(xstride)
    # x axis in ascending order
    else:
        x0 = p1 - 1
        x1 = p2 - 1
    # is y axis in descending order
    if p3 > p4:
        y0 = p4 - 1
        y1 = p3 - 1
        ystride = -abs(ystride)
    # y axis in ascending order
    else:
        y0 = p3 - 1
        y1 = p4 - 1
    # ensure no negative indices
    if x0 < 0:
        x0 = 0
    if y0 < 0:
        y0 = 0
    # package output with python axis ordering
    sec = (y0, y1, x0, x1)
    stride = (ystride, xstride)

    return sec, stride


def plotlabel(args):
    lab = "[ Img # %d " % args.ccddata.header['FRAMENO']
    lab += "(%s) " % args.illum
    lab += "Slicer: %s, " % args.ifuname
    lab += "Grating: %s, " % args.grating
    lab += "CWave: %d" % int(args.cwave)
    if 'BLUE' in args.ccddata.header['CAMERA'].upper():
        lab += ", Filter: %s " % args.filter
    lab += "] "
    return lab


class ingest_file(BasePrimitive):

    def __init__(self, action, context):
        """
        Constructor
        """
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def get_keyword(self, keyword):
        try:
            keyval = self.ccddata.header[keyword]
        except KeyError:
            keyval = None
        # return self.context.data_set.get_info_column(self.name, keyword)
        return keyval

    def camera(self):
        camera = self.get_keyword('CAMERA').upper()
        if 'BLUE' in camera:
            return 0
        elif 'RED' in camera:
            return 1
        else:
            return -1

    def camang(self):
        if self.camera() == 0:  # Blue
            key = 'BARTANG'
        elif self.camera() == 1:  # Red
            key = 'RARTANG'
        else:
            raise ValueError("unable to determine camera angle: "
                             "CAMERA undefined")
        return self.get_keyword(key)

    def filter(self):
        if self.camera() == 0:  # Blue
            filt = self.get_keyword('BFILTNAM')
        elif self.camera() == 1:  # Red
            filt = 'None'
        else:
            raise ValueError("unable to determine filter: "
                             "CAMERA undefined")
        return filt

    def grangle(self):
        if self.camera() == 0:  # Blue
            key = 'BGRANGLE'
        elif self.camera() == 1:  # Red
            key = 'RGRANGLE'
        else:
            raise ValueError("unable to determine grating angle: "
                             "CAMERA undefined")
        return self.get_keyword(key)

    def grating(self):
        if self.camera() == 0:  # Blue
            key = 'BGRATNAM'
        elif self.camera() == 1:  # Red
            key = 'RGRATNAM'
        else:
            raise ValueError("unable to determine grating: CAMERA undefined")
        return self.get_keyword(key)

    def adjang(self):
        if 'BH' in self.grating() or 'RH' in self.grating():
            return 180.0
        if 'BM' in self.grating() or 'RM' in self.grating():
            return 0.0
        if 'BL' in self.grating() or 'RL' in self.grating():
            return 0.0

    def atsig(self):
        if 'H' in self.grating():
            if self.ifunum() > 2:
                return 0.75
            else:
                return 2.5
        elif 'M' in self.grating():
            if self.ifunum() >= 2:
                return 4.0
            else:
                return 8.0
        elif 'L' in self.grating():
            if self.ifunum() == 2:
                return 10.0
            elif self.ifunum() == 3:
                return 5.0
            else:
                return 14.0
        else:
            if 'BIAS' in self.imtype():
                return 0.
            else:
                raise ValueError(
                    "unable to determine atlas sigma: GRATING undefined")

    def rho(self):
        if 'BH1' in self.grating():
            return 3.6
        elif 'BH2' in self.grating():
            return 3.2805
        elif 'BH3' in self.grating():
            return 2.8017
        elif 'RH1' in self.grating():
            return 2.360
        elif 'RH2' in self.grating():
            return 2.020
        elif 'RH3' in self.grating():
            return 1.735
        elif 'RH4' in self.grating():
            return 1.490
        elif 'BM' in self.grating():
            return 1.901
        elif 'RM1' in self.grating():
            return 1.300
        elif 'RM2' in self.grating():
            return 0.945
        elif 'BL' in self.grating():
            return 0.870
        elif 'RL' in self.grating():
            return 0.530
        else:
            return None

    def cwave(self):
        if self.camera() == 0:  # Blue
            key = 'BCWAVE'
        elif self.camera() == 1:  # Red
            key = 'RCWAVE'
        else:
            raise ValueError("unable to determine central wavelength: "
                             "CAMERA undefined")
        return self.get_keyword(key)

    def dich(self):
        if self.camera() == 0:  # Blue
            if self.get_keyword('RCWAVE'):
                return True
            else:
                return False
        elif self.camera() == 1:    # Red
            return True

    def resolution(self):
        """Return FWHM resolution in Angstroms for the given grating"""
        # get reference wavelength
        refwave = self.cwave()
        if refwave:
            rw = refwave
        else:
            if 'B' in self.grating():
                rw = 4500.
            else:
                rw = 7500.
        # Calc rez from grating resolution (d-lambda = lambda/R)
        # First, assume large slicer IFU
        if 'BH' in self.grating():
            rez = rw / 5000.
        elif 'RH' in self.grating():
            rez = rw / 4600.
        elif 'BM' in self.grating():
            rez = rw / 2500.
        elif 'RM' in self.grating():
            rez = rw / 1900.
        elif 'BL' in self.grating():
            rez = rw / 1250.
        elif 'RL' in self.grating():
            rez = rw / 800.
        else:
            # Resolution not needed for biases
            if 'BIAS' in self.imtype():
                rez = 0.
            else:
                raise ValueError("unable to compute atlas resolution: "
                                 "grating undefined")
        # Adjust for slicer
        if self.ifunum() == 2:  # Medium slicer
            rez /= 2.
        elif self.ifunum() == 3:  # Small slicer
            rez /= 4.

        return rez

    def delta_wave_out(self):
        """Return output delta lambda in Angstroms for the given grating"""
        # Calc delta wave out from grating
        if 'BH' in self.grating():
            dw = 0.125 * float(self.ybinsize())
        elif 'RH' in self.grating():
            dw = 0.125 * float(self.ybinsize())
        elif 'BM' in self.grating():
            dw = 0.25 * float(self.ybinsize())
        elif 'RM' in self.grating():
            dw = 0.25 * float(self.ybinsize())
        elif 'BL' in self.grating():
            dw = 0.5 * float(self.ybinsize())
        elif 'RL' in self.grating():
            dw = 0.5 * float(self.ybinsize())
        else:
            # Delta wave not needed for biases
            if 'BIAS' in self.imtype():
                dw = 0.
            else:
                raise ValueError("unable to compute output delta lambda: "
                                 "grating undefined")
        return dw

    def namps(self):
        return self.get_keyword('NVIDINP')

    def nasmask(self):
        if self.camera() == 0:  # Blue
            if 'Mask' in self.get_keyword('BNASNAM'):
                return True
            else:
                return False
        elif self.camera() == 1:  # Red
            if 'Mask' in self.get_keyword('RNASNAM'):
                return True
            else:
                return False
        else:
            raise ValueError("unable to determine mask: CAMERA undefined")

    def numopen(self):
        return self.get_keyword('NUMOPEN')

    def shufrows(self):
        return self.get_keyword('SHUFROWS')

    def ampmode(self):
        return self.get_keyword('AMPMODE')

    def xbinsize(self):
        return int(self.get_keyword('BINNING').split(',')[0])

    def ybinsize(self):
        return int(self.get_keyword('BINNING').split(',')[-1])

    def plotlabel(self):
        lab = "[ Img # %d " % self.get_keyword('FRAMENO')
        lab += "(%s) " % self.illum()
        lab += "Slicer: %s, " % self.ifuname()
        lab += "Grating: %s, " % self.grating()
        lab += "CWave: %d, " % int(self.cwave())
        lab += "Filter: %s " % self.filter()
        lab += "] "
        return lab

    def stdlabel(self):
        lab = "[ Img # %d " % self.get_keyword('FRAMENO')
        lab += "(%s) " % self.get_keyword('TARGNAME')
        lab += "Slicer: %s, " % self.ifuname()
        lab += "Grating: %s, " % self.grating()
        lab += "CWave: %d, " % int(self.cwave())
        lab += "Filter: %s " % self.filter()
        lab += "] "
        return lab

    def ifuname(self):
        return self.get_keyword('IFUNAM')

    def ifunum(self):
        return self.get_keyword('IFUNUM')

    def imtype(self):
        return self.get_keyword('IMTYPE')

    def illum(self):
        # set ILLUM keyword
        # ARCS
        if self.get_keyword('IMTYPE') == 'ARCLAMP':
            if self.get_keyword('LMP0STAT') == 1 and \
                    self.get_keyword('LMP0SHST') == 1:
                illum = self.get_keyword('LMP0NAM')
            elif self.get_keyword('LMP1STAT') == 1 and \
                    self.get_keyword('LMP1SHST') == 1:
                illum = self.get_keyword('LMP1NAM')
            else:
                illum = 'Test'
        # Internal FLATS
        elif self.get_keyword('IMTYPE') == 'FLATLAMP':
            if self.get_keyword('LMP3STAT') == 1:
                illum = 'Contin'
            else:
                illum = 'Test'
        # DOMES
        elif self.get_keyword('IMTYPE') == 'DOMEFLAT':
            if self.get_keyword('FLIMAGIN') == 'on' or \
                    self.get_keyword('FLSPECTR') == 'on':
                illum = 'Dome'
            else:
                illum = 'Test'
        # Twilight FLATS
        elif self.get_keyword('IMTYPE') == 'TWIFLAT':
            illum = 'Twilit'
        # BARS
        elif self.get_keyword('IMTYPE') == 'CONTBARS':
            if self.get_keyword('LMP3STAT') == 1:
                illum = 'Contin'
            else:
                illum = 'Test'
        # OBJECT
        elif self.get_keyword('IMTYPE') == 'OBJECT':
            obnam = self.get_keyword('OBJECT')
            if obnam is None:
                obnam = self.get_keyword('TARGNAME')
                if obnam is None:
                    obnam = 'Object'
            # clean up the string
            illum = obnam.replace("/", "_").replace(" ", "").replace(".", "_")
        else:
            illum = 'Test'
        return illum

    def calibration_lamp(self):
        if self.get_keyword('IMTYPE') != 'ARCLAMP':
            return None
        else:
            lamps_dictionary = {
                0: "FeAr",
                1: "ThAr",
                2: "Aux",
                3: "Continuum A"
            }
            for key in lamps_dictionary.keys():
                status = self.get_keyword('LMP%dSTAT' % key)
                shutter = self.get_keyword('LMP%dSHST' % key)
                if status == 1 and shutter == 1:
                    return lamps_dictionary[key]

    def map_ccd(self, xbin, ybin):
        """Return CCD section variables useful for processing

        Uses FITS keyword NVIDINP to determine how many amplifiers were used
        to read out the CCD.  Then reads the corresponding BSECn, and
        DSECn keywords, where n is the amplifier number.  The indices are
        converted to Python (0-biased, y axis first) indices and an array
        is constructed for each of the two useful sections of the CCD as
        follows:

        Bsec[0][0] - First amp, y lower limit
        Bsec[0][1] - First amp, y upper limit
        Bsec[0][2] - First amp, x lower limit
        Bsec[0][3] - First amp, x upper limit
        Bsec[1][0] - Second amp, y lower limit
        etc.

        Bsec is the full overscan region for the given amplifier and is used
        to calculate and perform the overscan subtraction.

        Dsec is the full CCD region for the given amplifier and is used to
        trim the image after overscan subtraction has been performed.

        Tsec accounts for trimming the image according to Dsec.

        Amps are assumed to be organized as follows:
                  BLUE                          RED
        (0,ny)	--------- (nx,ny)    (0,ny)  --------- (nx,ny)
                | 2 | 3 |                    | 0 | 2 |
                ---------                    ---------
                | 0 | 1 |                    | 1 | 3 |
        (0,0)	--------- (nx, 0)     (0,0)  --------- (nx, 0)

        Args:
        -----
            self: instance of CcdPrimitive class (automatic)

        Returns:
        --------
            list: (int) y0, y1, x0, x1 for bias section
            list: (int) y0, y1, x0, x1 for data section
            list: (int) y0, y1, x0, x1 for trimmed section
            list: (bool) y-direction, x-direction, True if forward, else False
        """

        namps = self.namps()    # int(self.get_keyword('NVIDINP'))
        # TODO: check namps
        camera = self.get_keyword('CAMERA').upper()
        ampmode = self.get_keyword('AMPMODE')
        # section lists
        bsec = []
        dsec = []
        tsec = []
        strides = []
        amps = []
        if 'BLUE' in camera:
            nb = 1    # numbering bias (0 or 1)
            # loop over amps
            for i in range(namps):
                ia = i + 1
                amps.append(ia)
                section = self.get_keyword('BSEC%d' % ia)
                sec, stride = parse_imsec(section)
                bsec.append(sec)
                section = self.get_keyword('DSEC%d' % ia)
                sec, stride = parse_imsec(section)
                dsec.append(sec)
                section = self.get_keyword('CSEC%d' % ia)
                sec, stride = parse_imsec(section)
                strides.append(stride)
                if i == 0:
                    y0 = 0
                    y1 = int((sec[1] - sec[0]) / ybin)
                    x0 = 0
                    x1 = int((sec[3] - sec[2]) / xbin)
                elif i == 1:
                    y0 = 0
                    y1 = int((sec[1] - sec[0]) / ybin)
                    x0 = tsec[0][3] + 1
                    x1 = x0 + int((sec[3] - sec[2]) / xbin)
                elif i == 2:
                    y0 = tsec[0][1] + 1
                    y1 = y0 + int((sec[1] - sec[0]) / ybin)
                    x0 = 0
                    x1 = int((sec[3] - sec[2]) / xbin)
                elif i == 3:
                    y0 = tsec[0][1] + 1
                    y1 = y0 + int((sec[1] - sec[0]) / ybin)
                    x0 = tsec[0][3] + 1
                    x1 = x0 + int((sec[3] - sec[2]) / xbin)
                else:
                    # should not get here
                    y0 = -1
                    y1 = -1
                    x0 = -1
                    x1 = -1
                    # self.log.info("ERROR - bad amp number: %d" % i)
                tsec.append((y0, y1, x0, x1))
        elif 'RED' in camera:
            nb = 0    # numbering bias (0 or 1)
            amp_count = 0
            for amp in red_amp_dict.keys():
                if amp in ampmode:
                    amp_count += 1
                    i = red_amp_dict[amp]
                    amps.append(i)
                    section = self.get_keyword('BSEC%d' % i)
                    sec, stride = parse_imsec(section)
                    bsec.append(sec)
                    section = self.get_keyword('DSEC%d' % i)
                    sec, stride = parse_imsec(section)
                    dsec.append(sec)
                    section = self.get_keyword('CSEC%d' % i)
                    sec, stride = parse_imsec(section)
                    strides.append(stride)
                    y0, y1, x0, x1 = sec
                    y0 = int(y0 / abs(stride[0]))
                    y1 = int(y1 / abs(stride[0]))
                    x0 = int(x0 / abs(stride[1]))
                    x1 = int(x1 / abs(stride[1])) - (abs(stride[1]) - 1)
                    tsec.append((y0, y1, x0, x1))
                else:
                    bsec.append((0, 0, 0, 0))
                    dsec.append((0, 0, 0, 0))
                    strides.append((1, 1))
                    tsec.append((0, 0, 0, 0))
            if amp_count != namps:
                self.logger.warning("Didn't get all the amps: %d", amp_count)
        else:
            self.logger.warning("Unknown CAMERA: %s" % camera)
            nb = 0

        return bsec, dsec, tsec, strides, amps, nb

    def _perform(self):
        # if self.context.data_set is None:
        #    self.context.data_set = DataSet(None, self.logger, self.config,
        #    self.context.event_queue)
        # self.context.data_set.append_item(self.action.args.name)
        self.logger.info(
            "------------------- Ingesting file %s -------------------" %
            self.action.args.name)

        self.action.event._recurrent = True

        self.name = self.action.args.name
        out_args = Arguments()

        ccddata, table = kcwi_fits_reader(self.name)

        # save the ccd data into an object
        # that can be shared across the functions
        self.ccddata = ccddata

        out_args.ccddata = ccddata
        out_args.table = table

        imtype = self.get_keyword("IMTYPE")
        groupid = self.get_keyword("GROUPID")

        if groupid is None:
            groupid = "NONE"
        else:
            if len(groupid) <= 0:
                groupid = "NONE"
        if imtype is None:
            fname = os.path.basename(self.action.args.name)
            self.logger.warn(f"Unknown IMTYPE {fname}")
            self.action.event._reccurent = False

        if self.check_if_file_can_be_processed(imtype) is False:
            # self.logger.warn("Object frame cannot be reduced. Rescheduling")
            self.action.new_event = None
            return None
        else:
            self.action.event._recurrent = False

        out_args.name = self.action.args.name
        out_args.imtype = imtype
        out_args.groupid = groupid
        # CAMERA
        out_args.camera = self.camera()
        # DICH
        out_args.dich = self.dich()
        # CAMANGLE
        out_args.camangle = self.camang()
        # FILTER
        out_args.filter = self.filter()
        # GRANGLE
        out_args.grangle = self.grangle()
        # GRATING
        out_args.grating = self.grating()
        # ADJANGLE
        out_args.adjang = self.adjang()
        # RHO
        out_args.rho = self.rho()
        # CWAVE
        out_args.cwave = self.cwave()
        # RESOLUTION
        out_args.resolution = self.resolution()
        # ATSIG
        out_args.atsig = self.atsig()
        # DELTA WAVE OUT
        out_args.dwout = self.delta_wave_out()
        # NAMPS
        out_args.namps = self.namps()
        # NASMASK
        out_args.nasmask = self.nasmask()
        # SHUFROWS
        out_args.shufrows = self.shufrows()
        # NUMOPEN
        out_args.numopen = self.numopen()
        # AMPMODE
        out_args.ampmode = self.ampmode()
        # BINNING
        out_args.xbinsize, out_args.ybinsize = \
            map(int, self.get_keyword('BINNING').split(','))
        # IFUNUM
        out_args.ifunum = int(self.get_keyword('IFUNUM'))
        # IFUNAM
        out_args.ifuname = self.get_keyword('IFUNAM')
        # PLOTLABEL
        out_args.plotlabel = self.plotlabel()
        # STDLABEL
        out_args.stdlabel = self.stdlabel()
        # ILUM
        out_args.illum = self.illum()
        # MAPCCD
        out_args.map_ccd = self.map_ccd(out_args.xbinsize, out_args.ybinsize)
        # CALIBRATION LAMP
        out_args.calibration_lamp = self.calibration_lamp()
        # TTIME
        out_args.ttime = self.get_keyword('TTIME')
        # Are we already in proctab?
        out_args.in_proctab = self.context.proctab.in_proctab(
            frame=out_args.ccddata)

        return out_args

    def apply(self):
        if self._pre_condition():
            try:
                output = self._perform()
            except ValueError as e:
                self.logger.warn("UNABLE TO INGEST THE FILE")
                self.logger.warn("Reason: %s" % e)
                return None
            if self._post_condition():
                self.output = output
        return self.output

    def check_if_file_can_be_processed(self, imtype):

        if imtype == 'OBJECT':
            # bias frames
            bias_frames = self.context.proctab.search_proctab(
                frame=self.ccddata, target_type='MBIAS', nearest=True)
            # master flats
            masterflat_frames = self.context.proctab.search_proctab(
                frame=self.ccddata, target_type='MFLAT', nearest=True)
            # arclamp
            arclamp_frames = self.context.proctab.search_proctab(
                frame=self.ccddata, target_type='MARC', nearest=True)
            if len(bias_frames) > 0:
                '''\
                and len(masterflat_frames) > 0 \
                and len(arclamp_frames) > 0:'''
                return True
            else:
                self.logger.warn("Cannot reduce OBJECT frame. "
                                 "Rescheduling for later. Found:")
                self.logger.warn(f"\tMBIAS: {len(bias_frames)}")
                self.logger.warn(f"\tMFLAT: {len(masterflat_frames)}")
                self.logger.warn(f"\tARCLAMP: {len(arclamp_frames)}")
                return False

        if imtype == 'ARCLAMP':
            # continuum bars
            contbars_frames = self.context.proctab.search_proctab(
                frame=self.ccddata, target_type='MCBARS', nearest=True)
            if (len(contbars_frames)) > 0:
                return True
            else:
                self.logger.warn("Cannot reduce ARCLAMP frame. "
                                 "Missing master continuum bars. "
                                 "Rescheduling for later.")
                return False

        if imtype in ['FLATLAMP', 'TWIFLAT', 'DOMEFLAT']:
            # bias frames
            bias_frames = self.context.proctab.search_proctab(
                frame=self.ccddata, target_type='MBIAS', nearest=True)
            if len(bias_frames) > 0:
                return True
            else:
                self.logger.warn(f"Cannot reduce {imtype} frame. Missing "
                                 f"master bias. Rescheduling for later.")
                return False

        return True


def kcwi_fits_reader(file):
    """A reader for KeckData objects.
    Currently, this is a separate function, but should probably be
    registered as a reader similar to fits_ccddata_reader.
    Arguments:
    file -- The filename (or pathlib.Path) of the FITS file to open.
    """
    try:
        hdul = fits.open(file)
    except (FileNotFoundError, OSError) as e:
        print(e)
        raise e
    read_imgs = 0
    read_tabs = 0
    # primary image
    ccddata = CCDData(hdul['PRIMARY'].data, meta=hdul['PRIMARY'].header,
                      unit='adu')
    read_imgs += 1
    # check for other legal components
    if 'UNCERT' in hdul:
        ccddata.uncertainty = hdul['UNCERT'].data
        read_imgs += 1
    if 'FLAGS' in hdul:
        ccddata.flags = hdul['FLAGS'].data
        read_imgs += 1
    if 'MASK' in hdul:
        ccddata.mask = hdul['MASK'].data
        read_imgs += 1
    if 'Exposure Events' in hdul:
        table = hdul['Exposure Events']
        read_tabs += 1
    else:
        table = None
    # prepare for floating point
    ccddata.data = ccddata.data.astype(np.float64)
    # Fix red headers
    fix_red_header(ccddata)
    # Check for CCDCFG keyword
    if 'CCDCFG' not in ccddata.header:
        ccdcfg = ccddata.header['CCDSUM'].replace(" ", "")
        ccdcfg += "%1d" % ccddata.header['CCDMODE']
        ccdcfg += "%02d" % ccddata.header['GAINMUL']
        ccdcfg += "%02d" % ccddata.header['AMPMNUM']
        ccddata.header['CCDCFG'] = ccdcfg

    if ccddata:
        if 'BUNIT' in ccddata.header:
            ccddata.unit = ccddata.header['BUNIT']
            if ccddata.uncertainty:
                ccddata.uncertainty.unit = ccddata.header['BUNIT']
            # print("setting image units to " + ccddata.header['BUNIT'])

    logger.info("<<< read %d imgs and %d tables out of %d hdus in %s" %
                (read_imgs, read_tabs, len(hdul), file))
    return ccddata, table


def write_table(output_dir=None, table=None, names=None, comment=None,
                keywords=None, output_name=None, clobber=False):
    output_file = os.path.join(output_dir, output_name)
    # check if file already exists
    if os.path.exists(output_file) and clobber is False:
        logger.warning("Table %s already exists and overwrite is set to False"
                       % output_file)
        return
    if os.path.exists(output_file) and clobber is True:
        logger.warning("Table %s already exists and will be overwritten"
                       % output_file)
        logger.info("Removing file: %s" % output_file)
        os.remove(output_file)

    t = Table(table, names=names)
    if comment:
        t.meta['COMMENT'] = comment
    if keywords:
        for k, v in keywords.items():
            t.meta[k] = v
    try:
        t.write(output_file, format='fits')
        logger.info("Output table: %s" % output_file)
    except (FileExistsError, OSError):
        logger.error("Something went wrong creating the table: "
                     "table already exists")


def read_table(input_dir=None, file_name=None):
    # Set up return table
    input_file = os.path.join(input_dir, file_name)
    logger.info("Trying to read table: %s" % input_file)
    try:
        retab = Table.read(input_file, format='fits')
    except FileNotFoundError:
        logger.warning("No table to read")
        retab = None
    return retab


def kcwi_fits_writer(ccddata, table=None, output_file=None, output_dir=None,
                     suffix=None):
    
    # Determine if the version info is already in the header
    contains_version = False
    for h in ccddata.header["HISTORY"]:
        if "kcwidrp version" in h:
            contains_version = True

    if not contains_version:
        # Add setup.py version number to header
        version = pkg_resources.get_distribution('kcwidrp').version
        ccddata.header.add_history(f"kcwidrp version={version}")

        # Get string filepath to .git dir, relative to this primitive
        primitive_loc = os.path.dirname(os.path.abspath(__file__))
        git_loc = primitive_loc[:-18] + ".git"

        # Attempt to gather git version information
        git1 = subprocess.run(["git", "--git-dir", git_loc, "describe",
                               "--tags", "--long"], capture_output=True)
        git2 = subprocess.run(["git", "--git-dir", git_loc, "log", "-1",
                               "--format=%cd"], capture_output=True)
        
        # If all went well, save to the header
        if not bool(git1.stderr) and not bool(git2.stderr):
            git_v = git1.stdout.decode('utf-8')[:-1]
            git_d = git2.stdout.decode('utf-8')[:-1]
            ccddata.header.add_history(f"git version={git_v}")
            ccddata.header.add_history(f"git date={git_d}")
        else:
            logger.debug("Package not installed from a git repo, skipping")

    # If there is a data array, and the type of that array is a 64-bit float,
    # force it to 32 bits.
    if ccddata.data is not None and ccddata.data.dtype == np.float64:
        ccddata.data = ccddata.data.astype(np.float32)
    # If there is an uncertainty array, and the values within (the .array property), make it
    # 32 bits.
    if ccddata.uncertainty is not None and ccddata.uncertainty.array.dtype == np.float64:
        ccddata.uncertainty.array = ccddata.uncertainty.array.astype(np.float32)
    
    out_file = os.path.join(output_dir, os.path.basename(output_file))
    if suffix is not None:
        (main_name, extension) = os.path.splitext(out_file)
        out_file = main_name + "_" + suffix + extension
    hdus_to_save = ccddata.to_hdu()
    # check for flags
    flags = getattr(ccddata, "flags", None)
    if flags is not None:
        hdus_to_save.append(fits.ImageHDU(flags, name='FLAGS',
                                          do_not_scale_image_data=True))
    # something about the way the original table is written out is wrong
    # and causes problems.  Leaving it off for now.
    # if table is not None:
    #    hdus_to_save.append(table)
    # log
    logger.info(">>> Saving %d hdus to %s" % (len(hdus_to_save), out_file))
    hdus_to_save.writeto(out_file, overwrite=True)


def strip_fname(filename):
    if not filename:
        logger.error(f"Failed to strip file {filename}")
        return
    strip = Path(filename).stem
    return strip


def get_master_name(tab, target_type, loc=0):
    res = Path(strip_fname(tab['filename'][loc]) + '_' +
               target_type.lower() + ".fits").name
    return res


def master_bias_name(ccddata, target_type='MBIAS'):
    # Delivers a mbias filename that is unique for each CCD configuration
    # Any KCWI frame with a shared CCD configuration can use the same bias
    name = target_type.lower() + '_' + ccddata.header['CCDCFG'] + '.fits'
    return name


def master_flat_name(ccddata, target_type):
    # Delivers a name that is unique across an observing block
    name = target_type.lower() + '_' + ccddata.header['STATEID'] + '.fits'
    return name


def fix_red_header(ccddata):
    # are we red?
    if 'RED' in ccddata.header['CAMERA'].upper():
        # Fix red headers during Caltech AIT
        if 'TELESCOP' not in ccddata.header:
            # Add DCS keywords
            ccddata.header['TARGNAME'] = ccddata.header['OBJECT']
            dateend = ccddata.header['DATE-END']
            de = datetime.fromisoformat(dateend)
            day_frac = de.hour / 24. + de.minute / 1440. + de.second / 86400.
            jd = datetime.date(de).toordinal() + day_frac + 1721424.5
            mjd = jd - 2400000.5
            ccddata.header['MJD'] = mjd
        # Add NVIDINP
        ccddata.header['NVIDINP'] = ccddata.header['TAPLINES']
        # Add GAINMUL
        ccddata.header['GAINMUL'] = 1
        # Add CCDMODE
        if 'CDSSPEED' in ccddata.header:
            ccddata.header['CCDMODE'] = ccddata.header['CDSSPEED']
        else:
            ccddata.header['CCDMODE'] = 0
        # Add AMPMNUM
        if 'AMPMODE' in ccddata.header:
            ampmode = ccddata.header['AMPMODE']
            if ampmode == 'L2U2L1U1':
                ampnum = 0
            elif ampmode == 'L2U2':
                ampnum = 1
            elif ampmode == 'L1U1':
                ampnum = 2
            elif ampmode == 'L2L1':
                ampnum = 3
            elif ampmode == 'U2U1':
                ampnum = 4
            elif ampmode == 'L2':
                ampnum = 5
            elif ampmode == 'U2':
                ampnum = 6
            elif ampmode == 'L1':
                ampnum = 7
            elif ampmode == 'U1':
                ampnum = 8
            else:
                ampnum = 0
            ccddata.header['AMPMNUM'] = ampnum

            for amp in red_amp_dict.keys():
                if amp in ampmode:
                    gkey = 'GAIN%d' % red_amp_dict[amp]
                    gain = red_amp_gain[amp]
                    ccddata.header[gkey] = gain
