from keckdrpframework.models.arguments import Arguments
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.table import Table
import numpy as np

from keckdrpframework.primitives.base_primitive import BasePrimitive
import os
import logging

logger = logging.getLogger('KCWI')


def parse_imsec(section=None):

    xfor = True
    yfor = True

    p1 = int(section[1:-1].split(',')[0].split(':')[0])
    p2 = int(section[1:-1].split(',')[0].split(':')[1])
    p3 = int(section[1:-1].split(',')[1].split(':')[0])
    p4 = int(section[1:-1].split(',')[1].split(':')[1])
    # tests for individual axes
    if p1 > p2:
        x0 = p2 - 1
        x1 = p1 - 1
        xfor = False
    else:
        x0 = p1 - 1
        x1 = p2 - 1
    if p3 > p4:
        y0 = p4 - 1
        y1 = p3 - 1
        yfor = False
    else:
        y0 = p3 - 1
        y1 = p4 - 1
    # package output
    sec = (y0, y1, x0, x1)
    rfor = (yfor, xfor)
    # use python axis ordering
    return sec, rfor


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
        camera = self.get_keyword('CAMERA')
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

    def rho(self):
        if 'BH1' in self.grating():
            return 3.751
        elif 'BH2' in self.grating():
            return 3.255
        elif 'BH3' in self.grating():
            return 2.800
        elif 'RH1' in self.grating():
            return 2.420
        elif 'RH2' in self.grating():
            return 2.030
        elif 'RH3' in self.grating():
            return 1.705
        elif 'RH4' in self.grating():
            return 1.435
        elif 'BM' in self.grating():
            return 1.900
        elif 'RM1' in self.grating():
            return 1.220
        elif 'RM2' in self.grating():
            return 0.921
        elif 'BL' in self.grating():
            return 0.870
        elif 'RL' in self.grating():
            return 0.514
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
            rez = rw / 5000.
        elif 'BM' in self.grating():
            rez = rw / 2500.
        elif 'RM' in self.grating():
            rez = rw / 2500.
        elif 'BL' in self.grating():
            rez = rw / 1250.
        elif 'RL' in self.grating():
            rez = rw / 1250.
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

    def shufrows(self):
        return self.get_keyword('SHUFROWS')

    def ampmode(self):
        return self.get_keyword('AMPMODE')

    def xbinsize(self):
        return int(self.get_keyword('BINNING').split(',')[0])

    def ybinsize(self):
        return int(self.get_keyword('BINNING').split(',')[-1])

    def plotlabel(self):
        lab = "[Img # %d " % self.get_keyword('FRAMENO')
        lab += "(%s) " % self.illum()
        lab += "Slicer: %s " % self.ifuname()
        lab += "Filter: %s " % self.filter()
        lab += "Grating: %s" % self.grating()
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
            illum = 'Object'
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

    def map_ccd(self):
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

        (0,ny)	--------- (nx,ny)
                | 3 | 4 |
                ---------
                | 1 | 2 |
        (0,0)	--------- (nx, 0)

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

        namps = int(self.get_keyword('NVIDINP'))
        # TODO: check namps
        # section lists
        bsec = []
        dsec = []
        tsec = []
        direc = []
        # loop over amps
        for i in range(namps):
            section = self.get_keyword('BSEC%d' % (i + 1))
            sec, rfor = parse_imsec(section)
            bsec.append(sec)
            section = self.get_keyword('DSEC%d' % (i + 1))
            sec, rfor = parse_imsec(section)
            dsec.append(sec)
            direc.append(rfor)
            if i == 0:
                y0 = 0
                y1 = sec[1] - sec[0]
                x0 = 0
                x1 = sec[3] - sec[2]
            elif i == 1:
                y0 = 0
                y1 = sec[1] - sec[0]
                x0 = tsec[0][3] + 1
                x1 = x0 + sec[3] - sec[2]
            elif i == 2:
                y0 = tsec[0][1] + 1
                y1 = y0 + sec[1] - sec[0]
                x0 = 0
                x1 = sec[3] - sec[2]
            elif i == 3:
                y0 = tsec[0][1] + 1
                y1 = y0 + sec[1] - sec[0]
                x0 = tsec[0][3] + 1
                x1 = x0 + sec[3] - sec[2]
            else:
                # should not get here
                y0 = -1
                y1 = -1
                x0 = -1
                x1 = -1
                # self.log.info("ERROR - bad amp number: %d" % i)
            tsec.append((y0, y1, x0, x1))

        return bsec, dsec, tsec, direc

    def _perform(self):
        # if self.context.data_set is None:
        #    self.context.data_set = DataSet(None, self.logger, self.config,
        #    self.context.event_queue)
        # self.context.data_set.append_item(self.action.args.name)
        self.logger.info(
            "------------------- Ingesting file %s -------------------" %
            self.action.args.name)
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
        if imtype is None:
            fname = os.path.basename(self.action.args.name)
            self.logger.warn(f"Unknown IMTYPE {fname}")

        out_args.name = self.action.args.name
        out_args.imtype = imtype
        out_args.groupid = groupid
        # CAMERA
        out_args.camera = self.camera()
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
        # DELTA WAVE OUT
        out_args.dwout = self.delta_wave_out()
        # NAMPS
        out_args.namps = int(self.get_keyword('NVIDINP'))
        # NASMASK
        out_args.nasmask = self.nasmask()
        # SHUFROWS
        out_args.shufrows = self.shufrows()
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
        # ILUM
        out_args.illum = self.illum()
        # MAPCCD
        out_args.map_ccd = self.map_ccd()
        # CALIBRATION LAMP
        out_args.calibration_lamp = self.calibration_lamp()
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


def kcwi_fits_reader(file):
    """A reader for KeckData objects.
    Currently this is a separate function, but should probably be
    registered as a reader similar to fits_ccddata_reader.
    Arguments:
    file -- The filename (or pathlib.Path) of the FITS file to open.
    """
    try:
        hdul = fits.open(file)
    except FileNotFoundError as e:
        print(e)
        raise e
    except OSError as e:
        print(e)
        raise e

    if len(hdul) == 1:
        # Simple input, no table
        ccddata = CCDData(hdul[0].data, meta=hdul[0].header, unit='adu')
        table = None
    elif len(hdul) == 2:
        # raw data
        # 1- read the first extension into a ccddata
        ccddata = CCDData(hdul[0].data, meta=hdul[0].header, unit='adu')
        # 2- read the table
        table = hdul[1]
        # 3- prepare for floating point
        ccddata.data = ccddata.data.astype(np.float64)
        # 4- Check for CCDCFG keyword
        if 'CCDCFG' not in ccddata.header:
            ccdcfg = ccddata.header['CCDSUM'].replace(" ", "")
            ccdcfg += "%1d" % ccddata.header['CCDMODE']
            ccdcfg += "%02d" % ccddata.header['GAINMUL']
            ccdcfg += "%02d" % ccddata.header['AMPMNUM']
            ccddata.header['CCDCFG'] = ccdcfg
    elif len(hdul) == 3:
        # pure ccd data
        ccddata = CCDData(hdul['PRIMARY'].data, meta=hdul['PRIMARY'].header,
                          unit='adu')
        ccddata.mask = hdul['MASK'].data
        ccddata.uncertainty = hdul['UNCERT'].data
        table = None
    elif len(hdul) == 4:
        # pure ccd data
        ccddata = CCDData(hdul['PRIMARY'].data, meta=hdul['PRIMARY'].header,
                          unit='adu')
        ccddata.mask = hdul['MASK'].data
        ccddata.uncertainty = hdul['UNCERT'].data
        table = Table(hdul[3])
    else:
        logger.warning("Wrong number of HDUnits in %s: should be 1-4, but is %d"
                       % (file, len(hdul)))
        ccddata = None
        table = None

    if ccddata:
        if 'BUNIT' in ccddata.header:
            ccddata.unit = ccddata.header['BUNIT']
            # print("setting image units to " + ccddata.header['BUNIT'])

    return ccddata, table


class kcwi_fits_ingest(BasePrimitive):
    """
    classdocs
    """

    def __init__(self, action, context):
        """
        Constructor
        """
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        """
        Expects action.args.name as fits file name
        Returns HDUs or (later) data model
        """
        name = self.action.args.name
        self.logger.info(f"Reading {name}")
        out_args = Arguments()
        out_args.name = name
        ccddata, table = kcwi_fits_reader(name)
        out_args.ccddata = ccddata
        out_args.table = table
        out_args.imtype = out_args.hdus.header['IMTYPE']

        return out_args


def write_table(output_dir=None, table=None, names=None, comment=None,
                keywords=None, output_name=None):
    output_file = os.path.join(output_dir, output_name)

    t = Table(table, names=names)
    if comment:
        t.meta['COMMENT'] = comment
    if keywords:
        for k, v in keywords.items():
            t.meta[k] = v
    try:
        t.write(output_file, format='fits')
    except FileExistsError:
        logger.warning("Table already exists")
    print("output file: %s" % output_file)


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


def kcwi_fits_writer(ccddata, table=None, output_file=None, suffix=None):
    output_file = os.path.join(os.path.dirname(output_file), 'redux',
                               os.path.basename(output_file))
    if suffix is not None:
        (main_name, extension) = os.path.splitext(output_file)
        output_file = main_name + "_" + suffix + extension
    hdus_to_save = ccddata.to_hdu()
    # if table is not None:
    #    hdus_to_save.append(table)
    # hdus_to_save.info()
    logger.info(">>> Saving to %s" % output_file)
    hdus_to_save.writeto(output_file, overwrite=True)
