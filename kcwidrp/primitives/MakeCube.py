from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer, \
    kcwi_fits_reader

import time
import os
import math
import pickle
import numpy as np
from skimage import transform as tf
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.nddata import CCDData
from kcwidrp.core.bokeh_plotting import bokeh_plot
from bokeh.plotting import figure
from multiprocessing import Pool


def make_cube_helper(argument):
    """Warp each slice"""
    logger = argument['logger']
    logger.info("Transforming image slice %d" % (argument['slice_number']+1))
    # input params
    slice_number = argument['slice_number']
    xsize = argument['xsize']
    ysize = argument['ysize']
    tform = argument['geom']['tform'][slice_number]
    xl0 = argument['geom']['xl0'][slice_number]
    xl1 = argument['geom']['xl1'][slice_number]
    # slice data
    slice_img = argument['img'][:, xl0:xl1]
    slice_var = argument['std'][:, xl0:xl1]
    slice_msk = argument['msk'][:, xl0:xl1]
    slice_flg = argument['flg'][:, xl0:xl1]
    if 'obj' in argument:
        slice_obj = argument['obj'][:, xl0:xl1]
    else:
        slice_obj = None
    if 'sky' in argument:
        slice_sky = argument['sky'][:, xl0:xl1]
    else:
        slice_sky = None
    if 'del' in argument:
        slice_del = argument['del'][:, xl0:xl1]
    else:
        slice_del = None
    # do the warping
    warped = tf.warp(slice_img, tform, order=3,
                     output_shape=(ysize, xsize))
    varped = tf.warp(slice_var, tform, order=3,
                     output_shape=(ysize, xsize))
    marped = tf.warp(slice_msk, tform, order=3,
                     output_shape=(ysize, xsize))
    farped = tf.warp(slice_flg, tform, order=3,
                     output_shape=(ysize, xsize), preserve_range=True)

    if slice_obj is not None:
        oarped = tf.warp(slice_obj, tform, order=3,
                         output_shape=(ysize, xsize))
    else:
        oarped = None

    if slice_sky is not None:
        sarped = tf.warp(slice_sky, tform, order=3,
                         output_shape=(ysize, xsize))
    else:
        sarped = None

    if slice_del is not None:
        darped = tf.warp(slice_del, tform, order=3,
                         output_shape=(ysize, xsize), preserve_range=True)
    else:
        darped = None

    return argument['slice_number'], warped, varped, marped, farped, \
        oarped, sarped, darped


class MakeCube(BasePrimitive):
    """Transform 2D images to 3D data cubes"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Creating data cube")

        log_string = MakeCube.__module__

        # Are we interactive?
        if self.config.instrument.plot_level >= 3:
            do_inter = True
        else:
            do_inter = False
        self.logger.info("Generating data cube")
        # Find and read geometry transformation
        tab = self.context.proctab.n_proctab(frame=self.action.args.ccddata,
                                             target_type='ARCLAMP',
                                             nearest=True)
        if not len(tab):
            self.logger.error("No reference geometry, cannot make cube!")
            self.action.args.ccddata.header['GEOMCOR'] = (False,
                                                          'Geometry corrected?')
            self.logger.info(log_string)
            return self.action.args

        self.logger.info("%d arc frames found" % len(tab))
        ofname = tab['OFNAME'][0]
        geom_file = os.path.join(self.config.instrument.output_directory,
                                 ofname.split('.')[0] + '_geom.pkl')
        if os.path.exists(geom_file):
            self.logger.info("Reading %s" % geom_file)
            with open(geom_file, 'rb') as ifile:
                geom = pickle.load(ifile)
            # Slice size
            xsize = geom['xsize']
            ysize = geom['ysize']
            out_cube = np.zeros((ysize, xsize, 24), dtype=np.float64)
            out_vube = np.zeros((ysize, xsize, 24), dtype=np.float64)
            out_mube = np.zeros((ysize, xsize, 24), dtype=np.uint8)
            out_fube = np.zeros((ysize, xsize, 24), dtype=np.uint8)
            out_oube = np.zeros((ysize, xsize, 24), dtype=np.float64)
            out_sube = np.zeros((ysize, xsize, 24), dtype=np.float64)
            out_dube = np.zeros((ysize, xsize, 24), dtype=np.float64)
            # Store original data
            data_img = self.action.args.ccddata.data
            data_std = self.action.args.ccddata.uncertainty.array
            data_msk = self.action.args.ccddata.mask
            data_flg = self.action.args.ccddata.flags
            # check for obj, sky images
            obj = None
            data_obj = None
            sky = None
            data_sky = None
            if self.action.args.nasmask and self.action.args.numopen > 1:
                ofn = self.action.args.ccddata.header['OFNAME']

                objfn = ofn.split('.')[0] + '_objf.fits'
                full_path = os.path.join(
                    os.path.dirname(self.action.args.name),
                    self.config.instrument.output_directory, objfn)
                if os.path.exists(full_path):
                    obj = kcwi_fits_reader(full_path)[0]
                    data_obj = obj.data

                skyfn = ofn.split('.')[0] + '_skyf.fits'
                full_path = os.path.join(
                    os.path.dirname(self.action.args.name),
                    self.config.instrument.output_directory, skyfn)
                if os.path.exists(full_path):
                    sky = kcwi_fits_reader(full_path)[0]
                    data_sky = sky.data
            # check for geometry maps
            dew = None
            data_dew = None
            if 'ARCLAMP' in self.action.args.imtype:
                ofn = self.action.args.ccddata.header['OFNAME']
                dewfn = ofn.split('.')[0] + '_delmap.fits'
                full_path = os.path.join(
                    os.path.dirname(self.action.args.name),
                    self.config.instrument.output_directory, dewfn)
                if os.path.exists(full_path):
                    dew = kcwi_fits_reader(full_path)[0]
                    data_dew = dew.data
            # Loop over 24 slices
            my_arguments = []
            for isl in range(0, 24):
                arguments = {
                    'slice_number': isl,
                    'geom': geom,
                    'img': data_img,
                    'std': data_std,
                    'msk': data_msk,
                    'flg': data_flg,
                    'xsize': xsize,
                    'ysize': ysize,
                    'logger': self.logger
                }
                if obj is not None:
                    arguments['obj'] = data_obj
                if sky is not None:
                    arguments['sky'] = data_sky
                if dew is not None:
                    arguments['del'] = data_dew
                my_arguments.append(arguments)

            p = Pool()
            results = p.map(make_cube_helper, list(my_arguments))
            p.close()

            self.logger.info("Building cube")
            for partial_cube in results:
                slice_number = partial_cube[0]
                out_cube[:, :, slice_number] = partial_cube[1]
                out_vube[:, :, slice_number] = partial_cube[2]
                out_mube[:, :, slice_number] = partial_cube[3]
                out_fube[:, :, slice_number] = partial_cube[4]
                if obj is not None:
                    out_oube[:, :, slice_number] = partial_cube[5]
                if sky is not None:
                    out_sube[:, :, slice_number] = partial_cube[6]
                if dew is not None:
                    out_dube[:, :, slice_number] = partial_cube[7]

            if self.config.instrument.plot_level >= 3:
                for isl in range(0, 24):
                    warped = out_cube[:, :, isl]
                    ptitle = self.action.args.plotlabel + \
                        "WARPED Slice %d" % isl
                    p = figure(tooltips=[("x", "$x"), ("y", "$y"),
                                         ("value", "@image")],
                               title=ptitle,
                               x_axis_label="X (px)", y_axis_label="Y (px)",
                               plot_width=self.config.instrument.plot_width,
                               plot_height=self.config.instrument.plot_height)
                    p.x_range.range_padding = p.y_range.range_padding = 0
                    p.image([warped], x=0, y=0, dw=xsize, dh=ysize,
                            palette="Spectral11", level="image")
                    bokeh_plot(p, self.context.bokeh_session)
                    if do_inter:
                        q = input("Next? <cr>, q to quit: ")
                        if 'Q' in q.upper():
                            do_inter = False
                    else:
                        time.sleep(self.config.instrument.plot_pause)

            # Calculate some WCS parameters
            # Get object pointing
            try:
                if self.action.args.nasmask:
                    rastr = self.action.args.ccddata.header['RABASE']
                    decstr = self.action.args.ccddata.header['DECBASE']
                else:
                    rastr = self.action.args.ccddata.header['RA']
                    decstr = self.action.args.ccddata.header['DEC']
            except KeyError:
                try:
                    rastr = self.action.args.ccddata.header['TARGRA']
                    decstr = self.action.args.ccddata.header['TARGDEC']
                except KeyError:
                    rastr = ''
                    decstr = ''
            if len(rastr) > 0 and len(decstr) > 0:
                try:
                    coord = SkyCoord(rastr, decstr, unit=(u.hourangle, u.deg))
                except ValueError:
                    coord = None
            else:
                coord = None
            # Get rotator position
            if 'ROTPOSN' in self.action.args.ccddata.header:
                rpos = self.action.args.ccddata.header['ROTPOSN']
            else:
                rpos = 0.
            if 'ROTREFAN' in self.action.args.ccddata.header:
                rref = self.action.args.ccddata.header['ROTREFAN']
            else:
                rref = 0.
            skypa = rpos + rref
            crota = math.radians(-(skypa + self.config.instrument.ROTOFF))
            cdelt1 = -geom['slscl']
            cdelt2 = geom['pxscl']
            if coord is None:
                ra = 0.
                dec = 0.
                crota = 1
            else:
                ra = coord.ra.degree
                dec = coord.dec.degree
            cd11 = cdelt1 * math.cos(crota)
            cd12 = abs(cdelt2) * np.sign(cdelt1) * math.sin(crota)
            cd21 = -abs(cdelt1) * np.sign(cdelt2) * math.sin(crota)
            cd22 = cdelt2 * math.cos(crota)
            crpix1 = 12.
            crpix2 = xsize / 2.
            crpix3 = 1.
            porg = self.action.args.ccddata.header['PONAME']
            ifunum = self.action.args.ifunum
            if 'IFU' in porg:
                if ifunum == 1:
                    off1 = 1.0
                    off2 = 4.0
                elif ifunum == 2:
                    off1 = 1.0
                    off2 = 5.0
                elif ifunum == 3:
                    off1 = 0.05
                    off2 = 5.6
                else:
                    self.logger.warning("Unknown IFU number: %d" % ifunum)
                    off1 = 0.
                    off2 = 0.
                off1 = off1 / float(self.action.args.xbinsize)
                off2 = off2 / float(self.action.args.ybinsize)
                crpix1 += off1
                crpix2 += off2
            # Update header
            # Geometry corrected?
            self.action.args.ccddata.header['GEOMCOR'] = (
                True, 'Geometry corrected?')
            #
            # Spatial geometry
            self.action.args.ccddata.header['BARSEP'] = (
                geom['barsep'], 'separation of bars (binned pix)')
            self.action.args.ccddata.header['BAR0'] = (
                geom['bar0'], 'first bar pixel position')
            # Wavelength ranges
            if self.action.args.nasmask:
                self.action.args.ccddata.header['WAVALL0'] = (
                    geom['wavensall0'], 'Low inclusive wavelength')
                self.action.args.ccddata.header['WAVALL1'] = (
                    geom['wavensall1'], 'High inclusive wavelength')
                self.action.args.ccddata.header['WAVGOOD0'] = (
                    geom['wavensgood0'], 'Low good wavelength')
                self.action.args.ccddata.header['WAVGOOD1'] = (
                    geom['wavensgood1'], 'High good wavelength')
                self.action.args.ccddata.header['WAVMID'] = (
                    geom['wavensmid'], 'middle wavelength')
            else:
                self.action.args.ccddata.header['WAVALL0'] = (
                    geom['waveall0'], 'Low inclusive wavelength')
                self.action.args.ccddata.header['WAVALL1'] = (
                    geom['waveall1'], 'High inclusive wavelength')
                self.action.args.ccddata.header['WAVGOOD0'] = (
                    geom['wavegood0'], 'Low good wavelength')
                self.action.args.ccddata.header['WAVGOOD1'] = (
                    geom['wavegood1'], 'High good wavelength')
                self.action.args.ccddata.header['WAVMID'] = (
                    geom['wavemid'], 'middle wavelength')
            # Dichroic fraction
            try:
                dichroic_fraction = geom['dich_frac']
            except AttributeError:
                dichroic_fraction = 1.
            self.action.args.ccddata.header['DICHFRAC'] = (
                dichroic_fraction, 'Dichroic Fraction')
            # Wavelength fit statistics
            self.action.args.ccddata.header['AVWVSIG'] = (
                geom['avwvsig'], 'Avg. bar wave sigma (Ang)')
            self.action.args.ccddata.header['SDWVSIG'] = (
                geom['sdwvsig'], 'Stdev. var wave sigma (Ang)')
            # Pixel scales
            self.action.args.ccddata.header['PXSCL'] = (
                geom['pxscl'], 'Pixel scale along slice (deg)')
            self.action.args.ccddata.header['SLSCL'] = (
                geom['slscl'], 'Pixel scale perp. to slices (deg)')
            # Geometry origins
            self.action.args.ccddata.header['CBARSNO'] = (
                geom['cbarsno'], 'Continuum bars image number')
            self.action.args.ccddata.header['CBARSFL'] = (
                geom['cbarsfl'], 'Continuum bars image filename')
            self.action.args.ccddata.header['ARCNO'] = (
                geom['arcno'], 'Arc image number')
            self.action.args.ccddata.header['ARCFL'] = (
                geom['arcfl'], 'Arc image filename')
            self.action.args.ccddata.header['GEOMFL'] = (
                geom_file.split('/')[-1], 'Geometry file')
            # WCS
            self.action.args.ccddata.header['IFUPA'] = (
                skypa, 'IFU position angle (degrees)')
            self.action.args.ccddata.header['IFUROFF'] = (
                self.config.instrument.ROTOFF, 'IFU-SKYPA offset (degrees)')
            self.action.args.ccddata.header['WCSDIM'] = (
                3, 'number of dimensions in WCS')
            self.action.args.ccddata.header['WCSNAME'] = 'KCWI'
            self.action.args.ccddata.header['EQUINOX'] = 2000.
            self.action.args.ccddata.header['RADESYS'] = 'FK5'
            self.action.args.ccddata.header['CTYPE1'] = 'RA---TAN'
            self.action.args.ccddata.header['CTYPE2'] = 'DEC--TAN'
            self.action.args.ccddata.header['CTYPE3'] = ('AWAV',
                                                         'Air Wavelengths')
            self.action.args.ccddata.header['CUNIT1'] = ('deg', 'RA units')
            self.action.args.ccddata.header['CUNIT2'] = ('deg', 'DEC units')
            self.action.args.ccddata.header['CUNIT3'] = ('Angstrom',
                                                         'Wavelength units')
            self.action.args.ccddata.header['CNAME1'] = ('KCWI RA', 'RA name')
            self.action.args.ccddata.header['CNAME2'] = ('KCWI DEC', 'DEC name')
            self.action.args.ccddata.header['CNAME3'] = ('KCWI Wavelength',
                                                         'Wavelength name')
            self.action.args.ccddata.header['CRVAL1'] = (ra, 'RA zeropoint')
            self.action.args.ccddata.header['CRVAL2'] = (dec, 'DEC zeropoint')
            self.action.args.ccddata.header['CRVAL3'] = (geom['wave0out'],
                                                         'Wavelength zeropoint')
            self.action.args.ccddata.header['CRPIX1'] = (crpix1,
                                                         'RA reference pixel')
            self.action.args.ccddata.header['CRPIX2'] = (crpix2,
                                                         'DEC reference pixel')
            self.action.args.ccddata.header['CRPIX3'] = (
                crpix3, 'Wavelength reference pixel')
            self.action.args.ccddata.header['CD1_1'] = (
                cd11, 'RA degrees per column pixel')
            self.action.args.ccddata.header['CD2_1'] = (
                cd21, 'DEC degrees per column pixel')
            self.action.args.ccddata.header['CD1_2'] = (
                cd12, 'RA degrees per row pixel')
            self.action.args.ccddata.header['CD2_2'] = (
                cd22, 'DEC degrees per row pixel')
            self.action.args.ccddata.header['CD3_3'] = (
                geom['dwout'], 'Wavelength Angstroms per pixel')
            self.action.args.ccddata.header['LONPOLE'] = (
                180.0, 'Native longitude of Celestial pole')
            self.action.args.ccddata.header['LATPOLE'] = (
                0.0, 'Native latitude of Celestial pole')
            # write out cube
            self.action.args.ccddata.header['HISTORY'] = log_string
            self.action.args.ccddata.data = out_cube
            self.action.args.ccddata.uncertainty.array = out_vube
            self.action.args.ccddata.mask = out_mube
            self.action.args.ccddata.flags = out_fube

            # write out int image
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name,
                             output_dir=self.config.instrument.output_directory,
                             suffix="icube")
            self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                                suffix="icube")
            self.context.proctab.write_proctab()

            # check for obj, sky outputs
            if obj is not None:
                out_obj = CCDData(out_oube,
                                  meta=self.action.args.ccddata.header,
                                  unit=self.action.args.ccddata.unit)
                kcwi_fits_writer(
                    out_obj, output_file=self.action.args.name,
                    output_dir=self.config.instrument.output_directory,
                    suffix="ocube")
            if sky is not None:
                out_sky = CCDData(out_sube,
                                  meta=self.action.args.ccddata.header,
                                  unit=self.action.args.ccddata.unit)
                kcwi_fits_writer(
                    out_sky, output_file=self.action.args.name,
                    output_dir=self.config.instrument.output_directory,
                    suffix="scube")
            # check for dew outputs
            if dew is not None:
                out_dew = CCDData(out_dube,
                                  meta=self.action.args.ccddata.header,
                                  unit=dew.unit)
                kcwi_fits_writer(
                    out_dew, output_file=self.action.args.name,
                    output_dir=self.config.instrument.output_directory,
                    suffix="dcube")
        else:
            self.logger.error("Geometry file not found: %s" % geom_file)

        self.logger.info(log_string)

        return self.action.args
    # END: class MakeCube()
