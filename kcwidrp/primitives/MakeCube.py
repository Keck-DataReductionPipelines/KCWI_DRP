from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer

import os
import math
import pickle
import numpy as np
from skimage import transform as tf
from astropy.coordinates import SkyCoord
from astropy import units as u


class MakeCube(BasePrimitive):
    """Transform 2D images to 3D data cubes"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Creating data cube")

        log_string = MakeCube.__module__ + "." + MakeCube.__qualname__

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
            # Store original data
            data_img = self.action.args.ccddata.data
            data_var = self.action.args.ccddata.uncertainty.array
            data_msk = self.action.args.ccddata.mask
            # Loop over 24 slices
            for isl in range(0, 24):
                tform = geom['tform'][isl]
                xl0 = geom['xl0'][isl]
                xl1 = geom['xl1'][isl]
                self.logger.info("Transforming image slice %d" % isl)
                slice_img = data_img[:, xl0:xl1]
                slice_var = data_var[:, xl0:xl1]
                slice_msk = data_msk[:, xl0:xl1]
                # wmed = np.nanmedian(slice_img)
                # wstd = np.nanstd(slice_img)
                # pl.clf()
                # pl.imshow(slice_img, vmin=(wmed - wstd * 2.),
                #          vmax=(wmed + wstd * 2.))
                # pl.ylim(0, 2056)
                # pl.title('raw slice %d' % isl)
                # if do_inter:
                #    q = input("Next? <cr>, q to quit: ")
                #    if 'Q' in q.upper():
                #        do_inter = False
                #        pl.ioff()
                # else:
                #    pl.pause(self.action.args.ccddata.plotpause())
                warped = tf.warp(slice_img, tform, order=3,
                                 output_shape=(ysize, xsize))
                varped = tf.warp(slice_var, tform, order=3,
                                 output_shape=(ysize, xsize))
                marped = tf.warp(slice_msk, tform, order=3,
                                 output_shape=(ysize, xsize))
                for iy in range(ysize):
                    for ix in range(xsize):
                        out_cube[iy, ix, isl] = warped[iy, ix]
                        out_vube[iy, ix, isl] = varped[iy, ix]
                        out_mube[iy, ix, isl] = int(marped[iy, ix])
                wmed = np.nanmedian(warped)
                wstd = np.nanstd(warped)
                if self.config.instrument.plot_level >= 2:
                    print(wmed, wstd)
                    ptitle = self.action.args.plotlabel + "WARPED Slice %d" \
                        % isl
                    p = figure(tooltips=[("x", "$x"), ("y", "$y"),
                                         ("value", "@image")],
                               title=ptitle,
                               x_axis_label="X (px)", y_axis_label="Y (px)",
                               plot_width=self.config.instrument.plot_width,
                               plot_height=self.config.instrument.plot_height)
                    p.x_range.range_padding = p.y_range.range_padding = 0
                    p.image([warped], x=0, y=0, dw=xsize, dh=ysize,
                            palette="Spectral11", level="image")
                    # vmin=(wmed - wstd * 2.), vmax=(wmed + wstd * 2.))
                    # pl.ylim(0, ysize)
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
                coord = SkyCoord(rastr, decstr, unit=(u.hourangle, u.deg))
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
            # Wavelength fit statistics
            self.action.args.ccddata.header['AVWVSIG'] = (
                geom['avwvsig'], 'Avg. bar wave sigma (Ang)')
            self.action.args.ccddata.header['SDWVSIG'] = (
                geom['sdwvsig'], 'Stdev. var wave sigma (Ang)')
            # Pixel scales
            self.action.args.ccddata.header['PXSCL'] = (
                geom['pxscl'], 'Pixel scale along slice')
            self.action.args.ccddata.header['SLSCL'] = (
                geom['slscl'], 'Pixel scale perpendicular to slices')
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

            # write out int image
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name, suffix="icube")
            self.context.proctab.update_proctab(frame=self.action.args.ccddata,
                                                suffix="icube")
            self.context.proctab.write_proctab()
        else:
            self.logger.error("Geometry file not found: %s" % geom_file)

        self.logger.info(log_string)

        return self.action.args
    # END: class MakeCube()


