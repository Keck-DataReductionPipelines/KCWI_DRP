from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer

import os
import numpy as np
import pickle
from astropy.nddata import CCDData
from astropy import units as u


class GenerateMaps(BasePrimitive):
    """Generate map images"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Generating geometry maps")

        log_string = GenerateMaps.__module__

        if self.action.args.geometry_file is not None and \
                os.path.exists(self.action.args.geometry_file):
            with open(self.action.args.geometry_file, 'rb') as ifile:
                geom = pickle.load(ifile)
            # get geom params
            xl0s = geom['xl0']  # lower slice pos limit
            xl1s = geom['xl1']  # upper slice pos limit
            invtf_list = geom['invtf']
            wave0 = geom['wave0out']
            dw = geom['dwout']
            xsize = geom['xsize']
            # Store original data
            data_img = self.action.args.ccddata.data
            ny = data_img.shape[0]  # number of wavelength pixels
            # Create map images
            wave_map_img = np.full_like(data_img, fill_value=-1.)
            xpos_map_img = np.full_like(data_img, fill_value=-1.)
            slice_map_img = np.full_like(data_img, fill_value=-1.)
            delta_map_img = np.full_like(data_img, fill_value=-1.)
            # loop over slices
            for isl in range(0, 24):
                itrf = invtf_list[isl]
                xl0 = xl0s[isl]
                xl1 = xl1s[isl]
                # loop over slice position
                for ix in range(xl0, xl1):
                    coords = np.zeros((ny, 2))
                    # loop over wavelength pixels
                    for iy in range(0, ny):
                        coords[iy, 0] = ix - xl0
                        coords[iy, 1] = iy
                    ncoo = itrf(coords)
                    # loop over wavelength
                    for iy in range(0, ny):
                        if 0 <= ncoo[iy, 0] <= xsize:
                            slice_map_img[iy, ix] = isl
                            xpos_map_img[iy, ix] = ncoo[iy, 0]
                            wave_map_img[iy, ix] = ncoo[iy, 1] * dw + wave0
                            if iy > 0:
                                delta_map_img[iy, ix] = abs(
                                    (ncoo[iy, 1] - ncoo[iy-1, 1]) * dw)

            # update header
            self.action.args.ccddata.header['HISTORY'] = log_string
            # Spatial geometry
            self.action.args.ccddata.header['BARSEP'] = (
                geom['barsep'], 'separation of bars (binned pix)')
            self.action.args.ccddata.header['BAR0'] = (
                geom['bar0'], 'first bar pixel position')
            # Dichroic fraction
            try:
                dichroic_fraction = geom['dich_frac']
            except AttributeError:
                dichroic_fraction = 1.
            self.action.args.ccddata.header['DICHFRAC'] = (
                dichroic_fraction, 'Dichroic Fraction')
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
                self.action.args.geometry_file.split('/')[-1], 'Geometry file')

            # output maps
            header = self.action.args.ccddata.header

            kcwi_fits_writer(CCDData(wave_map_img, meta=header,
                                     unit=u.angstrom),
                             output_file=self.action.args.name,
                             output_dir=self.config.instrument.output_directory,
                             suffix="wavemap")
            kcwi_fits_writer(CCDData(xpos_map_img, meta=header,
                                     unit=u.pix),
                             output_file=self.action.args.name,
                             output_dir=self.config.instrument.output_directory,
                             suffix="posmap")
            kcwi_fits_writer(CCDData(slice_map_img, meta=header,
                                     unit=u.pix),
                             output_file=self.action.args.name,
                             output_dir=self.config.instrument.output_directory,
                             suffix="slicemap")
            kcwi_fits_writer(CCDData(delta_map_img, meta=header,
                                     unit=u.angstrom),
                             output_file=self.action.args.name,
                             output_dir=self.config.instrument.output_directory,
                             suffix="delmap")

        else:
            self.logger.error("Geom file not accessible")

        self.logger.info(log_string)

        return self.action.args
    # END: class GenerateMaps()
