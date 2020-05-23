from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer

import numpy as np
from astropy.nddata import CCDData


class CubeImage(BasePrimitive):
    """Transform 2D images to 3D data cubes"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _pre_condition(self):
        """Checks if we have a cube to convert"""
        self.logger.info("Checking precondition for CubeImage")
        sz = self.action.args.ccddata.data.shape
        if len(sz) == 3:
            if sz[2] == 24:
                self.logger.info("CubeImage precondition met.")
                self.action.args.cube_size = sz
                return True
            else:
                self.logger.warning("Mal-formed cube!")
                return False
        else:
            self.logger.warning("Not a 3D data cube!")
            return False

    def _perform(self):
        self.logger.info("Creating 2D image from 3D cube")

        log_string = CubeImage.__module__

        # get wavelength ranges
        wall0 = self.action.args.ccddata.header['WAVALL0']
        wall1 = self.action.args.ccddata.header['WAVALL1']
        w0 = self.action.args.ccddata.header['CRVAL3']
        dw = self.action.args.ccddata.header['CD3_3']
        crpixw = self.action.args.ccddata.header['CRPIX3']
        y0 = int((wall0 - w0) / dw)
        if y0 < 0:
            y0 = 0
        out_w0 = y0 * dw + w0
        y1 = int((wall1 - w0) / dw)
        if y1 > self.action.args.cube_size[0]:
            y1 = self.action.args.cube_size[0]
        out_y = y1 - y0
        # create output image
        cub_x = self.action.args.cube_size[1]
        out_x = int(24 * cub_x)
        out_img = np.zeros((out_y, out_x), dtype=np.float)
        # set spatial scale
        s0 = 0.
        ds = 24.0 / out_x
        crpixs = 1.
        self.logger.info("cube dims: y, x, z: %d, %d, %d" %
                         self.action.args.cube_size)
        self.logger.info("y0, y1 = %d, %d: out_x, out_y = %d, %d" %
                         (y0, y1, out_x, out_y))

        # loop over slices
        for isl in range(24):
            x0 = isl * cub_x
            x1 = (isl + 1) * cub_x
            out_img[:, x0:x1] = self.action.args.ccddata.data[y0:y1, :, isl]
        # output CCDData structure
        out_ccd = CCDData(out_img, meta=self.action.args.ccddata.header,
                          unit=self.action.args.ccddata.unit)

        # update header
        out_ccd.header['HISTORY'] = log_string
        del out_ccd.header['RADESYS']
        del out_ccd.header['LONPOLE']
        del out_ccd.header['LATPOLE']
        del out_ccd.header['CTYPE1']
        del out_ccd.header['CTYPE2']
        del out_ccd.header['CTYPE3']
        del out_ccd.header['CUNIT1']
        del out_ccd.header['CUNIT2']
        del out_ccd.header['CUNIT3']
        del out_ccd.header['CNAME1']
        del out_ccd.header['CNAME2']
        del out_ccd.header['CNAME3']
        del out_ccd.header['CRVAL1']
        del out_ccd.header['CRVAL2']
        del out_ccd.header['CRVAL3']
        del out_ccd.header['CRPIX1']
        del out_ccd.header['CRPIX2']
        del out_ccd.header['CRPIX3']
        del out_ccd.header['CD1_1']
        del out_ccd.header['CD1_2']
        del out_ccd.header['CD2_1']
        del out_ccd.header['CD2_2']
        del out_ccd.header['CD3_3']

        out_ccd.header['WCSDIM'] = 2

        out_ccd.header['CTYPE1'] = ('SPATIAL', 'slice number')
        out_ccd.header['CUNIT1'] = ('slu', 'slice units')
        out_ccd.header['CNAME1'] = ('KCWI Slice', 'slice name')
        out_ccd.header['CRVAL1'] = (s0, 'slice zeropoint')
        out_ccd.header['CRPIX1'] = (crpixs, 'slice reference pixel')
        out_ccd.header['CDELT1'] = (ds, 'slice units per pixel')

        out_ccd.header['CTYPE2'] = ('AWAV', 'Air Wavelengths')
        out_ccd.header['CUNIT2'] = ('Angstrom', 'Wavelength units')
        out_ccd.header['CNAME2'] = ('KCWI Wavelength', 'Wavelength name')
        out_ccd.header['CRVAL2'] = (out_w0, 'Wavelength zeropoint')
        out_ccd.header['CRPIX2'] = (crpixw, 'Wavelength reference pixel')
        out_ccd.header['CDELT2'] = (dw, 'Wavelength Angstroms per pixel')

        # write out image
        kcwi_fits_writer(out_ccd, output_file=self.action.args.name,
                         output_dir=self.config.instrument.output_directory,
                         suffix="icube_2d")

        self.logger.info(log_string)

        return self.action.args
    # END: class CubeImage()
