from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer

import os
import numpy as np
import pickle


class GenerateMaps(BasePrimitive):
    """Generate map images"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Generating geometry maps")

        log_string = GenerateMaps.__module__ + "." + GenerateMaps.__qualname__

        if self.action.args.geom_file is not None and \
                os.path.exists(self.action.args.geom_file):
            with open(self.action.args.geom_file, 'rb') as ifile:
                geom = pickle.load(ifile)
            # get geometry params
            xl0s = geom['xl0']
            xl1s = geom['xl1']
            invtf_list = geom['invtf']
            wave0 = geom['wave0out']
            dw = geom['dwout']
            xsize = geom['xsize']
            # Store original data
            data_img = self.action.args.ccddata.data
            ny = data_img.shape[0]
            # Create map images
            wave_map_img = np.full_like(data_img, fill_value=-1.)
            xpos_map_img = np.full_like(data_img, fill_value=-1.)
            slice_map_img = np.full_like(data_img, fill_value=-1.)
            # loop over slices
            for isl in range(0, 24):
                itrf = invtf_list[isl]
                xl0 = xl0s[isl]
                xl1 = xl1s[isl]
                for ix in range(xl0, xl1):
                    coords = np.zeros((ny, 2))
                    for iy in range(0, ny):
                        coords[iy, 0] = ix - xl0
                        coords[iy, 1] = iy
                    ncoo = itrf(coords)
                    for iy in range(0, ny):
                        if 0 <= ncoo[iy, 0] <= xsize:
                            slice_map_img[iy, ix] = isl
                            xpos_map_img[iy, ix] = ncoo[iy, 0]
                            wave_map_img[iy, ix] = ncoo[iy, 1] * dw + wave0

            self.action.args.ccddata.header['HISTORY'] = log_string

            # output maps
            self.action.args.ccddata.data = wave_map_img
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name,
                             suffix="wavemap")
            self.action.args.ccddata.data = xpos_map_img
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name,
                             suffix="posmap")
            self.action.args.ccddata.data = slice_map_img
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name,
                             suffix="slicemap")
            self.action.args.ccddata.data = data_img

        else:
            self.logger.error("Geom file not accessible")

        self.logger.info(log_string)

        return self.action.args
    # END: class GenerateMaps()


