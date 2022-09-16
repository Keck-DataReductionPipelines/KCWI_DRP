from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer

import numpy as np
import pkg_resources
import os
import pandas as pd


class CorrectDefects(BasePrimitive):
    """Remove known bad columns"""

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Correcting detector defects")

        # Header keyword to update
        key = 'BPCLEAN'
        keycom = 'cleaned bad pixels?'

        # Create flags for bad columns fixed
        flags = np.zeros(self.action.args.ccddata.data.shape, dtype=np.uint8)

        # Nod and Shuffle?
        if self.action.args.nasmask and self.action.args.numopen > 1:
            nastr = "_nas"
        else:
            nastr = ""

        # Does the defect file exist?
        path = "data/defect_%s_%dx%d%s.dat" % (self.action.args.ampmode.strip(),
                                               self.action.args.xbinsize,
                                               self.action.args.ybinsize, nastr)
        package = __name__.split('.')[0]
        full_path = pkg_resources.resource_filename(package, path)
        number_of_bad_pixels = 0   # count of defective pixels cleaned
        if os.path.exists(full_path):
            self.logger.info("Reading defect list in: %s" % full_path)
            defect_table = pd.read_csv(full_path, sep=r'\s+')
            # range of pixels for calculating good value
            pixel_range_for_good_value = 5
            for index, row in defect_table.iterrows():
                # Get coords and adjust for python zero bias
                x0 = row['X0'] - 1
                x1 = row['X1']
                y0 = row['Y0'] - 1
                y1 = row['Y1']
                # Loop over y range
                for by in range(y0, y1):
                    # sample on low side of bad area
                    values = list(self.action.args.ccddata.data[by,
                                  x0-pixel_range_for_good_value:x0])
                    # sample on high side
                    values.extend(self.action.args.ccddata.data[by,
                                  x1+1:x1+pixel_range_for_good_value+1])
                    # get replacement value
                    good_values = np.nanmedian(np.asarray(values))
                    # Replace baddies with gval
                    for bx in range(x0, x1):
                        self.action.args.ccddata.data[by, bx] = good_values
                        flags[by, bx] += 2
                        number_of_bad_pixels += 1
            self.action.args.ccddata.header[key] = (True, keycom)
            self.action.args.ccddata.header['BPFILE'] = (path, 'defect list')
        else:
            self.logger.error("Defect list not found for %s" % full_path)
            self.action.args.ccddata.header[key] = (False, keycom)

        self.logger.info("Cleaned %d bad pixels" % number_of_bad_pixels)
        self.action.args.ccddata.header['NBPCLEAN'] = \
            (number_of_bad_pixels, 'number of bad pixels cleaned')

        log_string = CorrectDefects.__module__
        self.action.args.ccddata.header['HISTORY'] = log_string
        self.logger.info(log_string)

        # add flags array
        self.action.args.ccddata.mask = flags
        self.action.args.ccddata.flags = flags

        if self.config.instrument.saveintims:
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name,
                             output_dir=self.config.instrument.output_directory,
                             suffix="def")

        return self.action.args
    # END: class CorrectDefects()
