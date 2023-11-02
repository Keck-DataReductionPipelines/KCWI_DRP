from keckdrpframework.primitives.base_primitive import BasePrimitive
from kcwidrp.primitives.kcwi_file_primitives import kcwi_fits_writer

import numpy as np
import pkg_resources
import os
import pandas as pd


class CorrectDefects(BasePrimitive):
    """
    Remove known bad columns.

    Looks for a defect list file in the data directory of kcwidrp based on the
    CCD ampmode and x and y binning.  Records the defect correction in the
    FITS header with the following keywords:

        * BPFILE: the bad pixel file used to correct defects
        * NBPCLEAN: the number of bad pixels cleaned

    Uses the following configuration parameter:

        * saveintims: if set to ``True`` write out a \*_def.fits file with defects corrected.  Default is ``False``.

    Updates image in returned arguments.

    """

    def __init__(self, action, context):
        BasePrimitive.__init__(self, action, context)
        self.logger = context.pipeline_logger

    def _perform(self):
        self.logger.info("Correcting detector defects")

        # Header keyword to update
        key = 'BPCLEAN'
        keycom = 'cleaned bad pixels?'

        # Create flags for bad columns fixed
        if self.action.args.ccddata.flags is None:
            self.action.args.ccddata.flags = np.zeros(
                self.action.args.ccddata.data.shape, dtype=np.uint8)

        flags = self.action.args.ccddata.flags

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
            pixel_range_for_good_value = 2
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
                                  x1:x1+pixel_range_for_good_value])
                    # get replacement value
                    good_values = np.nanmedian(np.asarray(values))
                    # Replace baddies with good_values
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
        # DN 2023-may-28: commenting out because it causes bad things later on
        # self.action.args.ccddata.mask = flags
        self.action.args.ccddata.flags = flags

        if self.config.instrument.saveintims:
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name,
                             output_dir=self.config.instrument.output_directory,
                             suffix="def")

        return self.action.args
    # END: class CorrectDefects()
