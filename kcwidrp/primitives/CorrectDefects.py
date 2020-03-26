from keckdrpframework.primitives.base_primitive import BasePrimitive
from keckdrpframework.models.arguments import Arguments

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

        # Does the defect file exist?
        path = "data/defect_%s_%dx%d.dat" % (self.action.args.ampmode.strip(),
                                             self.action.args.xbinsize,
                                             self.action.args.ybinsize)
        pkg = __name__.split('.')[0]
        defpath = pkg_resources.resource_filename(pkg, path)
        nbpix = 0   # count of defective pixels cleaned
        if os.path.exists(defpath):
            self.logger.info("Reading defect list in: %s" % defpath)
            deftab = pd.read_csv(defpath, sep=r'\s+')
            bcdel = 5   # range of pixels for calculating good value
            for indx, row in deftab.iterrows():
                # Get coords and adjust for python zero bias
                x0 = row['X0'] - 1
                x1 = row['X1']
                y0 = row['Y0'] - 1
                y1 = row['Y1']
                # Loop over y range
                for by in range(y0, y1):
                    # sample on low side of bad area
                    vals = list(self.action.args.ccddata.data[by,
                                x0-bcdel:x0])
                    # sample on high side
                    vals.extend(self.action.args.ccddata.data[by,
                                x1+1:x1+bcdel+1])
                    # get replacement value
                    gval = np.nanmedian(np.asarray(vals))
                    # Replace baddies with gval
                    for bx in range(x0, x1):
                        self.action.args.ccddata.data[by, bx] = gval
                        flags[by, bx] += 2
                        nbpix += 1
            self.action.args.ccddata.header[key] = (True, keycom)
            self.action.args.ccddata.header['BPFILE'] = (path, 'defect list')
        else:
            self.logger.error("Defect list not found for %s" % defpath)
            self.action.args.ccddata.header[key] = (False, keycom)

        self.logger.info("Cleaned %d bad pixels" % nbpix)
        self.action.args.ccddata.header['NBPCLEAN'] = \
            (nbpix, 'number of bad pixels cleaned')

        logstr = CorrectDefects.__module__ + "." + CorrectDefects.__qualname__
        self.action.args.ccddata.header['HISTORY'] = logstr
        self.logger.info(logstr)

        # add flags array
        self.action.args.ccddata.mask = flags
        self.action.args.ccddata.flags = flags

        if self.config.instrument.saveintims:
            kcwi_fits_writer(self.action.args.ccddata,
                             table=self.action.args.table,
                             output_file=self.action.args.name, suffix="def")

        return self.action.args
    # END: class CorrectDefects()

