from astropy.io import fits as pf
import pkg_resources
import os

import numpy as np
from scipy.interpolate import interp1d


def kcwi_correct_extin(img, hdr, logger=None):
    """Atmospheric extinction correction"""
    # get airmass
    air = hdr['AIRMASS']
    # read extinction data
    path = 'data/extin/snfext.fits'
    package = __name__.split('.')[0]
    full_path = pkg_resources.resource_filename(package, path)
    if os.path.exists(full_path):
        hdul = pf.open(full_path)
        exwl = hdul[1].data['LAMBDA']
        exma = hdul[1].data['EXT']
        # get object wavelengths
        sz = img.shape
        dw = hdr['CD3_3']
        w0 = hdr['CRVAL3']
        owls = np.arange(sz[0]) * dw + w0
        # linear interpolation
        exint = interp1d(exwl, exma, kind='cubic', bounds_error=False,
                         fill_value='extrapolate')
        # resample extinction curve
        oexma = exint(owls)
        # convert to flux ratio
        flxr = 10.**(oexma * air * 0.4)
        if len(sz) == 3:
            # apply to cube
            for ix in range(sz[2]):
                for iy in range(sz[1]):
                    img[:, iy, ix] *= flxr
        else:
            # apply to vector
            img *= flxr

        flrmn = np.nanmean(flxr)
        hdr['HISTORY'] = 'kcwi_correct_extin'
        hdr['EXTCOR'] = (True, 'extinction corrected?')
        hdr['AVEXCOR'] = (flrmn, 'average extin. correction (flux ratio)')
        if logger:
            logger.info("Extinction corrected")
        else:
            print("Extinction corrected")
    else:
        if logger:
            logger.warning("Extinction data file (%s) not found!" % full_path)
        else:
            print("Extinction data file (%s) not found!" % full_path)
