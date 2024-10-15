import os
from astropy.io.fits import file
import pkg_resources
from astropy.io import fits


def kcwi_get_std(targname, logger=None):
    """Checks if object is a standard star, returns std file path and name"""

    stdfile = None
    stdname = None
    obname = targname.lower()
    obname = "".join(obname.split()) # remove all whitespace
    path = 'data/stds/%s.fits' % obname
    package = __name__.split('.')[0]
    full_path = pkg_resources.resource_filename(package, path)
    if os.path.exists(full_path):
        logger.info("Found std file: %s" % full_path)
        stdfile = full_path
        stdname = obname
    else:
        logger.info("Not found in data/stds: %s" % full_path)
    return stdfile, stdname


def is_file_kcwi_std(filename, logger=None):
    try:
        with fits.open(filename) as hdul:
            imtype = hdul[0].header['IMTYPE']
            targname = hdul[0].header['TARGNAME']
            if 'OBJECT' in imtype:
                stdfile, stdname = kcwi_get_std(targname, logger=logger)
                if stdfile is not None:
                    return True
        return False
    except FileNotFoundError:
        logger.error(f"Could not find file at {filename}")
        return False
