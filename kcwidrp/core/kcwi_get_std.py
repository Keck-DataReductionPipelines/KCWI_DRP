import os
import pkg_resources


def kcwi_get_std(targname, logger=None):
    """Checks if object is a standard star, returns std file path and name"""

    stdfile = None
    stdname = None
    obname = targname.lower()
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
