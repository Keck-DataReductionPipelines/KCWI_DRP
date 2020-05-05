"""Creates mask image from ds9 region file.

Args:
    imagename (string): The name of a *_intf.fits image
    regionname (string): The name of a ds9 region file

Returns:
    None

Note:
    To use this routine, process your data with kcwi_stage4flat.pro.
    Then display the target *_intf.fits file in ds9.
    Use region shapes to indicate non-sky pixels in image (box, circle, etc.).
    Write out ds9 region file (*.reg).
    Run this routine:

    python ~/kderp/devel/kcwi_masksky_ds9.py kb180101_00111_intf.fits ds9.reg

    (replace paths/filenames with your local paths/filenames)

    This should create kb180101_00111_smsk.fits, which will be used when you
    run kcwi_stage5sky.pro.
"""
from astropy.io import fits as pf
import pyregion

import numpy as np
import sys
import os


def main():
    # check args
    narg = len(sys.argv)

    # should be three (including routine name)
    if narg != 3:
        print("Usage: kcwi_masksky_ds9 <imagename> <regionname>")
        print("imagename : used for array dimensions and filename purposes, ")
        print("            must be an _intf image.")
        print("regionname: name of region file containing ds9 mask regions")
        print("            (typically a .reg)")
        exit()

    # read arg values
    imfname = sys.argv[1]
    regfname = sys.argv[2]

    # make sure it's an _intf image
    if '_intf.fits' not in imfname:
        print("imagename must be _intf.fits image")
        exit()

    # do inputs exist?
    if not os.path.exists(imfname):
        print("Image file does not exist: "+imfname)
        exit()

    if not os.path.exists(regfname):
        print("Region file does not exist: "+regfname)
        exit()

    # create output mask image name
    outfile = imfname.replace("_intf", "_smsk")
    print("Creating: "+outfile)

    # load the header from the pointed-to image.
    hdu_list = pf.open(imfname)
    header = hdu_list[0].header

    # load in the region file
    r = pyregion.open(regfname).as_imagecoord(header)
    m = pyregion.get_mask(r, hdu_list[0])

    # write out the mask
    hdu = pf.PrimaryHDU(np.uint8(m))
    hdu.writeto(outfile, overwrite=True)

    print("Done.")


if __name__ == "__main__":
    main()
