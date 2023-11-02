#!/usr/bin/env python
"""Modify FITS header keyword"""

from astropy.io import fits
import argparse
import sys
import os


def _parse_arguments(in_args: list) -> argparse.Namespace:

    description = "modify FITS header keyword"
    parser = argparse.ArgumentParser(prog=f"{in_args[0]}",
                                     description=description)

    parser.add_argument('infile', type=str, default=None,
                        help="Input FITS file")
    parser.add_argument('-k', '--keyword', type=str, default=None,
                        help="Keyword to modify")
    parser.add_argument('-v', '--value', type=str, default=None,
                        help="Value for keyword")
    parser.add_argument('-s', '--safe', action='store_true', default=False,
                        help="Keep original file")

    out_args = parser.parse_args(in_args[1:])
    return out_args


def main():
    args = _parse_arguments(sys.argv)

    if args.safe:
        safe_file = args.infile + '~'
        os.rename(args.infile, safe_file)
        with fits.open(safe_file) as hdul:
            hdul[0].header[args.keyword] = args.value
            hdul.writeto(args.infile)
    else:
        fits.setval(args.infile, args.keyword, value=args.value)


if __name__ == "__main__":
    main()
