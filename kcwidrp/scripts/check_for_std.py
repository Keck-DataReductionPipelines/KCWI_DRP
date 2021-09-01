"""Script that takes an input list of fits files, and prints a list of which
ones exist in the standard star list
"""

from astropy.io import fits
import argparse
import sys

from kcwidrp.core.kcwi_get_std import is_file_kcwi_std as chk_std


def parse_arguments(in_args: list) -> argparse.Namespace:
    description = "KCWI pipeline CLI"

    # this is a simple case where we provide a frame and a configuration file
    parser = argparse.ArgumentParser(prog=f"{in_args[0]}",
                                     description=description)

    parser.add_argument("infiles", nargs="*", help="Files to be checked")

    out_args = parser.parse_args(in_args[1:])
    return out_args

class log():
    """Dummy class to silence the logger in is_file_kcwi_std
    """
    def info(*args):
        pass
    def error(*args):
        pass

def main():

    args = parse_arguments(sys.argv)
    
    for file in args.infiles:
        if chk_std(file, logger=log()):
            print(f"{file} is a standard star observation")