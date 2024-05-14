#!/usr/bin/env python

"""
This script scans a given directory and returns a report on whether the cals
inside the directory meet the minimum requirements needed to reduce the science
frames found inside. It follows the following logic:

#. Finds all OBJECT frames. It uses these to determine what instrument configurations are needed. This is saved in an internal "Proc table"

#. Searches for BIAS and CONTBARS frames for needed setups (1x1 and 2x2)

#. Searches for ARCS and FLATS

#. Searches for matching standard stars for the setup

"""

from pathlib import Path
import argparse
import logging
import warnings
import pkg_resources

from astropy.nddata import CCDData
from astropy.utils.exceptions import AstropyWarning

from kcwidrp.core.kcwi_proctab import Proctab
from keckdrpframework.config.framework_config import ConfigClass
from kcwidrp.core.kcwi_get_std import kcwi_get_std
from kcwidrp.primitives.kcwi_file_primitives import fix_header


warnings.simplefilter('ignore', category=AstropyWarning)


def parse_args():
    """Parse arguments passed into this script from the command line

    Returns
    -------
    argparse
        Dict-like object with parsed args

    :meta private:
    """

    parser = argparse.ArgumentParser(
        description="Checks a directory of data to see if the OBJECTs within "
                    "have the required cals.")

    parser.add_argument('filepaths', help="Files to inspect", nargs="+")
    parser.add_argument('-v', '--verbose', dest="verbose", action="store_true",
                        help="Print exhaustive information")
    parser.add_argument('-a', '--auto', dest="auto", action="store_true",
                        help='For autonomy mode, only returns PASSED or FAILED.')
    parser.add_argument('-c', '--config', dest="config", type=str,
                        help="KCWI configuration file", default=None)

    return parser.parse_args()


def check_cal_type(proctab, ccd_frame, setup_frame, targ_type, minimum, logger):
    """Checks the given proctab for calibrations matching the input frame

    Parameters
    ----------
    proctab : Proctab
        Proccessing table instance to search through
    ccd_frame : CCDData
        Astropy CCDData (or kcwidrp KCCDData) object we want to match cals to
    setup_frame : astropy.Table row
        Proctab row corresponding to an OBJECT setup we are looking for
    targ_type : str
        Type of calibration being looked for
    minimum : int
        Minimum number of this type of cal required
    logger : logging.Logger
        Logger for debuging

    Returns
    -------
    str
        String representing the results of the search. PASSED if enough matching cals were found, FAILED otherwise
    
    :meta private:
    """

    found_list = proctab.search_proctab(frame=ccd_frame, target_type=targ_type) #target_group=setup_frame['GRPID']
    logger.debug(f"For setup {setup_frame['CID']}, found {len(found_list)} {targ_type} frames, need {minimum}")
    if len(found_list) >= minimum:
        return "PASSED"
    else:
        return f"FAILED: found {len(found_list): <3}, needed {minimum: <3}"

def main():

    # Parse arguments
    args = parse_args()


    # Set up logging
    logger = logging.getLogger("Logger")
    logger.addHandler(logging.StreamHandler())

    if args.auto:
        logger.setLevel("ERROR")
    else:
        if args.verbose:
            logger.setLevel("DEBUG")
        else:
            logger.setLevel("INFO")
    

    # Load config
    pkg = 'kcwidrp'
    if args.config is None:
        kcwi_config_file = 'configs/kcwi.cfg'
        kcwi_config_fullpath = pkg_resources.resource_filename(pkg, kcwi_config_file)
        config = ConfigClass(kcwi_config_fullpath, default_section='KCWI')
    else:
        config = ConfigClass(args.config, default_section='KCWI')

    logger.debug(f"Loaded KCWI config")

    # Prepare to read files
    files = [Path(file) for file in args.filepaths]
    if len(files) == 0:
        logger.error("No files found! Exiting...")
        return
    logger.info(f"Found {len(files)} files to inspect")


    # Set up proc table. Use a dummy logger, since we don't want its messages
    dummy_logger = logging.getLogger("dummy_logger")
    proctab = Proctab(dummy_logger)
    proctab.read_proctab("temp.proc")


    # Load all files into the proctable
    frames = {}
    for file in files:
        try:
            frame = CCDData.read(file, unit='adu')
            fix_header(frame)
            if 'CCDCFG' not in frame.header:
                ccdcfg = frame.header['CCDSUM'].replace(" ", "")
                ccdcfg += "%1d" % frame.header['CCDMODE']
                ccdcfg += "%02d" % frame.header['GAINMUL']
                ccdcfg += "%02d" % frame.header['AMPMNUM']
                frame.header['CCDCFG'] = ccdcfg
        except FileNotFoundError as e:
            logger.error(f"Failed to open {file}")
        proctab.update_proctab(frame, filename=file.name)
        frames[Path(file).name] = frame

    # Get all the unique OBJECT setups:
    objects = proctab.proctab[proctab.proctab["TYPE"] == "OBJECT"]

    standards = {}
    CIDs = set()
    DIDs = set()
    unique_setups = []
    for obj in objects:
        # See if the setup is unique
        if obj["CID"] not in CIDs or obj["DID"] not in DIDs:
            CIDs.add(obj["CID"])
            DIDs.add(obj["DID"])
            unique_setups.append(obj)
        
        # While we're here, compile a dict with all standard stars ordered by CID for later
        stdfile, _ = kcwi_get_std(obj["TARGNAME"], logger=dummy_logger)
        if stdfile is not None:
            CID = standards.get(obj["CID"], None)
            if CID is None:
                standards[obj["CID"]] = {"files" : [], "names" : []}
                CID = standards[obj["CID"]]
            CID["files"].append(obj['filename'])
            CID["names"].append(obj['TARGNAME'])

    # Log results
    logger.info(f"\nFound {len(unique_setups)} individual setups:")

    logger.info(f"{'Config ID': <25}{'Detector ID': <15}{'Camera': <10}{'IFU': <10}{'Grating': <10}{'Grating Angle': <15}{'Central Wavelength': <20}{'Binning': <10}{'Standards': <25}")
    for obj in unique_setups:
        matching_standards = standards.get(obj["CID"], None)
        if matching_standards:
            standards_str = ", ".join(set(matching_standards["names"]))
        else:
            standards_str = "NA"
        logger.info(f"{obj['CID']: <25}{obj['DID']: <15}{obj['CAM']: <10}{obj['IFU']: <10}{obj['GRAT']: <10}{obj['GANG']: <15}{obj['CWAVE']: <20}{obj['BIN']: <10}{standards_str: <25}")


    # Check cals for each setup:
    report = {}
    for setup_frame in unique_setups:
        frame_config = config[setup_frame['CAM']]

        bias_min_frames = int(frame_config['bias_min_nframes'])
        cont_min_frames = int(frame_config['contbars_min_nframes'])
        arc_min_frames = int(frame_config['arc_min_nframes'])
        flat_min_frames = int(frame_config['flat_min_nframes'])


        ccd_frame = frames[setup_frame['filename']]
        results = {
            "BIAS" : "UNCHECKED",
            "CONTBARS" : "UNCHECKED",
            "ARCS" : "UNCHECKED",
            "FLATS" : "UNCHECKED",
            "all_pass": False,
            "STANDARDS" : "UNCHECKED"
        }
        results["BIAS"] = check_cal_type(proctab, ccd_frame, setup_frame, "BIAS", bias_min_frames, logger)
        results["CONTBARS"] = check_cal_type(proctab, ccd_frame, setup_frame, "CONTBARS", cont_min_frames, logger)
        results["ARCS"] = check_cal_type(proctab, ccd_frame, setup_frame, "ARCLAMP", arc_min_frames, logger)
        results["FLATS"] = check_cal_type(proctab, ccd_frame, setup_frame, "FLATLAMP", flat_min_frames, logger)
        
        # Check the standards list from earlier for matching setups
        matching_standards = standards.get(setup_frame["CID"], None)
        if matching_standards is not None:
            out = {}
            for (name, file) in zip(matching_standards["names"], matching_standards["files"]):
                if out.get(name, None):
                    out[name].append(file)
                else:
                    out[name] = [file]
            out_str = ""
            for name in out.keys():
                out_str += f"{name} ({', '.join(out[name])})"
            # results["STANDARDS"] = ", ".join(matching_standards["names"])
            results["STANDARDS"] = out_str
            standards_result = True
        else:
            results["STANDARDS"] = "FAILED"
            standards_result = False
        
        # Collate a final "are we good here" boolean for the report
        results["all_pass"] = bool(results["BIAS"] == "PASSED" and 
                            results["CONTBARS"] == "PASSED" and
                            results["ARCS"] == "PASSED" and
                            results["FLATS"] == "PASSED" and
                            standards_result)
        
        report[setup_frame['CID']] = results


    # Print results, passes then fails
    passes = []
    fails = []
    logger.info("")
    
    # Collect the passes and fails into their own lists
    for setup_id in report.keys():
        if report[setup_id]['all_pass']:
            passes.append(setup_id)
        else:
            fails.append(setup_id)
    
    # Log info on passes (in green text)
    if len(passes) > 0:
        logger.info("Following setups \033[32mPASSED\033[0m:")
        logger.info(", ".join(passes))
        if args.verbose:
            for p in passes:
                logger.debug(f"Setup {p} will use the following standard stars:")
                logger.debug(report[p]["STANDARDS"])
    
    # Log info on all the fails (in red text)
    if len(fails) > 0:
        logger.info("Following setups \033[31mFAILED\033[0m:")
        for fail in fails:
            logger.info(f"{fail}:")
            for key in report[fail].keys():
                if key == "all_pass" : continue
                logger.info(f"\t{key: <10}\t{report[fail][key]: <20}")
            logger.info("\n\tThis effects the following OBJECT frames:")
            logger.info('\t' + str(objects[objects["CID"] == fail]['filename', 'TARGNAME']).replace('\n', '\n\t'))
    else:
        logger.info("\033[32mNo failures to report.\033[0m")
    
    if args.auto:
        if len(passes) > 0 and len(fails) == 0:
            print("PASSED")
        else:
            print("FAILED")

if __name__ == "__main__":
    main()