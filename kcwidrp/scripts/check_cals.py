#!/usr/bin/env python

"""
This script scans a given directory and returns a report on whether the cals
inside the directory meet the minimum requirements needed to reduce the science
frames found inside. It follows the following logic:

1. Finds all OBJECT frames. It uses these to determine what instrument
configurations are needed. This is saved in an internal "Proc table"
2. Searches for BIAS and CONTBARS frames for needed setups (1x1 and 2x2)
3. Searches for ARCS in the right wavelength range
4. Optionally, checks for STANDARDS -- TODO
"""

from pathlib import Path
import argparse
import logging
import warnings
import pkg_resources
import pprint

from astropy.nddata import CCDData
from astropy.utils.exceptions import AstropyWarning

from kcwidrp.core.kcwi_proctab import Proctab
from keckdrpframework.config.framework_config import ConfigClass
from kcwidrp.core.kcwi_get_std import kcwi_get_std



warnings.simplefilter('ignore', category=AstropyWarning)

def parse_args():

    parser = argparse.ArgumentParser(description="Checks a directory of data to see if the OBJECTs within have the required cals.")

    parser.add_argument('-d', '--data_dir', dest="data_dir", help="Directory with data in it")
    parser.add_argument('-f', '--file_pattern', dest="file_glob", default="*.fits", help="File pattern to match (e.g. KB*.fits)")
    parser.add_argument('-v', '--verbose', dest="verbose", action="store_true", help="Print exhaustive information")
    parser.add_argument('-c', '--config', dest="config", type=str, help="KCWI configuration file", default=None)

    return parser.parse_args()

def check_cal_type(proctab, ccd_frame, setup_frame, targ_type, minimum, logger):
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


    # Find files
    logger.debug(f"Looking for files matching {args.data_dir}/{args.file_glob}")
    data_dir = Path(args.data_dir)

    # Set up proc table. Use a dummy logger, since we don't want its messages
    dummy_logger = logging.getLogger("dummy_logger")
    proctab = Proctab(dummy_logger)
    proctab.read_proctab("temp.proc")


    # Load all files into the proctable
    files = list(data_dir.glob(args.file_glob))
    frames = {}
    logger.info(f"Found {len(files)} files to inspect")
    for file in files:
        try:
            frame = CCDData.read(file, unit='adu')
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
            results["STANDARDS"] = ", ".join(matching_standards["names"])
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
            logger.info(f"This effects the following OBJECT frames:")
            logger.info(objects[objects["CID"] == fail]['filename', 'TARGNAME'])
    else:
        logger.info("\033[32mNo failures to report.\033[0m:")

if __name__ == "__main__":
    main()