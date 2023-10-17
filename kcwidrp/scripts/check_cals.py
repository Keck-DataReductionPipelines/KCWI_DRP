"""
This script scans a given directory and returns a report on whether the cals
inside the directory meet the minimum requirements needed to reduce the science
frames found inside. It follows the following logic:

1. Finds all OBJECT frames. It uses these to determine what instrument
configurations are needed. This is saved in an internal "Proc table"
2. Searches for BIAS and CONTBARS frames for needed setups (1x1 and 2x2)
3. Searches for ARCS in the right wavelength range
4. Optionally, checks for STANDARDS
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
from kcwidrp.core.kcwi_get_std import is_file_kcwi_std



warnings.simplefilter('ignore', category=AstropyWarning)

def parse_args():

    parser = argparse.ArgumentParser(description="DESC")

    parser.add_argument('-d', '--data_dir', dest="data_dir", help="Directory with data in it")
    parser.add_argument('-f', '--file_pattern', dest="file_glob", default="*.fits")
    parser.add_argument('-v', '--verbose', dest="verbose", action="store_true")
    parser.add_argument('-s', '--silent', dest="silent", action="store_true")
    parser.add_argument('-c', '--config', dest="config", type=str, help="KCWI configuration file", default=None)

    return parser.parse_args()

def check_cal_type(ccd_frame, setup_frame, targ_type, minimum, logger):
    found_list = proctab.search_proctab(frame=ccd_frame, target_type=targ_type) #target_group=setup_frame['GRPID']
    logger.debug(f"For setup {setup_frame['CID']}, found {len(found_list)} {targ_type} frames, need {minimum}")
    if len(found_list) >= minimum:
        return "PASSED"
    else:
        return f"{len(found_list)}, NEEDED {minimum}"

if __name__ == "__main__":

    # Parse arguments
    args = parse_args()


    # Set up logging
    logger = logging.getLogger("Logger")
    logger.addHandler(logging.StreamHandler())

    if args.verbose:
        logger.setLevel("DEBUG")
    elif args.silent:
        logger.setLevel("WARNING")
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
    #pprint.pprint({section: dict(config[section]) for section in config.sections()})



    logger.debug(f"Looking for files matching {args.data_dir}/{args.file_glob}")
    data_dir = Path(args.data_dir)

    proc_logger = logging.getLogger("proc_logger")
    proctab = Proctab(proc_logger)
    proctab.read_proctab("temp.proc")


    # Load all files into a proctable
    files = list(data_dir.glob(args.file_glob))
    frames = {}
    standards = {}
    logger.info(f"Found {len(files)} files to process")
    for file in files:
        try:
            frame = CCDData.read(file, unit='adu')
        except FileNotFoundError as e:
            logger.error(f"Failed to open {file}")
        proctab.update_proctab(frame, filename=file.name)
        frames[Path(file).name] = frame

    # Get all the unique OBJECT setups:
    objects = proctab.proctab[proctab.proctab["TYPE"] == "OBJECT"]

    CIDs = set()
    DIDs = set()
    unique_setups = []
    for obj in objects:
        if obj["CID"] not in CIDs or obj["DID"] not in DIDs:
            CIDs.add(obj["CID"])
            DIDs.add(obj["DID"])
            unique_setups.append(obj)

    # Log results
    logger.info(f"Found {len(unique_setups)} individual setups")

    logger.debug("Found setups:")
    logger.debug(f"{'Config ID': <25}{'Detector ID': <15}{'Camera': <10}{'IFU': <10}{'Grating': <10}{'Grating Angle': <15}{'Central Wavelength': <20}{'Binning': <10}")
    for obj in unique_setups:
        logger.debug(f"{obj['CID']: <25}{obj['DID']: <15}{obj['CAM']: <10}{obj['IFU']: <10}{obj['GRAT']: <10}{obj['GANG']: <15}{obj['CWAVE']: <20}{obj['BIN']: <10}")

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
            "all_pass": False
            # "STANDARDS" : "UNCHECKED"
        }
        results["BIAS"] = check_cal_type(ccd_frame, setup_frame, "BIAS", bias_min_frames, logger)
        results["CONTBARS"] = check_cal_type(ccd_frame, setup_frame, "CONTBARS", cont_min_frames, logger)
        results["ARCS"] = check_cal_type(ccd_frame, setup_frame, "ARCLAMP", arc_min_frames, logger)
        results["FLATS"] = check_cal_type(ccd_frame, setup_frame, "FLATLAMP", flat_min_frames, logger)
        results["all_pass"] = (results["BIAS"] == "PASSED" and 
                            results["CONTBARS"] == "PASSED" and
                            results["ARCS"] == "PASSED" and
                            results["FLATS"] == "PASSED")
        
        # Check standards:
        # results["STANDARDS"] = ""
        report[setup_frame['CID']] = results

    # Print results, passes then fails
    passes = []
    fails = []
    for setup_id in report.keys():
        if report[setup_id]['all_pass']:
            passes.append(setup_id)
        else:
            fails.append(setup_id)
    if len(passes) > 0:
        logger.info("Following setups PASSED:")
        logger.info(", ".join(passes))
    if len(fails) > 0:
        logger.info("Following setups FAILED:")
        for fail in fails:
            logger.info(pprint.pprint(fail))

