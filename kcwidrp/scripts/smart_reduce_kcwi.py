"""
Created on Jul 19, 2019

Test Fits to PNG pipeline with HTTP server.

@author: skwok
"""

from keckdrpframework.core.framework import Framework
from keckdrpframework.config.framework_config import ConfigClass
from keckdrpframework.models.arguments import Arguments
from keckdrpframework.utils.drpf_logger import getLogger
from kcwidrp.core.bokeh_plotting import check_bokeh_server

import subprocess
import time
import argparse
import sys
import traceback
import os
import pkg_resources

from kcwidrp.pipelines.kcwi_pipeline import Kcwi_pipeline
from kcwidrp.core.kcwi_proctab import Proctab
import logging.config


def _parse_arguments(in_args: list) -> argparse.Namespace:
    description = "KCWI pipeline CLI"

    # this is a simple case where we provide a frame and a configuration file
    parser = argparse.ArgumentParser(prog=f"{in_args[0]}",
                                     description=description)
    parser.add_argument('-c', '--config', dest="kcwi_config_file", type=str,
                        help="KCWI configuration file", default=None)
    parser.add_argument('-f', '--frames', nargs='*', type=str,
                        help='input image files (full path, list ok)',
                        default=None)
    parser.add_argument('-l', '--list', dest='file_list',
                        help='File containing a list of files to be processed',
                        default=None)

    # kcwi specific parameter
    parser.add_argument("-p", "--proctab", dest='proctab', help='Proctab file',
                        default='kcwi.proc')

    out_args = parser.parse_args(in_args[1:])
    return out_args


def check_directory(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)
        print("Directory %s has been created" % directory)

def main():

    def process(subset):
        for frame in subset.index:
            arguments = Arguments(name=frame)
            framework.append_event('next_file', arguments)

    args = _parse_arguments(sys.argv)

    # START HANDLING OF CONFIGURATION FILES ##########
    pkg = 'kcwidrp'

    # check for the logs diretory
    check_directory("logs")
    # check for the plots directory
    check_directory("plots")

    framework_config_file = "configs/framework.cfg"
    framework_config_fullpath = pkg_resources.resource_filename(pkg, framework_config_file)

    framework_logcfg_file = 'configs/logger.cfg'
    framework_logcfg_fullpath = pkg_resources.resource_filename(pkg, framework_logcfg_file)

    # add kcwi specific config files # make changes here to allow this file
    # to be loaded from the command line
    if args.kcwi_config_file is None:
        kcwi_config_file = 'configs/kcwi.cfg'
        kcwi_config_fullpath = pkg_resources.resource_filename(pkg, kcwi_config_file)
        kcwi_config = ConfigClass(kcwi_config_fullpath, default_section='KCWI')
    else:
        kcwi_config_fullpath = os.path.abspath(args.kcwi_config_file)
        kcwi_config = ConfigClass(args.kcwi_config_file, default_section='KCWI')


    # END HANDLING OF CONFIGURATION FILES ##########

    try:
        framework = Framework(Kcwi_pipeline, framework_config_fullpath)
        # add this line ONLY if you are using a local logging config file
        logging.config.fileConfig(framework_logcfg_fullpath)
        framework.config.instrument = kcwi_config
    except Exception as e:
        print("Failed to initialize framework, exiting ...", e)
        traceback.print_exc()
        sys.exit(1)
    framework.context.pipeline_logger = getLogger(framework_logcfg_fullpath, name="KCWI")

    # check for the REDUX directory
    check_directory(framework.config.instrument.output_directory)

    # start the bokeh server is requested by the configuration parameters
    if framework.config.instrument.enable_bokeh is True:
        if check_bokeh_server() is False:
            subprocess.Popen('bokeh serve --session-ids=unsigned --session-token-expiration=86400', shell=True)
            time.sleep(5)
        subprocess.Popen('open http://localhost:5006?bokeh-session-id=kcwi', shell=True)


    # initialize the proctab and read it
    framework.context.proctab = Proctab(framework.logger)
    framework.context.proctab.read_proctab(tfil=args.proctab)

    framework.logger.info("Framework initialized")

    # add a start_bokeh event to the processing queue, if requested by the configuration parameters
    if framework.config.instrument.enable_bokeh:
        framework.append_event('start_bokeh', None)

    # set the default ingestion event to None
    framework.config.default_ingestion_event = "no_event"

    # single frame processing
    if args.frames:
        framework.ingest_data(None, args.frames, False)

    # processing of a list of files contained in a file
    elif args.file_list:
        frames = []
        with open(args.file_list) as file_list:
            for frame in file_list:
                if "#" not in frame:
                    frames.append(frame.strip('\n'))
        framework.ingest_data(None, frames, False)

    data_set = framework.context.data_set

    # remove focus images and ccd clearing images from the dataset
    focus_frames = data_set.data_table[data_set.data_table.OBJECT == "focus"].index
    ccdclear_frames = data_set.data_table[data_set.data_table.OBJECT == "Clearing ccd"].index
    data_set.data_table.drop(focus_frames, inplace=True)
    data_set.data_table.drop(ccdclear_frames, inplace=True)

    # processing
    imtypes = ['BIAS', 'CONTBARS', 'ARCLAMP', 'FLATLAMP', 'OBJECT']

    for imtype in imtypes:
        subset = data_set.data_table[framework.context.data_set.data_table.IMTYPE == imtype]
        process(subset)

    framework.start(False, False, True, True)

if __name__ == "__main__":
    main()