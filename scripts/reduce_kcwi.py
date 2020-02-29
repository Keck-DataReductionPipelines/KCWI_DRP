"""
Created on Jul 19, 2019

Test Fits to PNG pipeline with HTTP server.

@author: skwok
"""

from keckdrpframework.core.framework import Framework
from keckdrpframework.config.framework_config import ConfigClass
from keckdrpframework.models.arguments import Arguments
from kcwidrp.logger import start_logger

import subprocess
import time
import argparse
import sys
import traceback
import os
import pkg_resources
import threading

from kcwidrp.pipelines.kcwi_pipeline import Kcwi_pipeline
from kcwidrp.core.kcwi_proctab import Proctab


def _parseArguments(in_args: list) -> argparse.Namespace:
    description = "KCWI pipeline CLI"

    # this is a simple case where we provide a frame and a configuration file
    parser = argparse.ArgumentParser(prog=f"{in_args[0]}",
                                     description=description)
    parser.add_argument('-c', '--config', dest="kcwi_config_file", type=str,
                        help="KCWI configuration file", default=None)
    parser.add_argument('-f', '--frames', nargs='*', type=str,
                        help='input image files (full path, list ok)',
                        default=None)

    # in this case, we are loading an entire directory,
    # and ingesting all the files in that directory
    parser.add_argument('-i', '--infiles', dest="infiles", help="Input files",
                        nargs="*")
    parser.add_argument('-d', '--directory', dest="dirname", type=str,
                        help="Input directory", nargs='?', default=None)
    # after ingesting the files,
    # do we want to continue monitoring the directory?
    parser.add_argument('-m', '--monitor', dest="monitor",
                        help='Continue monitoring the directory '
                             'after the initial ingestion',
                        action='store_true', default=False)

    # special arguments, ignore
    parser.add_argument("-I", "--ingest_data_only", dest="ingest_data_only",
                        action="store_true",
                        help="Ingest data and terminate")
    parser.add_argument("-w", "--wait_for_event", dest="wait_for_event",
                        action="store_true", help="Wait for events")
    parser.add_argument("-W", "--continue", dest="continuous",
                        action="store_true",
                        help="Continue processing, wait for ever")
    parser.add_argument("-s", "--start_queue_manager_only",
                        dest="queue_manager_only", action="store_true",
                        help="Starts queue manager only, no processing",)

    # kcwi specific parameter
    parser.add_argument("-p", "--proctab", dest='proctab', help='Proctab file',
                        default='kcwi.proc')

    args = parser.parse_args(in_args[1:])
    return args


def check_redux_dir():
    if not os.path.isdir(framework.config.instrument.output_directory):
        os.makedirs(framework.config.instrument.output_directory)
        print("Output directory created: %s" %
              framework.config.instrument.output_directory)

def check_logs_dir():
    if not os.path.isdir('logs'):
        os.makedirs('logs')
        print("Logs directory created")


if __name__ == "__main__":

    args = _parseArguments(sys.argv)

    # configuration files
    pkg='kcwidrp'

    # check for the logs diretory
    check_logs_dir()

    framework_config_file = "configs/framework.cfg"
    framework_config_path = pkg_resources.resource_filename(pkg, framework_config_file)

    framework_logcfg_file = 'configs/framework_logger.cfg'
    framework_logcfg_path = pkg_resources.resource_filename(pkg, framework_logcfg_file)

    # add kcwi specific config files # make changes here to allow this file to be loaded from the command line
    if args.kcwi_config_file is None:
        kcwi_config_file = 'configs/kcwi.cfg'
        kcwi_config_path = pkg_resources.resource_filename(pkg, kcwi_config_file)
        kcwi_config = ConfigClass(kcwi_config_path)
    else:
        kcwi_config = ConfigClass(args.kcwi_config_file)

    try:
        framework = Framework(Kcwi_pipeline, framework_config_path)
        framework.logger = start_logger('DRPFrame', framework_logcfg_path)
        framework.config.instrument = kcwi_config
    except Exception as e:
        print("Failed to initialize framework, exiting ...", e)
        traceback.print_exc()
        sys.exit(1)


    framework.context.pipeline_logger = start_logger(pkg, kcwi_config_path)

    # check for the REDUX directory
    check_redux_dir()

    if framework.config.enable_bokeh is True:
        subprocess.Popen('bokeh serve', shell=True)
        time.sleep(5)


    # initialize the proctab and read it
    framework.context.proctab = Proctab(framework.logger)
    framework.context.proctab.read_proctab(tfil=args.proctab)

    framework.logger.info("Framework initialized")
    if framework.config.enable_bokeh:
        framework.append_event('start_bokeh', None)

    # start queue manager only (useful for RPC)
    if args.queue_manager_only:
        # The queue manager runs for ever.
        framework.logger.info("Starting queue manager only, no processing")
        framework.start_queue_manager()

    # single frame processing
    elif args.frames:
        for frame in args.frames:
            framework.ingest_data(None, [frame], False)
        for frame in args.frames:
            arguments = Arguments(name=frame)
            framework.append_event('next_file', arguments)

    # ingest an entire directory, trigger "next_file" on each file,
    # optionally continue to monitor if -m is specified
    elif (len(args.infiles) > 0) or args.dirname is not None:
        framework.ingest_data(args.dirname, args.infiles, args.monitor)
        # framework.context.data_set.start_monitor()

    framework.start(args.queue_manager_only, args.ingest_data_only,
                    args.wait_for_event, args.continuous)
