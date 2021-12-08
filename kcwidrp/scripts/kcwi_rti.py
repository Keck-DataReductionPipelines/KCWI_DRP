"""
Created on May 05, 2021

@author: skwok, mbrodheim
"""

from keckdrpframework.core.framework import Framework
from keckdrpframework.config.framework_config import ConfigClass
from keckdrpframework.models.arguments import Arguments
from keckdrpframework.utils.drpf_logger import getLogger
from kcwidrp.core.bokeh_plotting import check_running_process
from kcwidrp.core.kcwi_get_std import is_file_kcwi_std
from kcwidrp.core.notify_api import send_data_api

import subprocess
import time
import argparse
import sys
import traceback
import os
import pkg_resources
import types
import psutil
import getpass
import pandas

from kcwidrp.pipelines.keck_rti_pipeline import Keck_RTI_Pipeline
from kcwidrp.core.kcwi_proctab import Proctab
import logging.config


def _parse_arguments(in_args: list) -> argparse.Namespace:
    description = "KCWI pipeline CLI"

    # this is a simple case where we provide a frame and a configuration file
    parser = argparse.ArgumentParser(prog=f"{in_args[0]}",
                                     description=description)
    parser.add_argument('-c', '--config', dest="kcwi_config_file", type=str,
                        help="KCWI configuration file", default=None)
    parser.add_argument('-r', '--rti_config', dest="rti_config_file", type=str,
                        help="KCWI RTI configuration file", default=None)
    parser.add_argument('-f', '--frames', nargs='*', type=str,
                        help='input image files (full path, list ok)',
                        default=None)
    parser.add_argument('-l', '--list', dest='file_list',
                        help='File containing a list of files to be processed',
                        default=None)
    parser.add_argument('-g', '--groups', dest='group_mode',
                        help='Use group mode: separate files by image type and '
                             'reduce in the correct order',
                        default=False, action='store_true')
    parser.add_argument('-t', '--taperfrac', dest='taperfrac', type=float,
                        help='Taper fraction for wavelength fitting',
                        default=None)
    parser.add_argument('-a', '--atlas_line_list', dest='atlas_line_list',
                        type=str, help="Atlas line list file", default=None)

    # in this case, we are loading an entire directory,
    # and ingesting all the files in that directory
    parser.add_argument('-i', '--infiles', dest="infiles",
                        help="Input files, or pattern to match", nargs="?")
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
    parser.add_argument("--rti-level", dest="rti_lev", type=str,
                        choices=['lev1', 'lev2'], help="[lev1, lev2]")

    # kcwi specific parameter
    parser.add_argument("-p", "--proctab", dest='proctab', help='Proctab file',
                        default='kcwi.proc')

    parser.add_argument("--rti-level", dest="rti_lev", help="[lev1, lev2]")

    out_args = parser.parse_args(in_args[1:])
    return out_args


def main():

    args = _parse_arguments(sys.argv)

    def _on_exit(_, exit_status):

        if rti_config.rti_ingesttype == 'lev2':
            _send_complete()

        stop_processes(rti_config.extra_process, framework.logger)

        os._exit(exit_status)

    def _send_complete():
        data_date = None
        for direct in kcwi_config.cwd.split('/'):
            try:
                data_date = int(direct)
                assert len(str(data_date)) == 8
            except (ValueError, AssertionError):
                pass

        if not data_date:
            framework.logger.error("Could not determine the date to send "
                                   "the complete status to koa")
            return

        data = {'instrument': 'KCWI', 'utdate': data_date, 'ingesttype': 'lev2'}

        result = send_data_api(data, rti_config, framework.logger)
        framework.logger.info(f"Complete status sent to koa_status table, "
                              f"with result {result}")


    # START HANDLING OF CONFIGURATION FILES ##########
    pkg = 'kcwidrp'

    # check for the logs directory
    check_directory("logs")

    # check for the plots directory
    check_directory("plots")

    framework_config_file = "configs/framework.cfg"
    framework_config_fullpath = \
        pkg_resources.resource_filename(pkg, framework_config_file)

    framework_logcfg_file = 'configs/logger.cfg'
    framework_logcfg_fullpath = \
        pkg_resources.resource_filename(pkg, framework_logcfg_file)

    # add kcwi specific config files # make changes here to allow this file
    # to be loaded from the command line
    if not args.kcwi_config_file:
        kcwi_config_file = 'configs/kcwi_koarti.cfg'
        kcwi_config_fullpath = pkg_resources.resource_filename(
            pkg, kcwi_config_file)
    else:
        kcwi_config_fullpath = args.kcwi_config_file

    kcwi_config = ConfigClass(kcwi_config_fullpath, default_section='KCWI')

    # RTI specific configuration file
    if not args.rti_config_file:
        rti_config_file = "configs/rti.cfg"
        rti_config_fullpath = pkg_resources.resource_filename(
            pkg, rti_config_file)
    else:
        rti_config_fullpath = args.rti_config_file

    dfault_section = f'RTI_LEVEL{kcwi_config.level}'
    rti_config = ConfigClass(rti_config_fullpath, default_section=dfault_section)

    if args.rti_lev:
        rti_config.rti_ingesttype = args.rti_lev

    # END HANDLING OF CONFIGURATION FILES ##########

    # Add current working directory to config info
    kcwi_config.cwd = os.getcwd()

    # check for the output directory
    check_directory(kcwi_config.output_directory)

    try:
        framework = Framework(Keck_RTI_Pipeline, framework_config_fullpath)
        # add this line ONLY if you are using a local logging config file
        logging.config.fileConfig(framework_logcfg_fullpath)
        framework.config.instrument = kcwi_config
        framework.config.rti = rti_config
    except Exception as e:
        print("Failed to initialize framework, exiting ...", e)
        traceback.print_exc()
        sys.exit(1)

    framework.context.pipeline_logger = getLogger(framework_logcfg_fullpath,
                                                  name="KCWI")
    framework.logger = getLogger(framework_logcfg_fullpath, name="DRPF")

    if args.infiles is not None:
        framework.config.file_type = args.infiles

    # check for taperfrac argument
    if args.taperfrac:
        def_tf = getattr(framework.config.instrument, 'TAPERFRAC', None)
        if def_tf is not None:
            framework.context.pipeline_logger.info(
                "Setting new taperfrac = %.3f" % args.taperfrac)
            framework.config.instrument.TAPERFRAC = args.taperfrac

    # check for atlas line list argument
    if args.atlas_line_list:
        def_ll = getattr(framework.config.instrument, 'LINELIST', None)
        if def_ll is not None:
            framework.context.pipeline_logger.info(
                "Using line list %s instead of generated list" %
                args.atlas_line_list)
            framework.config.instrument.LINELIST = args.atlas_line_list

    # start the bokeh server is requested by the configuration parameters
    if framework.config.instrument.enable_bokeh:
        if not check_running_process(process='bokeh'):
            subprocess.Popen('bokeh serve', shell=True)
            # --session-ids=unsigned --session-token-expiration=86400',
            # shell=True)
            time.sleep(5)
        # subprocess.Popen('open http://localhost:5006?bokeh-session-id=kcwi',
        # shell=True)

    # initialize the proctab and read it
    framework.context.proctab = Proctab(framework.logger)
    framework.context.proctab.read_proctab(tfil=args.proctab)

    framework.logger.info("Framework initialized")
    framework.logger.info(f"RTI url is {framework.config.rti.rti_url}")

    # add a start_bokeh event to the processing queue,
    # if requested by the configuration parameters
    if framework.config.instrument.enable_bokeh:
        framework.append_event('start_bokeh', None)

    # important: to be able to use the grouping mode, we need to reset
    # the default ingestion action to no-event,
    # otherwise the system will automatically trigger the next_file event
    # which initiates the processing
    if args.group_mode:
        # set the default ingestion event to None
        framework.config.default_ingestion_event = "add_only"

    # start queue manager only (useful for RPC)
    if args.queue_manager_only:
        # The queue manager runs for ever.
        framework.logger.info("Starting queue manager only, no processing")
        framework.start_queue_manager()

    # in the next two ingest_data command, if we are using standard mode,
    # the first event generated is next_file.
    # if we are in groups mode (aka smart mode), then the first event generated
    # is no_event then, manually, a next_file event is generated for each group
    # specified in the variable imtypes

    # single frame processing
    elif args.frames:
        framework.ingest_data(None, args.frames, False)

    # processing of a list of files contained in a file
    elif args.file_list:
        frames = []
        with open(args.file_list) as file_list:
            for frame in file_list:
                if "#" not in frame:
                    frames.append(frame.strip('\n'))
        framework.ingest_data(None, frames, False)

    # ingest an entire directory, trigger "next_file" (which is an option
    # specified in the config file) on each file,
    # optionally continue to monitor if -m is specified
    elif args.dirname is not None:

        framework.ingest_data(args.dirname, None, args.monitor)

    # implement the group mode
    if args.group_mode:
        process_group_mode(framework, args)

    framework.on_exit = types.MethodType(_on_exit, framework)
    framework.start(args.queue_manager_only, args.ingest_data_only,
                    args.wait_for_event, args.continuous)


def check_directory(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)
        print("Directory %s has been created" % directory)


def stop_processes(extra_process, logger):
    current_user = getpass.getuser()

    if not extra_process:
        logger.info('No extra processes defined to stop.')
        return

    extra_process = extra_process.split(',')

    for proc in psutil.process_iter():
        pinfo = proc.as_dict(attrs=['name', 'username', 'pid', 'cmdline'])

        if pinfo['username'] != current_user:
            continue

        for cmd in pinfo['cmdline']:
            for process in extra_process:
                if process in cmd:
                    logger.info(f"Terminating process: {pinfo['cmdline']}, "
                                f"pid: {pinfo['pid']}")
                    p = psutil.Process(pinfo['pid'])
                    p.terminate()


def process_group_mode(framework, args):
    """
    Group the files by imagetype and binning.  Then process in sets and in the
    sequence required to reduce the next set.  The order is defined by imtypes
    list,  and the processing is done one set at a time.  Each set is processed
    together in parallel and then returns before starting the next set.
    """

    def _on_grp_complete(_, exit_status):
        return

    def add_grp(grp_list, resched):
        for in_frame in grp_list:
            framework.logger.info(f"File to process: {in_frame}")
            arguments = Arguments(name=in_frame, reschedule=resched)
            framework.append_event('next_file', arguments, recurrent=True)

    def process_grp(grp_list, resched):
        add_grp(grp_list, resched)

        framework.on_exit = types.MethodType(_on_grp_complete, framework)
        framework.start(args.queue_manager_only, args.ingest_data_only,
                        args.wait_for_event, args.continuous)

    def process_by_expt(exptimes, imtype, binning):
        for expt in exptimes:
            subset = data_set.data_table[
                (framework.context.data_set.data_table.IMTYPE == imtype)
                & (framework.context.data_set.data_table.BINNING.str.contains(binning))
                & (framework.context.data_set.data_table.TTIME == expt)]

            framework.logger.info(f"Processing type={imtype}, binning={binning}"
                                  f", exptime={expt}")

            process_grp(subset.index, False)

    def divide_and_process():
        if imtype == 'BIAS':
            subset = data_set.data_table[
                (framework.context.data_set.data_table.IMTYPE == imtype) &
                (framework.context.data_set.data_table.BINNING.str.contains(binning))
                & (framework.context.data_set.data_table.TTIME == 0.0)]
        else:
            subset = data_set.data_table[
                (framework.context.data_set.data_table.IMTYPE == imtype) &
                (framework.context.data_set.data_table.BINNING.str.contains(binning))]

        if len(subset) == 0:
            return

        framework.logger.info(f"Processing type={imtype}, binning={binning}")

        if 'OBJECT' in imtype:  # Ensure that standards are processed first
            for frame in subset.index:
                # this only checks for standards with same name as KCWI standards
                if is_file_kcwi_std(frame, logger=framework.context.logger):
                    standard.append(frame)
                else:
                    science.append(frame)

        elif 'DARK' in imtype:
            exptimes = set()
            for expt in subset.TTIME:
                exptimes.add(expt)
            process_by_expt(exptimes)

        else:
            process_grp(subset.index, False)

    data_set = framework.context.data_set

    # remove focus images and ccd clearing images from the dataset
    try:
        focus_frames = data_set.data_table[data_set.data_table.OBJECT == "focus"].index
    except AttributeError:
        framework.logger.warn("Empty dataset,  exiting group mode.")
        return

    ccdclear_frames = data_set.data_table[data_set.data_table.OBJECT == "Clearing ccd"].index
    data_set.data_table.drop(focus_frames, inplace=True)
    data_set.data_table.drop(ccdclear_frames, inplace=True)

    # list order is important to process the calibrations in the correct order.
    imtypes = ['BIAS', 'DARK', 'CONTBARS', 'ARCLAMP', 'FLATLAMP',
               'DOMEFLAT', 'TWIFLAT', 'OBJECT']

    allowed_binning = ('1,1', '2,2')

    science = []
    standard = []

    for imtype in imtypes:
        for binning in allowed_binning:
            divide_and_process()

    process_grp(standard, False)
    add_grp(science, False)


if __name__ == "__main__":
    main()
