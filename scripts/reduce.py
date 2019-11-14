'''
Created on Jul 19, 2019

Test Fits to PNG pipeline with HTTP server.

@author: skwok
'''

from keckdrpframework.core.framework import Framework
from keckdrpframework.config.framework_config import ConfigClass
from keckdrpframework.models.arguments import Arguments
import subprocess
import time
import argparse
import sys

from kcwidrp.pipelines.kcwi_pipeline import Kcwi_pipeline


def _parseArguments(in_args):
    description = "KCWI pipeline CLI"

    parser = argparse.ArgumentParser(prog=f"{in_args[0]}", description=description)
    parser.add_argument('-c', dest="config_file", type=str, help="Configuration file")
    parser.add_argument('frame', nargs=1, type=str, help='input image file')
    # parser.add_argument("-d", "--directory", dest="dirname", type=str, help="Input directory")

    args = parser.parse_args(in_args[1:])
    return args




if __name__ == "__main__":

    args = _parseArguments(sys.argv)

    config = ConfigClass(args.config_file)
    if config.enable_bokeh is True:
        subprocess.Popen('bokeh serve', shell=True)
        time.sleep(2)
    try:
        framework = Framework(Kcwi_pipeline, config)
    except Exception as e:
        print("Failed to initialize framework, exiting ...", e)
        traceback.print_exc()
        sys.exit(1)

    framework.logger.info("Framework initialized")
    arguments = Arguments(name=args.frame[0])
    print(arguments)
    framework.append_event('next_file', arguments)

    framework.start()

