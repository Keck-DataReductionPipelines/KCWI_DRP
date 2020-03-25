'''
Created on Jul 19, 2019

Test Fits to PNG pipeline with HTTP server.

@author: skwok
'''
import json
import socket

from keckdrpframework.utils.easyHTTP import EasyHTTPHandler, EasyHTTPServer, EasyHTTPServerThreaded
#from utils.try_wrapper import tryEx

import traceback

import sys
#sys.path.append('/Users/lrizzi/Python_Projects/Framework/prototype')
import os.path
import glob

from keckdrpframework.core.framework import Framework
from keckdrpframework.config.framework_config import ConfigClass
from keckdrpframework import config
from keckdrpframework.pipelines.kcwi_pipeline import Kcwi_pipeline
from keckdrpframework.models.arguments import Arguments
import subprocess
import time

import argparse
import requests
import os


def reduce(file_name):

    dirname = os.path.dirname(file_name)
    if dirname == '':
        cwd = os.getcwd()
        file_name = os.path.join(cwd, file_name)
    r = requests.get('http://127.0.0.1:50100/add_next_file_event?file_name=%s' % file_name)
    

def main():

    parser = argparse.ArgumentParser(description='Process a single file.')
    parser.add_argument('frame', nargs=1, type=str, help='input image file')

    args = parser.parse_args()

    reduce(args.frame[0])

if __name__ == "__main__":
    main()
