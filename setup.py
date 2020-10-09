# !usr/bin/env python
# -*- coding: utf-8 -*-
#
# Licensed under a 3-clause BSD license.

from setuptools import setup, find_packages

# Get some values from the setup.cfg
try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser

conf = ConfigParser()
conf.read(['setup.cfg'])
metadata = dict(conf.items("metadata"))
options = dict(conf.items("options"))

NAME = 'kcwidrp'
VERSION = '0.11.1dev'
RELEASE = 'dev' not in VERSION
AUTHOR = metadata["author"]
AUTHOR_EMAIL = metadata["author_email"]
LICENSE = metadata["license"]
DESCRIPTION = metadata["description"]

# scripts = [fname for fname in glob.glob(os.path.join('scripts', '*'))
#            if os.path.basename(fname) != 'README.rst']
scripts = []
# Define entry points for command-line scripts
entry_points = {
    'console_scripts': [
        "reduce_kcwi = kcwidrp.scripts.reduce_kcwi:main",
        "kcwi_masksky_ds9 = kcwidrp.scripts.kcwi_masksky_ds9:main"
        "smart_reduce_kcwi = kcwidrp.scripts.smart_reduce_kcwi:main"
    ]}

setup(name=NAME,
      provides=NAME,
      version=VERSION,
      license=LICENSE,
      description=DESCRIPTION,
      long_description=open('README.rst').read(),
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      packages=find_packages(),
      scripts=scripts,
      entry_points=entry_points,
      install_requires=['ccdproc', 'bokeh', 'numpy', 'scipy', 'astropy']
      )
