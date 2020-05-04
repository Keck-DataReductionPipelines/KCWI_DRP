# !usr/bin/env python
# -*- coding: utf-8 -*-
#
# Licensed under a 3-clause BSD license.

from setuptools import setup

# Get some values from the setup.cfg
try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser

conf = ConfigParser()
conf.read(['setup.cfg'])
# metadata = dict(conf.items('metadata'))

NAME = 'kcwidrp'
VERSION = '0.11.1dev'
RELEASE = 'dev' not in VERSION

# scripts = [fname for fname in glob.glob(os.path.join('scripts', '*'))
#            if os.path.basename(fname) != 'README.rst']
scripts = []
# Define entry points for command-line scripts
entry_points = {
    'console_scripts': [
        "reduce_kcwi = kcwidrp.scripts.reduce_kcwi:main",
        "kcwi_masksky_ds9 = kcwidrp.scripts.kcwi_masksky_ds9:main"
        # "smart_reduce_kcwi = kcwidrp.scripts.smart_reduce_kcwi:main"
    ]}

setup(name=NAME,
      provides=NAME,
      version=VERSION,
      license='BSD3',
      description='KCWI DRP',
      long_description=open('README.txt').read(),
      author='Don Neill',
      author_email='neilljd@gmail.com',
      packages=['kcwidrp', ],
      scripts=scripts,
      entry_points=entry_points, install_requires=['ccdproc', 'bokeh', 'numpy',
                                                   'scipy', 'astropy']
      )
