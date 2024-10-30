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
VERSION = '1.2.0'
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
        "kcwi_masksky_ds9 = kcwidrp.scripts.kcwi_masksky_ds9:main",
        "smart_reduce_kcwi = kcwidrp.scripts.smart_reduce_kcwi:main",
        "modhead = kcwidrp.scripts.modhead:main",
        "start_kcwi_rti = kcwidrp.scripts.kcwi_rti:main",
        "wb = kcwidrp.scripts.wb:wb_main",
        "wr = kcwidrp.scripts.wr:wr_main",
        "check_cals = kcwidrp.scripts.check_cals:main"
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
      package_data={'kcwidrp': ['configs/*.cfg', 'data/*',
                                'data/extin/*', 'data/stds/*']},
      scripts=scripts,
      entry_points=entry_points,
      install_requires=['scikit-image~=0.24',
                        'astropy~=6.1.3',
                        'astroscrappy~=1.1.0',
                        'ccdproc~=2.4.2',
                        'numpy==1.26.4',
                        'scipy~=1.14.1',
                        'pyerfa',
                        'bokeh~=2.4.3',
                        'jinja2~=3.0.3',
                        'psutil~=5.7.0',
                        'pytest~=5.4.1',
                        'keckdrpframework',
                        'requests',
                        'pandas~=1.3.5',
                        'matplotlib~=3.9.2',
                        'ref_index~=1.0',
                        'pyregion~=2.0',
                        'cython',
                        'selenium',
                        'phantomjs'],
      python_requires="~=3.12"
      )
