==============================================
Keck Cosmic Web Images Data Reduction Pipeline
==============================================

This pipeline is in development and not ready for production.  It uses the KeckDRPFramework and will be a complete port of the IDL
pipeline to python.  For now use the IDL pipeline (KcwiDRP).

03-Jan-2020

Documentation
-------------

Documentation for this project is maintained with ReadTheDocs at the following link:

`https://kcwi-drp.readthedocs.io/en/latest/installing.html`


.. image:: https://readthedocs.org/projects/kcwi-drp/badge/?version=latest
   :target: https://kcwi-drp.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

Tests
-----

Status of tests:
`<https://travis-ci.com/github/Keck-DataReductionPipelines/KCWI_DRP>`_


.. image:: https://travis-ci.com/Keck-DataReductionPipelines/KCWI_DRP.svg?branch=master
   :target: https://travis-ci.com/Keck-DataReductionPipelines/KCWI_DRP


Dependencies
------------

* conda install bokeh
* conda install selenium
* conda install -c conda-forge phantomjs
* conda install -c astropy ccdproc
* pip install ref_index
* # LACOSMICX
* git clone https://github.com/cmccully/lacosmicx.git
* cd lacosmicx
* python setup.py install

It is possible to install the requirements by using:

.. code-block:: bash

   pip install -r requirements.txt

Also:

USE THE MASTER BRANCH of the KeckDRPFramework

Luca, March 6, 2020
