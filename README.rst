==============================================
Keck Cosmic Web Imager Data Reduction Pipeline
==============================================

This pipeline is the official DRP for the Keck Cosmic Web Imager at the W. M. 
Keck Observatory. It is built on the 
`KeckDRPFramework <https://github.com/Keck-DataReductionPipelines/KeckDRPFramework>`_,
and is a port of the old 
`IDL pipeline <https://github.com/Keck-DataReductionPipelines/KcwiDRP>`_.

The pipeline has only been thoroughly tested on MacOS, Ubuntu, and CentOS, although it
should run smoothly on any unix-like system that runs Python. Unfortunately, Windows 
is not currently supported.

01-Apr-2021

Documentation
-------------

Documentation for this project is maintained with ReadTheDocs at the following link:

`https://kcwi-drp.readthedocs.io`


.. image:: https://readthedocs.org/projects/kcwi-drp/badge/?version=latest
   :target: https://kcwi-drp.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

Tests
-----

The tests are run locally; see ``docs/testing.md`` for how to run them.


Dependencies
------------

Please see 'Installing KCWIDRP' in the documentation above for details about
dependencies and how to satisfy them. 



Developed by Luca Rizzi, Don Neill, Max Brodheim (04/2021)
