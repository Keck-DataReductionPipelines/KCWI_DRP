Keck Cosmic Web Imager Data Reduction Pipeline
==============================================

KCWI_DRP is a the official data reduction pipeline (DRP) for the Keck Cosmic Web Imager.
This DRP has been developed by the KCWI team at Caltech in
collaboration with the W. M. Keck Observatory Scientific Software Group.

Release 1.2
===========

We are excited to announce the release version of 1.2.0 of the KCWI DRP!

Like version 1.1, this version now with support for both the red and blue sides
of the instrument. This pipeline is based off of the original KCWIdrp IDL pipeline
developed by the KCWI team at Caltech, and this pipeline continues to be the only
pipeline supported by the W. M. Keck Observatory.


New Features
------------

* Simplified installation via pip and conda environment
* Vacuum to air and heliocentric or barycentric correction
* Ability of using KOA file names or original file names
* Full traceability of DRP versions and execution steps in HISTORY headers
* Versatile sky subtraction modes including using external sky frames, ability to manually mask regions, or skipping it entirely
* Formal support system via GitHub issues
* Added support for Mac intel M2 and M3 chips in Python 3.12

For older versions, see :doc:`versions`.

Users
=====

Follow the documentation in this section to reduce KCWI data.

.. toctree::
   :maxdepth: 1

   installing
   quickstart
   configuration
   running
   sky_subtraction
   calibration
   data_products
   scripts
   help
   versions


More information
================

.. toctree::
   :maxdepth: 1

   pipeline_concepts
   api

Developers
==========

The DRP was developed in collaboration between:

* Don Neill, Matt Matuszewki, Nik Prusinski, Chris Martin and the KCWI Team at Caltech
* Max Brodheim, Luca Rizzi, and Rosalie McGurk at W. M. Keck Observatory

Vacuum-to-air and helio/barycentric correction algorithms provided by Yuguang Chen at Caltech.

Maintenance and support of the DRP are provided as part of the Data Services Initiative (PI: John O'Meara), in collaboration with the National Aeronautics and Space Administration.

.. toctree::
   :maxdepth: 1

   updating_docs
   development
   issues


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
