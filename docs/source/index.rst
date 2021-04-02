.. Keck Cosmic Web Imager DRP documentation master file, created by
   sphinx-quickstart on Fri May 29 10:09:02 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Keck Cosmic Web Imager DRP's
============================

KCWI_DRP is a the official data reduction pipeline (DRP) for the Keck Cosmic Web Imager.
This DRP has been developed by the KCWI team at Caltech in
collaboration with the W. M. Keck Observatory Scientific Software Group.

Release 1.0
===========

We are happy to announce the first official release of the KCWI pipeline. 
This pipeline is based entirely on the glorious IDL pipeline developed by the KCWI team at Caltech and it is 
the only pipeline supported by the W. M. Keck Observatory.

What this version provides
--------------------------

* Simplified installation via pip and conda environment
* Vacuum to air and heliocentric or barycentric correction (the algorithms used here are courtesy of Yuguang Chen at Caltech)
* Ability of using KOA file names or original file names
* Better provenance and traceability of DRP versions and execution steps in the headers
* Versatile sky subtraction modes including using external sky frames and ability of masking regions
* Formal support system via GitHub issues
  
The pipeline is available for use at WMKO. We are in the process of automating the execution and ingestion of the reduced data into KOA.

For older versions, see :doc:`versions`.


Users
=====

Follow the documentation in this section to reduce KCWI data.

.. toctree::
   :maxdepth: 1

   installing
   quickstart
   configuration
   start_up
   sky_subtraction
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

* Don Neill, Matt Matuszewki, Chris Martin and the KCWI Team at Caltech
* Luca Rizzi, Max Brodheim at W. M. Keck Observatory

Maintenance and support of the DRP are provided as part of the Data Services Initiative (PI: John O'Meara), in collaboration with the National Aeronautics and Space Administration.

.. toctree::
   :maxdepth: 1

   updating_docs
   building


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
