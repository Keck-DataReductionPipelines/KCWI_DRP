.. Keck Cosmic Web Imager DRP documentation master file, created by
   sphinx-quickstart on Fri May 29 10:09:02 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Keck Cosmic Web Imager DRP's
============================

KCWI_DRP is a Python based data reduction pipeline (DRP) developed by the KCWI team at Caltech in
collaboration with the W. M. Keck Observatory Scientific Software Group.
While the algorithms and the reduction steps are based on the previous IDL pipeline (KDERP),
this pipeline uses a new frame-by-frame, event-driven approach, and relies on the Keck DRP Framework (provide link)

Release 0.1
===========

What this version provides
--------------------------

* First end-to-end version including all the reduction step of the IDL pipeline
* Three execution modes

   * Reduce all files in a directory in the order in which they appear
   * Reduce all files after grouping them by file type and in the correct order
   * Monitor a directory for new files and reduce them as they appear

* Multi-threading for CPU intensive tasks such as wavelength calibration
* Multi-processing for large datasets

What this version is missing
----------------------------

What are we working on?

Users
=====

Follow the documentation in this section to reduce KCWI data.

.. toctree::
   :maxdepth: 1

   installing
   quickstart
   configuration
   start_up


More information
================

.. toctree::
   :maxdepth: 1

   pipeline_concepts
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
