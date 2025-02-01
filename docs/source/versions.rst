.. _versions:

=================
Previous versions
=================



Version 1.1.0 (2023)
--------------------
Added support for red-side data reduction.

Version 1.0 (2021)
------------------

* Simplified installation via pip and conda environment
* Vacuum to air and heliocentric or barycentric correction (the algorithms used here are courtesy of Yuguang Chen at Caltech)
* Ability of using KOA file names or original file names
* Better provenance and traceability of DRP versions and execution steps in the headers
* Versatile sky subtraction modes including using external sky frames and ability of masking regions


Version 0.1 (2019)
------------------

* First end-to-end version including all the reduction step of the IDL pipeline
* Three execution modes

   * Reduce all files in a directory in the order in which they appear
   * Reduce all files after grouping them by file type and in the correct order
   * Monitor a directory for new files and reduce them as they appear

* Multi-threading for CPU intensive tasks such as wavelength calibration
* Multi-processing for large datasets

