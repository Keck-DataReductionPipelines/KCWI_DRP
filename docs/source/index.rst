Keck Cosmic Web Imager Data Reduction Pipeline
==============================================

KCWI_DRP is the official data reduction pipeline (DRP) for the Keck Cosmic Web Imager.
This DRP has been developed by the KCWI team at Caltech in
collaboration with the W. M. Keck Observatory Scientific Software Group.

Release 1.3
===========

The current release is version 1.3.3. Like previous versions, it supports both
the Blue and Red channels of the instrument. This pipeline is based on the
original KCWI IDL pipeline developed by the KCWI team at Caltech, and it is
the only pipeline supported by the W. M. Keck Observatory.

New in 1.3
----------

* Installation is defined by Poetry and ``pyproject.toml`` (``setup.py`` is gone)
* More spectrophotometric standards, including CALSPEC standards; HD 19445,
  HD 60753, HZ 29 and HZ 43 were removed
* Fixed the Tukey taper being applied repeatedly to the continuum bars during
  the dispersion calculation
* The ``NOSKYSUB`` extension is written as 32-bit floats and carries the WCS
* Options for troubleshooting Bokeh and Firefox (see :doc:`issues`)
* Fixed stacked flats (1.3.3): from 1.2.0 to 1.3.2, each stacked flat was made
  from the first input flat only. If flat-field noise matters for your data,
  consider re-reducing data reduced with those versions

For the details of each 1.3.x release, and for older versions, see
:doc:`versions`.

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

For more information on the internals of the pipeline, see the following sections:

.. toctree::
   :maxdepth: 1

   pipeline_concepts
   api

Citing the KCWI DRP
===================

If you use the KCWI DRP in published work, please cite the KCWI instrument
paper, `Morrissey et al. 2018, ApJ, 864, 93
<https://ui.adsabs.harvard.edu/abs/2018ApJ...864...93M>`_, and give the
version of the DRP you used. Every output file records it in a ``HISTORY``
card of the form ``kcwidrp version=1.3.3``.

Other Helpful Packages
======================

The following packages might be helpful to users working with reduced KCWI data.
Please note that these packages are not maintained or developed by W. M. Keck Observatory
or the KCWI team at Caltech. Rather, these are community-developed packages that
are listed only for convenience. **If you utilize any of these packages, please
ensure that you properly cite their developers in addition to the KCWI DRP's.**

KCWIKit
-------
`KCWIKit <https://github.com/yuguangchen1/KcwiKit>`_ is a community developed
package to extend the official KCWI DRP with a
variety of stacking tools and DRP improvements. The software offers masking and
median filtering scripts to be used while running the KCWI DRP, and a step-by-step
KCWI_DRP implementation for finer control over the reduction process. Once the DRP
has finished, KCWIKit can be used to stack the output cubes via the Montage package.
Various functions cross-correlate and mosaic the constituent cubes and the final
stacked cubes are WCS corrected. Helper functions can then be used to deproject
the stacked cube into lower-dimensional representations should the user desire.

Developers
==========

The KCWI DRP was developed in collaboration between:

* Don Neill, Matt Matuszewski, Nik Prusinski, Chris Martin and the KCWI Team at Caltech
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
