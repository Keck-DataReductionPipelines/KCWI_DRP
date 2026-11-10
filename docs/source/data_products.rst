=============
Data Products
=============

The pipeline writes its outputs to the ``redux`` directory (set by
``output_directory`` in the configuration file). Each output file is named
after the raw frame it was made from, with a suffix: for example
``kb230605_00098_icubes.fits`` is the flux-calibrated cube made from
``kb230605_00098.fits``. Diagnostic plots go to the ``plots`` directory and
logs to the ``logs`` directory.

Which file do I want?
---------------------

For most science work you want the ``icubes`` cube (flux calibrated), or the
``icubed`` cube if no standard star was reduced with your data. Both are
corrected for differential atmospheric refraction and use vacuum wavelengths
corrected to the heliocentric frame, unless you changed the
:ref:`wavelength correction parameters <configuration:Wavelength correction parameters>`.

Each cube has its sky-subtracted data in the primary extension and the same
data *without* sky subtraction in the ``NOSKYSUB`` extension (see
`FITS extensions`_), so you can do your own sky subtraction.

Products
--------

The products below are listed roughly in the order in which they are made.
Products made from science frames (image type ``OBJECT``) are made for
standard stars too.

.. csv-table:: Data products
    :file: outputs/data_products.csv
    :header-rows: 1
    :widths: 12, 18, 18, 52

The flux unit of ``icubes`` is 10\ :sup:`-16` erg/s/cm\ :sup:`2`/Å. In 1.3.x
the ``BUNIT`` keyword written for these files has the wrong sign in the
exponent (a factor of 1e+16 instead of 1e-16); read it as 10\ :sup:`-16`
(`issue #199 <https://github.com/Keck-DataReductionPipelines/KCWI_DRP/issues/199>`_).

Master calibration files
------------------------

Master calibrations are combined from several raw frames. Most are named after
the first raw frame of the set.

.. csv-table:: Master files
    :file: outputs/masters.csv
    :header-rows: 1
    :widths: 25, 20, 55

Nod-and-shuffle products
------------------------

Nod-and-shuffle observations (see :doc:`nod_and_shuffle`) produce these extra
files alongside the ``intk`` and ``icube*`` products:

.. csv-table:: Nod-and-shuffle products
    :file: outputs/nands_products.csv
    :header-rows: 1
    :widths: 22, 18, 60

In 1.3.x the radial velocity and air-to-vacuum corrections are applied to the
``icube`` cube only: ``ocubed``, ``scubed`` and ``dcubed`` are made from cubes
without those corrections.

Intermediate images
-------------------

When ``saveintims = True`` in the configuration file, the pipeline also writes
the image after each basic reduction step. These are for diagnosing problems
and are not needed for anything else.

.. csv-table:: Intermediate images
    :file: outputs/intermediate.csv
    :header-rows: 1
    :widths: 12, 28, 60

FITS extensions
---------------

Images and cubes are written with these extensions:

``PRIMARY``
    The data, in the units given by ``BUNIT``.

``UNCERT``
    The 1-sigma uncertainty of the data, in the same units.

``MASK``
    Pixel mask (0 = good).

``FLAGS``
    Bit flags describing what happened to each pixel; see `Flag values`_.

``NOSKYSUB``
    The data before sky subtraction (sky-subtracted products only, from
    ``intk`` on). It carries the same WCS as the primary extension.

Flag values
-----------

The ``FLAGS`` extension is the sum of the following values. For example, a
value of 6 means the pixel was both a known detector defect (2) and hit by a
cosmic ray (4).

.. list-table::
    :header-rows: 1
    :widths: 10, 90

    * - Value
      - Meaning
    * - 2
      - Detector defect, replaced by the median of the neighbouring good pixels
        (``CorrectDefects``)
    * - 4
      - Cosmic ray, cleaned (``RemoveCosmicRays``)
    * - 8
      - Saturated (``FlagSaturation``)
    * - 64
      - Nod-and-shuffle: outside the object or sky panel used
        (``NandshuffSubtractSky``)
    * - 128
      - Padding added by the DAR correction; the pixel has no data
        (``CorrectDar``, ``icubed`` and later)

In cubes, the flags of the image pixels are resampled onto the cube grid along
with the data, so treat any non-zero flag as suspect rather than relying on the
exact sum. The handling of flags and masks is being reviewed
(`issue #98 <https://github.com/Keck-DataReductionPipelines/KCWI_DRP/issues/98>`_).

Header keywords
---------------

Every product records its processing history in ``HISTORY`` cards: the
``kcwidrp`` version (``kcwidrp version=...``), the git version when the pipeline
was installed from a git checkout, and the name of each primitive that was
applied. Each step also sets a keyword saying whether it was done, for example:

* ``SKYCOR``: sky subtracted?
* ``STDCOR``: flux calibrated?
* ``VCORR``, ``VCORRTYP``: radial velocity correction (km/s) and its type
* ``NSATFLAG``: number of saturated pixels flagged
* ``NBPCLEAN``: number of defect pixels cleaned
* ``BPFILE``: the defect list used

The processing table
--------------------

The pipeline records every product it makes in the proc file (``kcwib.proc``
or ``kcwir.proc``) in the directory where it runs. This is a text table with
one row per frame and product type; the pipeline uses it to find the right
master calibrations for each frame and to skip frames that are already
reduced. See :ref:`Re-reducing frames <running:Re-reducing frames>`.
