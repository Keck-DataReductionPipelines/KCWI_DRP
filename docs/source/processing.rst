===================
Processing overview
===================

This page describes what the pipeline does with each kind of frame, and in
what order. For the internals (events, recipes and primitives), see
:doc:`pipeline_concepts`.

How frames are routed
---------------------

Each raw frame is classified by its ``IMTYPE`` header keyword and sent to the
recipe for that type. A frame that is already in the proc file is skipped
unless ``clobber = True`` (see
:ref:`Re-reducing frames <running:Re-reducing frames>`). Bias frames with an
exposure time above zero are discarded.

Master calibrations are made once enough frames of a kind have been processed
with the same configuration. The minimum numbers are set in the configuration
file (``bias_min_nframes``, ``flat_min_nframes`` and so on; see
:doc:`configuration`) and match the afternoon calibration scripts used at
WMKO. The ``check_cals`` script (see :doc:`scripts`) tells you before a
reduction whether each science frame has the calibrations it needs.

.. list-table::
    :header-rows: 1
    :widths: 16, 54, 30

    * - ``IMTYPE``
      - Steps
      - Main outputs
    * - ``BIAS``
      - Overscan subtraction and trimming; combined into a master bias
      - ``intb``, ``mbias_<CCDCFG>``
    * - ``DARK``
      - Basic reduction; combined into a master dark
      - ``int``, ``mdark_<STATEID>``
    * - ``CONTBARS``
      - Basic reduction; combined into a master; the continuum bars are
        found and traced to solve the slice geometry
      - ``mcbars``
    * - ``ARCLAMP``
      - Basic reduction; combined into a master; arc spectra extracted along
        the bars and matched to the arc lamp atlas (``ThAr`` or ``FeAr``) to
        solve the wavelength scale; geometry and wavelength solution combined
        into the maps used to build cubes
      - ``marc``, ``_geom.pkl``, ``wavemap``, ``slicemap``, ``posmap``,
        ``delmap``
    * - ``FLATLAMP``, ``DOMEFLAT``, ``TWIFLAT``
      - Basic reduction, dark and scattered light subtraction; stacked, then
        made into a master flat (illumination and response correction)
      - ``sflat``/``sdome``/``stwif``, ``mflat``/``mdome``/``mtwif``
    * - ``OBJECT``
      - Basic reduction, dark and scattered light subtraction, flat
        correction, sky subtraction, cube building, wavelength corrections,
        DAR correction; standard stars make an inverse sensitivity curve;
        all objects are flux calibrated if a curve is available
      - ``intf``, ``intk``, ``icube``, ``icubew``, ``icubed``, ``icubes``
        (and ``invsens`` for standards)
    * - ``OBJECT`` (nod and shuffle)
      - As ``OBJECT``, but the sky is the shuffled sky panel; see
        :doc:`nod_and_shuffle`
      - ``intk``, ``icube*``, ``ocube*``, ``scube*``

"Basic reduction" is: overscan subtraction and trimming, bias subtraction,
conversion to electrons, defect correction, cosmic ray rejection, creation of
the uncertainty image, and rotation of the image to a common orientation. Not
every frame type gets every step; :doc:`api` lists the primitives.

The products are described in :doc:`data_products`.

Group mode order
----------------

With ``-g`` (group mode), the pipeline first reads all the frames, then
processes them by type in this order, so that every master calibration exists
before the frames that need it:

#. ``BIAS``
#. ``CONTBARS``
#. ``ARCLAMP``
#. ``FLATLAMP``
#. ``DOMEFLAT``
#. ``TWIFLAT``
#. ``OBJECT``, with standard stars first, so that their inverse sensitivity
   curves are available to flux calibrate the other objects

Focus frames (``OBJECT = focus``) and CCD clearing frames are left out.

.. note::

    In 1.3.x, group mode does not process ``DARK`` frames. If you need a master
    dark, reduce the darks first with a file list (``-l``) or ``-f``.

Without ``-g``, frames are processed in the order they are given, so the
calibrations must come before the frames that use them.

Choosing the calibrations for a frame
-------------------------------------

For each frame, the pipeline looks in the proc file for master calibrations
from the same channel that match the frame:

* master biases and darks must have the same CCD configuration (``CCDCFG``:
  binning, readout mode, gain setting and amplifier mode); a master dark is
  scaled to the frame's exposure time;
* the other calibrations (bars, arcs, flats, inverse sensitivity) must have the
  same instrument state (``STATEID``, which is set for each instrument
  configuration).

When several match, the one nearest in time is used. This is why the
calibrations for a frame must be reduced in the same directory, with the same
proc file, as the frame itself.
