===========
Calibration
===========

In order to flux calibrate your data, an appropriate spectro-photometric
standard star must be chosen and observed.  The standard stars that are
compatible with the pipeline can be found in the kcwidrp/data/stds
directory in the pipeline install directory (KCWI_DRP).  In that directory, you
will find a WMKO starlist file called
`kcwi_stds_starlist.txt <https://github.com/Keck-DataReductionPipelines/KCWI_DRP/blob/master/kcwidrp/data/stds/kcwi_stds_starlist.txt>`_
with the stars in RA order.

You will also find in that directory plots of each standard star showing the
wavelength sampling and coverage for each standard star reference file.  Once
you have defined your configuration, it is recommended that you examine the
star list to narrow down the list of available standards during your observing
run, and then examine the plots of your narrowed list to be sure that any
standard you select is well sampled and covers the wavelength range of your
configuration.

You can also examine the standard star plots in the github repository by
following
`this link <https://github.com/Keck-DataReductionPipelines/KCWI_DRP/tree/master/kcwidrp/data/stds>`_.
Click on the \*.png files and they should display in your browser with flux
versus wavelength and the reference samples indicated by plus (+) signs.

Be aware that high resolution configurations usually cover smaller wavelength
ranges and in that case, you will need to find a densely sampled standard star.
You will need at least five points within your wavelength range, although this
will only provide a crude calibration (which may be appropriate for your needs).

If you copy the entries for the standards you want to observe from the starlist
here into your own starlist, then the observations will have the correct object
names, making sure the pipeline will recognize them.  This is highly
recommended.


How the pipeline recognizes a standard star
-------------------------------------------

An ``OBJECT`` frame is treated as a standard star when its ``OBJECT`` header
keyword (or, failing that, its ``TARGNAME`` keyword), in lower case and with
spaces removed, matches the name of a reference spectrum in
``kcwidrp/data/stds``. For example, ``OBJECT = 'Feige 34'`` matches
``feige34.fits``. The names in ``kcwi_stds_starlist.txt`` are the names the
pipeline expects, which is why copying the entries from that starlist is
recommended.

If a standard star was observed under another name, change its ``OBJECT``
keyword before reducing it, for example with the ``modhead`` script (see
:doc:`scripts`):

.. code-block:: bash

    modhead kb230605_00120.fits -k OBJECT -v feige34

How flux calibration works
--------------------------

#. The standard star is reduced like any other object, up to the DAR-corrected
   cube (``icubed``). The sky is masked around the star automatically (see
   :doc:`sky_subtraction`).
#. ``MakeInvsens`` extracts the star's spectrum, compares it with the reference
   spectrum and fits an inverse sensitivity curve. It writes
   ``_invsens.fits`` and ``_ea.fits`` (effective area), and plots of the fit,
   its residuals, and the effective area and efficiency. Check the residual
   plot in the ``plots`` directory.
#. ``FluxCalibrate`` multiplies each object cube by the inverse sensitivity
   curve and writes the flux-calibrated cube, ``icubes``.

For a science frame, ``FluxCalibrate`` uses the inverse sensitivity curve, from
the proc file, that has the same channel and instrument state (``STATEID``) as
the frame and is nearest in time. So:

* The standard star must be reduced in the same directory (with the same proc
  file) as the science frames, and before them. Group mode (``-g``) reduces
  standard stars before other objects for this reason (see :doc:`processing`).
* A standard star is only used for science frames taken with the same
  instrument state. If no curve matches, no ``icubes`` file is written for
  that frame; the ``icubed`` cube is still available, in electrons.
