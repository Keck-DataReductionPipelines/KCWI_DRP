===============
Nod and shuffle
===============

In a nod-and-shuffle observation, the telescope alternates between the target
and a nearby sky position while the charge on the detector is shuffled between
two panels. The result is one image with an object panel and a sky panel, with
interleaved exposures of the same total duration. Subtracting the sky panel
from the object panel removes the sky accurately, as long as the nodding was
fast enough to follow the changes in the sky.

How the pipeline recognizes it
------------------------------

An ``OBJECT`` frame is treated as nod and shuffle when the nod-and-shuffle
mask is in (header keyword ``BNASNAM`` for the Blue channel or ``RNASNAM`` for
the Red channel is ``Mask``) and the shutter opened more than once
(``NUMOPEN`` greater than 1). No configuration is needed. The panel positions
come from the ``SHUFROWS``, ``NSHFUP`` and ``NSHFDN`` header keywords.

What is different
-----------------

* The sky is the shuffled sky panel, not a sky model. ``kcwi.sky`` is not
  used, and dark and scattered light subtraction are not applied.
* The object and sky panels are written separately (``_obj.fits`` and
  ``_sky.fits``) and carried through the rest of the pipeline as their own
  cubes (``ocube*`` and ``scube*``). Use them to check that features in the
  sky-subtracted cube are not artefacts of the sky subtraction.
* Pixels outside the panel being used are flagged with the value 64 (see
  :ref:`data_products:Flag values`).
* Standard stars are not observed in nod and shuffle, so the inverse
  sensitivity curve must come from a standard star reduced in the usual way,
  with the same configuration (see :doc:`calibration`).

The products are listed in :ref:`data_products:Nod-and-shuffle products`.
