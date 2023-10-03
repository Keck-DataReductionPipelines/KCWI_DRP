===============
Sky Subtraction
===============

The pipeline performs automatic sky subtraction by identifying areas of the sky
that do not contain object flux. If this fails, it is possible to modify the sky
subtraction algorithm in two ways:

#. by using a different frame as the sky frame and

#. by specifying areas of the target or sky frame to exclude

The latter method is accomplished with a mask file (see `Building a mask file`_)
and will exclude the user-selected object-affected pixels from the calculation
of the sky.

To specify the preferred method for sky subtraction, create a file in the data
directory called ``kcwi.sky``. For each KCWI frame for which you want to modify
the sky subtraction algorithm, enter a row with this format:

``raw_object_file.fits raw_sky_frame.fits <mask.fits>``

For example:

``kb230605_00098.fits kb230605_00100.fits``

If you want to specify an external sky frame, only use the first two columns
and do not specify a mask file.

If you want to use a mask on the original file, make sure that the first two
columns contain the *same* file name and add the mask file.

For example:

``kb230605_00098.fits kb230605_00098.fits kb230605_00098_smsk.fits``

Building a mask file
--------------------

To build a mask file based on ``ds9`` regions, use the script
:func:`kcwi_masksky_ds9.py <kcwidrp.scripts.kcwi_masksky_ds9>` contained in the
``scripts`` directory.

This script uses the ``_intf.fits`` files which are generated as part of the
first execution of the pipeline. You should run the pipeline first and verify
the quality of the sky subtraction. If you are not satisfied with the sky
subtraction, use the procedure described here and run the pipeline again.

To start, display the ``_intf.fits`` frame in ``ds9`` and create regions around
areas of the frame that you would like to exclude. Save the regions to a file
and run the ``kcwi_masksky_ds9`` command indicating the file name
(``_intf.fits``) and the region file.  This will produce a mask file with a
filename that is the same as the input ``_intf.fits`` file, but replacing the
``intf`` with ``smsk``.  That filename can be added to the ``kcwi.sky`` file
as in the example above.

Automatic masking
-----------------

To perform automatic masking, use a different form for the entry in
``kcwi.sky``.  For a bright continuum source that facilitates automatic
object flux rejection add the entry:

``raw_object_file.fits cont``

where the word ``cont`` indicates a continuum source.  For example:

``kr230605_00120.fits cont``

For a fainter continuum source where automatic masking is not possible, add
the entry:

``raw_object_file.fits cont object_position width_of_object``

For example:

``kr230605_00120.fits cont 14.3 3.5``

In the first case, the DRP will find the location and width of the continuum
source automatically.  This case should be used for bright continuum sources.

In the second case, the third and fourth columns specify the location and
width of a fainter continuum source.  In both cases, the pixel positions
across the slice where the continuum source resides are all masked.  The sky is
determined from a small region just outside the continuum source mask.

To determine the pixel position of a faint continuum source, you can load the
``_intf.fits`` frame in ds9 and then load the corresponding position map image
``_posmap.fits`` into an adjacent frame.  If you blink between the two frames
you can determine the ``object_position`` and ``width_of_object`` values
to enter in the ``kcwi.sky`` file.

All standard star observations are automatically masked as bright continuum
sources.

Re-run the reduction of this file
---------------------------------

Rather than re-run the entire pipeline, it is possible to run it only on the
file that needs an improved sky subtraction.

The instructions are listed in :doc:`running`. In summary, you will need to 
specify just one file, the one for which you want to improve the sky
subtraction, and make sure that ``clobber`` is set to True in the config file
and be sure to remove the entry for this file in the proc file (`kcwib.proc` or
`kcwir.proc`).  For example:

``reduce_kcwi -b -f kb230605_00098.fits``

or

``reduce_kcwi -r -f kr230605_00120.fits``

depending on which channel was used for the observation.

Disabling sky subtraction
-------------------------

The DRP makes a general attempt at modeling the sky, however, no sky subtraction
algorithm will be perfect for all objects.  To acknowledge that fact and allow
the user to attempt their own sky subtraction, the DRP creates a new FITS
extension called ``NOSKYSUB`` that contains the data frame without sky
subtraction.

In addition, there are two methods to allow the user to skip sky subtraction.

#. If you enter the word ``skip`` on the line after the raw file in the
   ``kcwi.sky`` file, sky subtraction will be skipped for that file. For example:
   ``kr230605_00120.fits skip``.

#. If you use the command line parameter ``-k``, sky subtraction will be skipped
   for all images processed from that command line.
