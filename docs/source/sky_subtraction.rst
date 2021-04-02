===============
Sky Subtraction
===============

The pipeline performs automatic sky subtraction by identifying areas of the sky that do not contain objects.
If this fails, it is possible to modify the sky subtraction algorithm in two ways: (1) by using a different frame as the sky frame
and (2) by specifying areas of the target or sky frame to exclude from the calculation of the sky (because they are 
contaminated by an object, for example).

To specify the preferred method for sky subtraction, create a file in the data directory called ``kcwi.sky``. For each KCWI frame 
for which you want to modify the sky subtraction algorithm, enter a row with this format:

``raw_object_file.fits raw_sky_frame.fits <mask.fits>``

If you want to specify an external sky frame, only use the first two columns and do not specify a mask file.

If you want to use a mask on the original file, make sure that the first two columns contain the *same* file name and 
add the mask file.

Building a mask file
--------------------

To build a mask file based on ``ds9`` regions, use the script ``kcwi_masksky_ds9.py`` contained in the ``scripts`` directory.

This script uses the ``_intf.fits`` files which are generated as part of the first execution of the pipeline. You should run the pipeline first and verify the quality  
of the sky subtraction. If you are not satisfied with the sky subtraction, use the procedure described here and run the pipeline again.

To start, display the ``_intf.fits`` frame on ``ds9`` and create regions around areas of the frame that you would like to exclude. 
Save the regions to a file and run the kcwi_mask_sky command indicating the file name (``_intf.fits``) and the region file. 
This will produce a mask file that can be added to the ``kcwi.sky`` file. 


