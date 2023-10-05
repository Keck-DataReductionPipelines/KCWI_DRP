===========
Calibration
===========

In order to flux calibrate your data, an appropriate spectro-photometric
standard star must be chosen and observed.  The standard stars that are
compatible with the pipeline can be found in the kcwidrp/data/stds
directory in the pipeline install directory (KCWI_DRP).  In that directory, you
will find a WMKO starlist file called
`kcwi_stds_starlist.txt <https://github.com/Keck-DataReductionPipelines/KCWI_DRP/tree/kcrm_merge/kcwidrp/data/stds/kcwi_stds_starlist.txt>`_
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
`this link <https://github.com/Keck-DataReductionPipelines/KCWI_DRP/tree/kcrm_merge/kcwidrp/data/stds>`_.
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

