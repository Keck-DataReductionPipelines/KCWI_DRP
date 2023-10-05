===========
Calibration
===========

In order to flux calibrate your data, an appropriate spectro-photometric
standard star must be chosen and observed.  The standard stars that have been
incorporated into the pipeline can be found in the kcwidrp/data/stds
directory.  In that directory, you will find a file called
`kcwi_stds_starlist.txt <https://github.com/Keck-DataReductionPipelines/KCWI_DRP/tree/master/kcwidrp/data/stds/kcwi_stds_starlilst.txt>`_ in the WMKO starlist format.

You will also find in that directory plots of each standard star showing the
wavelength sampling and coverage for each standard star reference file.  Once
you have defined your configuration, it is recommended that you examine the
star list to narrow down the list of available standards during your observing
run, and then examine the plots to be sure that any standard you select
is well sampled and covers the wavelength range of your configuration.

You can also examine them from this web page by following `this link <https://github.com/Keck-DataReductionPipelines/KCWI_DRP/tree/master/kcwidrp/data/stds>`_.
Click on the \*.png files and they should display in your browser with flux
versus wavelength and the reference samples indicated by plus (+) signs.

Be aware that high resolution configurations usually cover smaller wavelength
ranges and in that case, you will need to find a densely sampled standard star.
You will need at least five points within your wavelength range, although this
will only provide a crude calibration (which may be appropriate for your needs).

