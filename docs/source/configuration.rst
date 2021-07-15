========================
Configuration Parameters
========================

A number of reduction parameters can be changed using entries in the configuration file.

If you installed the pipeline with ``pip``, the configuration file will not be
easy to find, since it will be stored with installed ``pip`` packages. If you 
need to modify the configuration, we advise that you download a copy of the 
`config file <https://github.com/Keck-DataReductionPipelines/KCWI_DRP/blob/master/kcwidrp/configs/kcwi.cfg>`_,
saving it to some easy to remember location, and using the ``-c`` option to point
to that file.

If you installed the package via ``git``, the master copy of the configuration
file is in the installation directory, in ``kcwidrp/config/kcwi.cfg``.  This
file can be modified in place or copied in another directory (the ``-c`` option
of the main reduction script is used to specify the configuration file).

The configuration file contains a number of parameters connected to the structure of the files
and to the specifications of the instrument, e.g. the number of continuum bars.
There is usually no reason to modify these parameters, and they are not described here.

The remaining parameters are used to control the processing algorithms and are described here.

Processing parameters
---------------------

.. code-block:: python

 output_directory = "redux"

This parameter specifies the output directory for the data products. If this directory
is missing, it will be created automatically.

.. code-block:: python

 bias_min_nframes = 7
 flat_min_nframes = 6
 dome_min_nframes = 3
 twiflat_min_nframes = 1
 dark_min_nframes = 3

These parameters control the minimum number of bias, internal/dome/twilight flats and darks
that the DRP expects before producing a master calibration.
The values shown here are synchronized with the calibration scripts that we use in the
afternoon.

.. code-block:: python

 clobber = False

This parameters controls the behaviour of the DRP if one of the data product has already
been generated: set ``clobber = True`` to overwrite existing products.

.. code-block:: python

 skipscat = False    # Skip subtracting scattered light?

This parameter disables the subtraction of scatteres light if set to ``True``. In some case
the subtraction of scattered light can produce unexpected results.

.. code-block:: python

 default_arc_lamp = 'ThAr'

KCWI has two calibration lamps, Thorium/Argon (ThAr) and Iron/Argon (FeAr). This parameter
specifies which of the two lamps should be used by the DRP. The default is to use the ThAr lamp.


Plotting parameters
-------------------

.. code-block:: python

 plot_pause = 1
 saveintims = False
 inter = 1
 plot_width=1000
 plot_height=600

 # BOKEH SERVER
 enable_bokeh = False
 plot_level = 0

These parameters control the plotting features of the DRP. Plotting is performed using
a combination of a Bokeh server running in the background and a browser front end.

To activate the plotting features, set ``enable_bokeh = True``. When the DRP starts, it will
check if there is an instance of the Bokeh server running or start one. A browser
window will be opened automatically when needed.

The ``plot_level`` parameter controls the level of interactivity. Setting it 0 will disable
interactive fetures: the DRP will produce plots when needed but it will not interact
with the user. A higher level will increase both the verbosity and the interactivity of the
plots. The highest level is 3 (CHECK). At this level, the user will be provided with a plot
of every arc line, for example, with a graphic representation of the fitting used to determine
the central position.

For general use, it is advisable to leave the plot level to 1.

The ``plot_pause`` parameter controls how long the DRP will pause between automatically generated
plots (in seconds).
Finally, the ``saventims`` parameter controls the generation of JPG diagnostics plots saved
in the current directory.

The size of the plotting window can be specified using ``plot_width`` and ``plot_height``.

Cosmic rays rejection parameters
--------------------------------

.. code-block:: python

 CRR_MINEXPTIME = 60.0
 CRR_PSSL = 0.0
 CRR_GAIN = 1.0
 CRR_READNOISE = 3.2
 CRR_SIGCLIP = 4.5
 CRR_SIGFRAC = 0.3
 CRR_OBJLIM = 4.0
 CRR_PSFFWHM = 2.5
 CRR_FSMODE = "median"
 CRR_PSFMODEL = "gauss"
 CRR_SATLEVEL = 60000.0
 CRR_VERBOSE = False
 CRR_SEPMED = False
 CRR_CLEANTYPE = "meanmask"
 CRR_NITER = 4

These parameters are used to control the CRR algorithms. See the documentation in
`astroscrappy <https://astroscrappy.readthedocs.io/en/latest/index.html>`_ for details (PROVIDE LINK)

Wavelength correction parameters:
---------------------------------

The ``radial_velocity_correction`` parameter controls what reference frame to use for
radial velocity corrections. The options are ``heliocentric``, ``barycentric``,
or ``none``

The ``air_to_vacuum`` parameter controls if the pipeline should convert
to vacuum wavelengths from air wavelengths.