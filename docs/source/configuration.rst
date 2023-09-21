========================
Configuration Parameters
========================

A number of reduction parameters can be changed using entries in the
configuration file.

If you installed the pipeline with ``pip``, the configuration file will not be
easy to find, since it will be stored with installed ``pip`` packages. If you 
need to modify the configuration, we advise that you download a copy of the 
`config file <https://github.com/Keck-DataReductionPipelines/KCWI_DRP/blob/master/kcwidrp/configs/kcwi.cfg>`_,
saving it to some easy to remember location, and using the ``-c`` option to
point to that file.

If you installed the package via ``git``, the master copy of the configuration
file is in the installation directory, in ``kcwidrp/config/kcwi.cfg``.  This
file can be modified in place or copied in another directory (the ``-c`` option
of the main reduction script is used to specify the configuration file).

The configuration file contains a number of parameters connected to the
structure of the files and to the specifications of the instrument, e.g., the
number of continuum bars. There is usually no reason to modify these parameters,
and they are not described here.

Note:
    If a parameter in the ``kcwi.cfg`` file is not described here, then you can
    assume that it should not be modified.

The remaining parameters are used to control the processing algorithms and are
described here.

Blue and Red sections of the configuration file
-----------------------------------------------

Now that the Red channel has been installed, there is a need to specify
different default parameters for each channel.  These are delineated in the
config file with ``[BLUE]`` and ``[RED]`` section headers.  For example, to
deal with cosmic rays, the Red channel uses a median stack of three continuum
bars images and a median stack of three arc lamp images, while the Blue channel
only requires one each of those images.

Processing parameters
---------------------

.. code-block:: python

 bias_min_nframes = 7
 flat_min_nframes = 6
 dome_min_nframes = 3
 twiflat_min_nframes = 1
 dark_min_nframes = 3
 arc_min_nframes = 1        # = 3 for [RED]
 contbars_min_nframes = 1   # = 3 for [RED]
 minoscanpix = 75           # = 20 for [RED]
 oscanbuf = 20              # = 5 for [RED]

These parameters control the minimum number of bias, internal/dome/twilight
flats and darks that the DRP expects before producing a master calibration. The
arcs and contbars minimum numbers are different for the Blue and Red channels as
described above.  The values shown here are synchronized with the calibration
scripts that we use in the afternoon.

The overscan parameters are based on the configuration of the Blue and Red
detectors and we do not recommend altering these parameters.

.. code-block:: python

 skipscat = False        # Skip subtracting scattered light?
 skipsky = False         # Skip sky subtraction?

These control scattered light subtraction and sky subtraction.  Setting either
of these to ``True`` will skip the subtraction for all subsequent runs of the
pipeline.  Skipping sky subtraction globally can also be invoked by using
``-k`` on the command line.

.. code-block:: python

 plot_pause = 1         # Pause time for each plot in seconds
 saveintims = False     # Save intermediate images during basic reduction
 verbose = 1            # Verbosity of output

These control various aspects of how the DRP runs: how long to pause at each
plot if not in interactive mode (see ``plot_level`` below), whether to save
intermediate images for diagnosis, and the verbosity level of text output.

.. code-block:: python

 TAPERFRAC = 0.2        # Adjusts edge taper for Atlas cross-correlation
 TUKEYALPHA = 0.2       # Tukey alpha value for cross-correlating bars
 FRACMAX = 0.5          # How much of line peak to use for fitting
 MIDFRAC = -1.0         # Middle fraction or -1 to use default calculation
 ATOFF = 0              # Atlas offset or 0 to use default calculation
 LINELIST = ""          # Optional line list to use instead of generated
 LINETHRESH = 100.      # Line threshhold for fitting

These adjust the way in which arc line fitting is performed.  In most cases, you
will not have to adjust these.  For the Red channel, we use these values:

.. code-block:: python

 TUKEYALPHA = 0.7
 FRACMAX = 0.25
 LINETHRESH = 10.

See the ``[RED]`` section to make changes for Red channel data.

.. code-block:: python

 default_arc_lamp = 'ThAr'

KCWI has two calibration lamps, Thorium/Argon (ThAr) and Iron/Argon (FeAr).
This parameter specifies which of the two lamps should be used by the DRP.
The default is to use the ThAr lamp.

Wavelength correction parameters:
---------------------------------

.. code-block:: python

 radial_velocity_correction = "heliocentric"
 air_to_vacuum = True   # Defaults to vacuum wavelengths

These control the refinement of the wavelength solution.  You can specify if you
want air wavelengths by setting ``air_to_vacuum`` to ``False``.  You can
specify the type of radial velocity correct as one of:

* heliocentric
* barycentric
* none


Plotting parameters
-------------------

.. code-block:: python

 # BOKEH SERVER
 enable_bokeh = False
 plot_level = 1

These parameters control the plotting features of the DRP. Plotting is
performed using a combination of a Bokeh server running in the background and a
browser front end.

To activate the plotting features, set ``enable_bokeh = True``. When the DRP
starts, it will check if there is an instance of the Bokeh server running or
start one. A browser window will be opened automatically when needed.

The ``plot_level`` parameter controls the level of interactivity. Setting it 0
will disable plotting.  Setting it to 1 will enable plotting, but the DRP will
not interact with the user. A higher level will increase both the verbosity and
the interactivity of the plots. The highest level is 3. At this level, the user
will be provided with a plot of every arc line, for example, with a graphic
representation of the fitting used to determine the central position.

For general use, it is advisable to leave the plot level to 1.

The size of the plotting window can be specified using ``plot_width`` and
``plot_height``.

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
