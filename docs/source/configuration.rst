========================
Configuration Parameters
========================

A number of reduction parameters can be changed using entries in the
configuration file.

If you installed the pipeline with ``pip``, the configuration file will not be
easy to find, since it will be stored with installed ``pip`` packages. We recommend
editing a copy of the config instead. You can create a copy of the config file
by invoking the pipeline with the ``--write_config`` option:

    .. code-block:: bash

        reduce_kcwi --write_config

This will create an editable file wherever you run the command from. You can
make whatever changes you need to this file, and then run the pipeline with the
``-c`` option (see :ref:`main reduction script <running:Other command line options>`)
pointing to that file. For example,

    .. code-block:: bash

        reduce_kcwi -b -c PATH/TO/CONFIG/kcwi.cfg -f kb0001.fits

The configuration file contains a number of parameters connected to the
structure of the files and to the specifications of the instrument, e.g., the
number of continuum bars. There is usually no reason to modify these parameters,
and they are not described here.

.. note::

    If a parameter in the ``kcwi.cfg`` file is not described here, then you can
    assume that it should not be modified.

The remaining parameters are used to control the processing algorithms and are
described here.

Blue and Red sections of the configuration file
-----------------------------------------------

Most parameters are in the ``[KCWI]`` section and apply to both channels. A
few have different defaults for each channel, and are set again in the
``[BLUE]`` and ``[RED]`` sections. For example, to deal with cosmic rays, the
Red channel uses a median stack of three continuum bars images and a median
stack of three arc lamp images, while the Blue channel only requires one each
of those images.

When you run with ``-b`` or ``-r``, the values in the ``[BLUE]`` or ``[RED]``
section replace those in ``[KCWI]`` for these parameters only:

``arc_min_nframes``, ``contbars_min_nframes``, ``object_min_nframes``,
``minoscanpix``, ``oscanbuf``, ``procfile``, ``LINETHRESH``, ``TUKEYALPHA``,
``FRACMAX``

.. warning::

    To change one of these parameters, edit it in the ``[BLUE]`` or ``[RED]``
    section: a value in ``[KCWI]`` is overwritten by the channel section.
    ``LINETHRESH``, ``TUKEYALPHA`` and ``FRACMAX`` can also be set on the
    command line (``-e``, ``-u`` and ``-F``), which takes precedence over both.
    Adding any other parameter to ``[BLUE]`` or ``[RED]`` has no effect.

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
scripts that are used at WMKO for afternoon calibrations.

The overscan parameters are based on the configuration of the Blue and Red
detectors and we do not recommend altering these parameters, but show them
here to illustrate that they differ between channels.

.. code-block:: python

 skipscat = False        # Skip subtracting scattered light?
 skipsky = False         # Skip sky subtraction?

These control scattered light subtraction and sky subtraction.  Setting either
of these to ``True`` will skip the subtraction for all subsequent runs of the
pipeline.  Skipping sky subtraction globally can also be invoked by using
``-k`` on the :ref:`command line <running:Other command line options>`.

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
 LINETHRESH = 100.      # Line threshold for fitting

These adjust the way in which arc line fitting is performed.  Most of these
parameters are also available on the
:ref:`command line <running:Other command line options>`.

In most cases, you will not have to adjust these.  For the Red channel, we use
these default values:

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

Output and re-reduction parameters
----------------------------------

.. code-block:: python

 clobber = False
 output_directory = "redux"
 procfile = kcwib.proc      # = kcwir.proc for [RED]

``clobber``: when ``False`` (the default), frames that are already in the proc
file are skipped. Set it to ``True`` to re-reduce frames; see
:ref:`Re-reducing frames <running:Re-reducing frames>`.

``output_directory``: where the products are written. Do not change it: in
1.3.x several steps read their inputs from ``redux`` whatever this is set to.

``procfile``: the name of the processing table, set per channel. It can also be
set with ``-p`` on the command line.

Wavelength correction parameters
--------------------------------

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

 enable_plot_viewer = True
 plot_viewer_port = 5006
 plot_viewer_open_browser = False
 plot_level = 1

Plots are drawn with Matplotlib. Two things are controlled separately:

* **Saved plots.** At ``plot_level`` 1 or more, every plot is saved as a PNG
  file in the ``plots`` directory, whatever the viewer settings are.
* **The live plot viewer.** With ``enable_plot_viewer = True``, each pipeline
  process serves its plots as a web page while it runs. At start-up it prints
  the address to open, for example::

     ========================================================================
       Plot viewer: http://127.0.0.1:5006/  (open this in a browser)
     ========================================================================

  The page shows each plot as it is drawn, with pan, zoom, a cursor readout
  and a download button. See :ref:`running:Viewing the plots`.

``plot_viewer_port`` is the first port tried. If it is taken (for example by
a second reduction running at the same time), the next 9 ports are tried, and
the address actually used is printed. Each process has its own viewer and its
own page.

``plot_viewer_open_browser = True`` opens the page in a browser tab
automatically; by default, open the printed address yourself.

If the viewer cannot start, the pipeline logs a warning and carries on: plots
are still saved, only the live page is missing.

The ``plot_level`` parameter controls the level of interactivity. Setting it 0
will disable plotting.  Setting it to 1 will enable plotting, but the DRP will
not interact with the user. A higher level will increase both the verbosity and
the interactivity of the plots. The highest level is 3. At this level, the user
will be provided with a plot of every arc line, for example, with a graphic
representation of the fitting used to determine the central position.

For general use, it is advisable to leave the plot level to 1.

The size of the plotting window can be specified using ``plot_width`` and
``plot_height``.

.. note::

   Configuration files from versions before 3 may contain ``enable_bokeh``,
   ``terminate_on_failed_bokeh_start``, ``plot_firefox_compat`` or
   ``plot_prewarm_firefox``. They still load without error. The last three are
   ignored. ``enable_bokeh`` is deprecated: it is used only when
   ``enable_plot_viewer`` is absent, and a warning is logged. Replace it with
   ``enable_plot_viewer``.

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

These parameters are used to control the CRR algorithms. See the documentation
in `astroscrappy <https://astroscrappy.readthedocs.io/en/latest/index.html>`_
for details.

Combining science images
------------------------

The following parameter controls what is done with science frames, including
standard star observations.  There is one in each of the ``[BLUE]`` and ``[RED]``
sections of the parameter file, so you can specify a non-default value for one
channel at a time.

.. code-block:: python

 object_min_nframes = 1

If you wish to combine images (to mitigate CRs, e.g.), you can set this to the
number of identically positioned science images that you want to combine.  There
is an additional step required.  The FITS header keyword GROUPID for the images
to be combined must be identical.  This can be achieved with any FITS header
keyword editor.  It is recommended that you combine at least three images to
have a good chance of mitigating CRs.  In the case that three images are
combined, the routine :class:`kcwidrp.primitives.MakeMasterObject` will use the
``median`` method for combining the images.  If there are more than three, it
will use the ``average`` method.  In both cases, it does a 2-sigma high clip.

A master object file will be output using the base name of the first file in
the group.  For example:

``kr230605_00119_mobj.fits``

This file will then be processed through the rest of the pipeline and result in
a combined, calibrated data cube:

``kr230605_00119_icubes.fits``

which should have fewer CRs than an individual frame.  Be aware, however, that if
you have only three images and use the ``median`` combine method, the S/N ratio
will be equivalent to a single exposure.  The other frames in the group will not
be processed beyond the point where they are combined (currently just after
illumination correction: ``_intf.fits`` creation).

Instrument and observatory constants
------------------------------------

These describe the instrument, the detector and the observatory. Do not change
them.

.. code-block:: python

 NBARS = 120              # number of continuum bars
 REFBAR = 57              # reference bar for the geometry solution
 PIXSCALE = 0.00004048    # unbinned pixel size in degrees
 SLICESCALE = 0.00037718  # Large slicer slice size in degrees
 PIX = 0.0150             # pixel size in mm
 FCAM = 305.              # focal length of the camera (mm)
 GAMMA = 4.               # mean out-of-plane angle for diffraction (deg)
 WAVEFID = 3000.          # fiducial wavelength for wavelength bins (Angstrom)
 ROTOFF = 0.              # rotator/IFU offset angle in degrees
 KNOTSPP = 1.25           # B-spline knots per image pixel for the sky model
 psfwid = 30              # nominal window for a point source (unbinned px)
 latitude = 19.82656      # Keck II, used for the radial velocity correction
 longitude = -155.4742
 altitude = 4160
