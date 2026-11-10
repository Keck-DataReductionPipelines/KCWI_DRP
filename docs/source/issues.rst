==================
FAQ / Known Issues
==================

This page contains a list of Frequently Asked Questions, and known issues related
the operations of the pipeline that are not likely to be fixed in the short term.
Only issues which will not be fixed with a new release of the pipeline within 
around a month will be added to this list. Wherever possible, a link to relevant
GitHub issues will be provided.

Frequently Asked Questions
==========================

Why are my output cubes reduced with this pipeline slightly offset from cubes made with the IDL pipeline?
---------------------------------------------------------------------------------------------------------

This version of the pipeline includes functionality that the IDL pipeline did not,
notably an air-to-vacuum wavelength correction, and a heliocentric (or
barycentric) wavelength correction. These, predictably, result in small but
noticeable offsets between new and old cubes. Both can be changed in the
:ref:`configuration <configuration:Wavelength correction parameters>`.

What are all of these output files?
-----------------------------------

Please reference the :doc:`data_products` page for information on what files are
produced by the pipeline.

I think I might be missing some calibration files, but I'm not sure. How do I check?
------------------------------------------------------------------------------------

There is a script installed with the pipeline called ``check_cals`` that can
help you here. It scans a data directory and pipeline config, and outputs a 
report about what calibrations it matched with what science frames. If any are
missing, it tells you which, and how many are needed. More details can be found
:ref:`here <scripts:check_cals>`.

Can I reduce files with KOA names (``KB.20230605.12345.67.fits``)?
-------------------------------------------------------------------

Yes. Since 1.2.0, input files can have any name; for example
``reduce_kcwi -b -f KB.*.fits -g``. The outputs are named after the input
files (``KB.20230605.12345.67_icubes.fits``), and the rows in ``kcwi.sky``
must use the same names as the files you give the pipeline.

Why was my frame not reduced?
-----------------------------

Look for these in the log (``logs`` directory):

* ``Already processed``: the frame is already in the proc file. See
  :ref:`Re-reducing frames <running:Re-reducing frames>`.
* ``BIAS frame with exposure time ... > 0. Discarding``: a bias with a
  non-zero exposure time (for example a focus image) is ignored.
* A missing calibration, for example no master arc or flat for the frame's
  configuration. A step whose calibration is missing is skipped (the log says
  so in the ``Checking precondition for ...`` messages), and the frame goes on
  without it, so later steps that need its output also produce nothing. Run
  ``check_cals`` to see which calibrations are missing.

Why is there no ``icubes`` file for my science frame?
-----------------------------------------------------

``icubes`` is only made when a standard star taken with the same instrument
configuration was reduced in the same directory, before the science frame. See
:ref:`calibration:How flux calibration works`.

Why did my ``kcwi.sky`` file have no effect?
--------------------------------------------

In 1.3.x, ``kcwi.sky`` is only read when you run ``reduce_kcwi`` in the
directory that contains it and the raw files, and give the raw files without a
directory path. See the warning in :doc:`sky_subtraction`.

Can I reduce Blue and Red data at the same time?
------------------------------------------------

Yes, since 1.2.0. Run one ``reduce_kcwi -b`` and one ``reduce_kcwi -r`` in
separate terminals, in the same directory: each channel has its own proc file
(``kcwib.proc`` and ``kcwir.proc``).

Which ``plot_level`` should I use?
----------------------------------

``1`` (the default) shows and saves the main diagnostic plots without stopping.
``0`` turns plots off, which is fastest for unattended reductions. ``2`` and
``3`` show more plots and stop for you to inspect and answer questions,
mainly during wavelength calibration; use them when a solution needs checking.
See :ref:`configuration:Plotting parameters`.

Known Issues
============

Plot viewer issues
------------------

The live plot viewer (see :ref:`running:Viewing the plots`) is optional:
when it cannot start, the pipeline logs a warning such as ``Plot viewer not
started, live plots disabled; plots are still saved`` and carries on. The PNG
files in ``plots`` are written either way.

* **No address printed.** Check that ``enable_plot_viewer = True`` in your
  configuration file. An old configuration file that only has
  ``enable_bokeh`` still works, with a deprecation warning in the log; replace
  the key with ``enable_plot_viewer``.
* **The page stays empty.** At ``plot_level = 0`` no plots are drawn, so the
  viewer has nothing to show; use 1 or more.
* **All ports taken.** The viewer tries ``plot_viewer_port`` (default 5006)
  and the next 9 ports. If another program uses that whole range, set
  ``plot_viewer_port`` to a free port, for example ``plot_viewer_port = 8100``.
  A leftover ``bokeh serve`` process from version 1.3.x can hold port 5006;
  the viewer skips it, but you can stop it with ``pkill -f "bokeh serve"``.
* **The page does not load from another machine.** The viewer only accepts
  connections from the pipeline machine itself. Use a remote desktop or SSH
  port forwarding, as described in :ref:`running:Viewing the plots`.

Version 3 needs no Bokeh server, Firefox, geckodriver, selenium or Qt
binding. The 1.3.x troubleshooting steps for those (``pkill bokeh``,
``plot_firefox_compat``, installing ``pyqt5``) no longer apply.

Massive Slowdown When Calculating Central Dispersion
----------------------------------------------------

This issue does not throw an error, but can be identified by the logs as it
happens. The logs will look something like ::

    2021-06-08 18:49:51:KCWI:INFO: Using TAPERFRAC = 0.200
    Bar#:   4, Cdisp: 0.2392
    Bar#:   0, Cdisp: 0.2391
    Bar#:   8, Cdisp: 0.2393
    Bar#:  12, Cdisp: 0.2393
    ...
    Bar#: 119, Cdisp: 0.2397

This step typically takes anywhere from 30 seconds to several minutes, depending
on the resources available to your computer. However, sometimes this step takes
upwards of 20 minutes, even on a powerful machine. This appears to be caused by
a conflict in thread allocation between various packages used by the pipeline,
although the specifics remain unknown. 

To fix the issue, you need to specify how threads are allocated directly. This
can be done directly from the command line by typing the following lines into
your terminal:

.. code-block:: bash

    export MKL_NUM_THREADS=16
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1

This will not persist between terminal sessions, so you should add it to your
:code:`.bashrc` file.
