===========
Quick Start
===========

For users who don't want to know too much about KCWI_DRP (and why not?), here
is a quickstart guide.

The assumption is that you have a directory containing KCWI data, and that the
names of the files are those assigned at the telescope, ``kb*.fits`` and
``kr*.fits``.

Give a quick look at the :doc:`configuration parameters <configuration>` for the
pipeline. To get an editable copy of the default configuration file
(``kcwi.cfg``) in the current directory, run ``reduce_kcwi --write_config``.

For a quick start, it is enough to decide if you want to see plots as the
pipeline runs or not.  This is controlled by
:ref:`two parameters <configuration:Plotting parameters>`,
``enable_plot_viewer`` and ``plot_level``. For no plotting, set the plot level
to 0. For plotting without interaction, set the plot level to 1 (this is the
default configuration): plots are saved in the ``plots`` directory, and with
``enable_plot_viewer = True`` you can also watch them in a browser, at the
address the pipeline prints when it starts. For interactive plotting (where you
have to hit <cr> to proceed after each step) set the plot level to 2 or more.

Before reducing, you can check that every science frame has the calibrations
it needs:

.. code-block:: bash

   cd mydata
   check_cals kb*.fits

Next, go to the data directory and run the startup script:

.. code-block:: bash

   cd mydata
   reduce_kcwi -b -f kb*.fits -g

The Keck DRP Framework will load and initialize the pipeline, ingest all the
Blue channel files, and then start processing according to the image types in
groups in the order required by the pipeline.  To do the same for the Red
channel data run:

.. code-block:: bash

   reduce_kcwi -r -f kr*.fits -g

Three directories will be created: a ``redux`` directory with the results of the
reduction, a ``logs`` directory with separate logs for the framework itself
and for the DRP, and a ``plots`` directory containing diagnostic plots.
The file you most likely want is the flux-calibrated cube, ``redux/*_icubes.fits``;
see :doc:`data_products`.



