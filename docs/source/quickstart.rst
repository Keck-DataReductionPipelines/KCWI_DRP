===========
Quick Start
===========

For users who don't want to know too much about KCWI_DRP (and why not?), here
is a quickstart guide.

The assumption is that you have a directory containing KCWI data, and that the
names of the files are those assigned at the telescope, ``kb*.fits`` and
``kr*.fits``.

Give a quick look at the configuration parameters for the pipeline, contained
in the file ``kcwidrp/config/kcwi.cfg``.

For a quick start, it is enough to decide if you want to see plots as the
pipeline runs or not.  This is controlled by two parameters: ``enable_bokeh``
and ``plot_level``. For no plotting, disable Bokeh and set the plot level to 0.
For plotting without interaction, enable Bokeh and set the plot level to 1.
For interactive plotting (where you have to hit <cr> to proceed after each step)
set the plot level to 2 or more.

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



