===========
Quick Start
===========

For users who don't want to know too much about KCWI_DRP, here is a quickstart
guide.

The assumption is that you have a directory containing KCWI data, and that the names
of the files are those assigned at the telescope, ``kb*.fits``.

Give a quick look at the configuration parameters for the pipeline, contained in the ``kcwidrp/config/kcwi.cfg``.

For a quick start, it is enough to decide if you want to see plots as the pipeline runs or not.
This is controlled by two parameters: ``enable_bokeh`` and ``plot_level``. For no plotting, disable Bokeh and
set the plot level to 0.

Next, go to the data directory and run the startup script:

.. code-block:: bash

   cd mydata
   reduce_kcwi -f kb*.fits

The Keck DRP Framework will load and initialize the pipeline, ingest all the files, and then
start processing in the order in which they appear on disk.

Two directories will be created: a ``redux`` directory with the results of the reduction, and a ``logs`` directory
with separate logs for the framework itself and for the DRP.



