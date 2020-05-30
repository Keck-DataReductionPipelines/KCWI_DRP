==================
The startup script
==================

The DRP startup script contains several command line options and different execution modes.

Process files
-------------

* all files in a directory in the order in which they appear

.. code-block:: shell

   reduce_kcwi -f kb*.fits

* a text file containing the list of files to be processed

.. code-block:: shell

   reduce_kcwi -l input_files.lst


The following option can be used to modify the behaviour of the DRP when processing files
in a directory:

``--groups``  This options forces the DRP to group the files by image type and then processes
them in the correct order so that master calibrations are produced before processing
science images, regardless of the order in which the files appear on disk or are specified
in the input list.

Monitor directories
-------------------

The DRP has the ability of monitoring a specified directory. When files appear, they are
ingested and processed. To start the DRP in this mode use:

.. code-block:: shell

   reduce_kcwi -d /home/mydata -i kb*.fits -m

The ``-i kb*.fits`` is the filter used to recognize the correct files. If it is not specified
the pipeline will ingest all files appearing in a directory, and might fail if those files
are not KCWI frames.

Other command line options
--------------------------

* ``-c config_file.cfg``  This options overrides the standard configuration file that is stored
in the installation directory in ``kcwidrp/config/kcwi.cfg``.

* ``-p proctable.proc``  When the DRP runs, it keeps track of the files processed using
a processing table. Normally that table is called ``kcwi.proc`` and is stored in the
current directory. This options is used to specify a different file if needed.

* ``-t taper_fraction``  TBD




