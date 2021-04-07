.. _running: 

====================
Running the pipeline
====================

The DRP is instantiated using a startup script offers contains several command line options and different execution modes.

Process files, file lists and entire directories
------------------------------------------------

* To reduce all files in a directory in the order in which they appear and group them correctly according to the logical sequence needed by the pipeline:
 

.. code-block:: shell

   reduce_kcwi -f kb*.fits -g

* To reduce only a subset of the files (for example, a single object out of an entire observing run):

.. code-block:: shell

   reduce_kcwi -l input_files.lst

* To reduce a single file:

.. code-block:: shell

   reduce_kcwi -f kb001.fits

For the "list" and "single file" cases, please note that the preliminary calibrations needed to reduce the file or the file list must be already on disk.
For example, if the file you are trying to reduce is a science frame, the wavelength calibration, flat field and bias frames must already exist.

Reducing a single file is a good way to re-reduce a target for which the pipeline didn't do a good job. An example case would be 
the sky subtraction: if you realize that the sky subtraction is not correct, and you modify the sky subtraction using an external file, then you can rerun 
the pipeline just on the trouble file. In this case it is advisable to use the ``clobber=True`` option.

The following option can be used to modify the behaviour of the DRP when processing files
in a directory:

``--groups``  This options forces the DRP to group the files by image type and then processes
them in the correct order so that master calibrations are produced before processing
science images, regardless of the order in which the files appear on disk or are specified
in the input list.

Running individual steps
------------------------

This is not a standard operating mode and is not supported but it is possible.
The pipeline has a number of high level procedures that determine the correct order for data reduction.

If you know what you are doing, you can in principle create individual file lists for the steps.

Example:

Create a file containing bias frames and call it ``bias.lst``.
You can can run 

.. code-block: shell

   reduce_kcwi -l bias.lst

Because the files are all and only bias, the pipeline will only proceed as far as generating the master bias.

For an even finer control, each file in the bias list could be run individually, and only when enough bias frames are present and reduced, the DRP will generate a master bias.



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




