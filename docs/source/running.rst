.. _running: 

====================
Running the pipeline
====================

The DRP is run using a startup script that offers several
command line options and different execution modes.  Now that the Red channel
has been added, the user must specify which channel to process with the
``-b`` or ``--blue`` command line flag, or ``-r`` or ``--red`` command line
flag.  If neither channel is specified, the script will not run.

Process files, file lists and entire directories
------------------------------------------------

- To reduce all Blue channel files in a directory in the order in which they
  appear and group them correctly according to the logical sequence needed by
  the pipeline:


.. code-block:: shell

   reduce_kcwi -b -f kb*.fits -g

- The corresponding Red channel command would be:

.. code-block:: shell

   reduce_kcwi -r -f kr*.fits -g

Here ``-b`` or ``-r`` specifies the blue or red channel, ``-f`` specifies file
input and is followed by a file specification, and ``-g`` specifies group mode,
which will group images according to what is needed by the pipeline.

* To reduce only a subset of the files (for example, a single object out of an
  entire observing run):

.. code-block:: shell

   reduce_kcwi -b -l input_files.lst

Here ``-l`` specifies list input mode followed by a file with a list of raw
image files for, in this case, the Blue channel (as indicated by the ``-b``
parameter), one per line in the file.

* To reduce a single file:

.. code-block:: shell

   reduce_kcwi -b -f kb001.fits

For the "list" and "single file" cases, please note that the preliminary
calibrations needed to reduce the file or the file list must be already on
disk.  For example, if the file you are trying to reduce is a science frame,
the wavelength calibration, flat field and bias frames must already exist.

Reducing a single file is a good way to re-reduce a target for which the
pipeline didn't do a good job. An example case would be the sky subtraction:
if you realize that the sky subtraction is not correct, and you modify the sky
subtraction using an external file, then you can rerun the pipeline just on the
trouble file. In this case it is advisable to remove the previously processed
image files in the `redux` directory, and remove that file's entry in the proc
file.

The following option can be used to modify the behaviour of the DRP when
processing files in a directory:

``--groups``  This options forces the DRP to group the files by image type and
then processes them in the correct order so that master calibrations are
produced before processing science images, regardless of the order in which the
files appear on disk or are specified in the input list. This is also accessible
with the ``-g`` flag.

Running individual steps
------------------------

If you know what you are doing, you can in principle create individual file
lists for each reduction step and process them step-wise.  This allows the user
to check the output and quality of the reduction steps prior to proceeding with
the next step.  This is facilitated by the `wr` and `wb` scripts.  If
you run these scripts on your images, they will produce a header listing and
a set of file lists that can be used to process each step individually.  For
example:

.. code-block:: shell

   wb kb*.fits > whatb.list

This will create a list of images (with extension .txt) that group images into
processing blocks.  If you had a set of Blue channel biases, it might produce a
list called `bias2x2TUP010.txt`, for example.  This can then be input as a list
to produce a master bias, which should be the first step in the pipeline
reduction:

.. code-block:: shell

   reduce_kcwi -b -l bias2x2TUP010.txt

Because the files are all and only bias, the pipeline will only proceed as far
as generating the master bias.  Next, you can reduce a list of continuum bars,
which will have a list for each configuration, for example:
`cbars2x2MedKBlueBL_4500_0.7.txt`.  This can be followed by arcs and then flats,
etc.

For an even finer control, each file in the bias list could be run individually,
and only when enough bias frames are present and reduced, the DRP will generate
a master bias.

Monitor directories
-------------------

The DRP has the ability to monitor a specified directory. When files appear,
they are ingested and processed. To start the DRP in this mode use:

.. code-block:: shell

   reduce_kcwi -b -d /home/mydata -i kb*.fits -m

The ``-i kb*.fits`` is the filter used to recognize the correct files. If it is
not specified, the pipeline will ingest all files in the directory, and
will fail if any of those files are not KCWI frames.  The ``-d`` and the
following parameter specify the directory to monitor, and ``-m`` specifies
monitor mode.

Other command line options
--------------------------

* ``-c config_file.cfg``  This options overrides the standard configuration
  file that is stored in the installation directory in
  ``kcwidrp/config/kcwi.cfg``.

* ``-p proctable.proc``  When the DRP runs, it keeps track of the files
  processed using a processing table. Normally that table is called
  ``kcwib.proc`` for the Blue channeld and ``kcwir.proc`` for the Red channel
  and is stored in the current directory. This options is used to
  specify a different file if needed (not recommended).

* ``-t taper_fraction``  This option allows the user to adjust the taper
  fraction that is used to cross-correlate with the atlas spectrum.  It defaults
  to 0.2 and can be increased if there is a strong line near the edge that is
  throwing off the cross-correlation.

* ``-a atlas_line_list`` Specify an input line list for the atlas instead of
  generating it on the fly.

* ``-M middle_fraction`` Specify what central fraction to use for the initial
  estimation of the central dispersion.  It defaults to 1/3 of the wavelength
  range, but can be increased if there aren't enough lines in the default range.

* ``-o atlas_offset`` Specify the atlas offset in pixels to line up the atlas
  and the observed spectrum.  This overrides the value calculated from
  cross-correlating the atlas and observed spectra.

* ``-e line_thresh`` Specify the line cleaning threshold in electrons below
  which to reject lines as too faint.

* ``-u tukey_alpha`` Specify the Tukey taper alpha that is used to cross-correlate
  the bars to each other.

* ``-F line_peak_fraction`` Specify the line fitting window threshold in units
  of the peak.  It defaults to 0.5 (Half-max), but can be either extended or
  narrowed as needed.
