=======
Scripts
=======

Installing the DRP gives you the following commands. The usage text on this
page is generated from each command's ``--help`` when the documentation is
built, so it matches the installed version.

reduce_kcwi
-----------

The main reduction script. :doc:`running` explains how to use it; the full
list of options is:

.. literalinclude:: _generated/reduce_kcwi.txt
    :language: none

The options ``-t``, ``-a``, ``-M``, ``-o``, ``-e``, ``-u`` and ``-F`` override
the wavelength fitting parameters in the configuration file (see
:ref:`running:Other command line options`). The options ``-i``, ``-d``,
``-m``, ``-I``, ``-w``, ``-W`` and ``-s`` control how the Keck DRP Framework
ingests files and waits for events; they are used for real-time reduction at
the telescope and are not needed for normal use.

check_cals
----------

Checks a set of raw files to see whether each ``OBJECT`` frame has the
calibrations it needs, given the minimum numbers of frames in the
configuration file. Run it on your data before a full reduction:

.. code-block:: bash

    check_cals kb*.fits

.. literalinclude:: _generated/check_cals.txt
    :language: none

.. automodule:: kcwidrp.scripts.check_cals
   :members:

wb
--

.. code-block:: console

    usage: wb <fspec>

.. automodule:: kcwidrp.scripts.wb
   :members:

wr
--

.. code-block:: console

    usage: wr <fspec>

The Red channel version of ``wb``.

.. automodule:: kcwidrp.scripts.wr
   :members:

kcwi_masksky_ds9
----------------

.. code-block:: console

    usage: kcwi_masksky_ds9 <image>_intf.fits <regions>.reg

Makes a sky mask from a ds9 region file. See
:ref:`sky_subtraction:Building a mask file`.

.. automodule:: kcwidrp.scripts.kcwi_masksky_ds9
   :members:

modhead
-------

Sets a FITS header keyword in a file. The value is always written as a
string. Without ``-s`` the file is changed in place; with ``-s`` the original
is kept as ``<infile>~``.

.. code-block:: bash

    modhead kb230605_00120.fits -k OBJECT -v feige34

.. literalinclude:: _generated/modhead.txt
    :language: none

start_kcwi_rti
--------------

Runs the real-time reduction used at WMKO, which sends its results to the Keck
Observatory Archive. It is not intended for general use.

.. literalinclude:: _generated/start_kcwi_rti.txt
    :language: none

smart_reduce_kcwi
-----------------

A test script for the Keck DRP Framework. In 1.3.x it does not run (it fails
with an ``ImportError``); use ``reduce_kcwi``.
