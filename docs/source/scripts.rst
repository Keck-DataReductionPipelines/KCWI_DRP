=======
Scripts
=======

Installing the DRP gives you access to several scripts beyond the main reduction
script, ``reduce_kcwi``. This page details their usage and purpose.


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

.. automodule:: kcwidrp.scripts.wr
   :members:

check_cals
----------

.. code-block:: console

    usage: check_cals [-h] [-v] [-c CONFIG] filepaths [filepaths ...]

    Checks a directory of data to see if the OBJECTs within have the required
    cals.

    positional arguments:
    filepaths             Files to inspect

    optional arguments:
    -h, --help            show this help message and exit
    -v, --verbose         Print exhaustive information
    -c CONFIG, --config CONFIG
                            KCWI configuration file

.. automodule:: kcwidrp.scripts.check_cals
   :members: