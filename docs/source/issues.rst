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
notably an air-to-vacuum wavelength correction, and a barycentric wavelength
correction. These, predictably, result in small but noticable offsets between 
new and old cubes.

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

Known Issues
============

Bokeh Issues
------------

The DRP relies on a package called ``Bokeh`` to render plots. This package has a
large number of dependencies (both Python and system) that can cause issues. If
you are having issues plotting, read on.

As a first check, try to kill any hanging ``bokeh`` processes:

.. code-block:: bash

    pkill bokeh

Try running the pipeline again, sometimes this is enough to clear any issues.

Failing that, add/update the following parameter to you ``kcwi.cfg`` config file:

.. code-block:: bash

    terminate_on_failed_bokeh_start = True

This tells the pipeline to immediately exit if Bokeh fails to initialize. When
you run the pipeline again, if it makes it immediately fails at ``StartBokeh``,
that confirms that Bokeh is the problem.

Next, try adding/updating the following to your config:

.. code-block:: bash

    plot_firefox_compat = True

This will try to automatically find a firefox installation and use it - this is
especially an issue if your machine is running firefox installed via snap on
Ubuntu, although it may well manifest on other machines. You may need to follow
steps 1 through 3 in the next section, as well.

The next time you run the pipeline, it may appear to hang at the first plotting
call - wait to see if it proceeds. It may take up to a minute or two
for firefox to sort itself out (there's a timeout somewhere that must be waited
out).

You may see a scary looking warning like

.. code-block:: bash

    [Parent 808986, Main Thread] WARNING: Failed to mkdir /home/kcwidrp/snap/firefox/7967/.config/ibus/bus: Not a directory: 'glib warning', file /build/firefox/parts/firefox/build/toolkit/xre/nsSigHandlers.cpp:201

This is not something to worry about.

Firefox/geckodriver cannot be found
+++++++++++++++++++++++++++++++++++

Error message:

.. code-block:: bash

    RuntimeError: Neither firefox and geckodriver nor a variant of chromium browser and chromedriver are available on system PATH. You can install the former with 'conda install -c conda-forge firefox geckodriver'.

This issue has been submitted by a small number of users. This error is thrown
by our plotting library, bokeh, when it can't automatically find the path to a
driver needed to generate plots. The following describes a one-off workaround.

These instructions assume you are using :code:`conda` to manage your environment.


#. If you installed the pipeline with :code:`pip`, uninstall with 
   :code:`pip uninstall kcwidrp`
#. Install the pipeline following the Installing for Development instructions on
   the :doc:`installing` page.
#. Install the following two packages: ::

    pip install selenium
    pip install geckodriver

Add or update the following in your config file:

.. code-block:: bash

    plot_firefox_compat = True

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


PyQt5
-----

Some users have reported an issue where ``PyQt5`` is required to run the DRP, which
looks like 

.. code-block:: console

    Failed to import any Qt binding

If this is the case, run

.. code-block:: bash

    pip install pyqt5

and the issue should be fixed.