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

Pipeline Ingests Blue Files as Red, and Vice Versa
--------------------------------------------------

There is a known bug in the Keck DRP framework that prevents two instances of
the DRP from running at the same time. This can be inconvenient, as this means
you cannot reduce red and blue side data at the same time on the same machine.
Until this bug gets fixed, the only solution is to run through your data first
for one set of blue/red data, then the other.

Firefox/geckodriver cannot be found
-----------------------------------

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
    
#. Open :code:`KCWI_DRP/kcwidrp/core/kcwi_plotting.py` in a text editor

#. Add the following import at the top of the file: ::

    from selenium import webdriver

#. Find your firefox installation by executing :code:`which firefox` from the
   command line. Make note of the output. It should look something like
   :code:`/PATH/TO/CONDA/envs/kcwidrp/bin/firefox`
#. Replace the function :code:`save_plot` with the following: ::

    def save_plot(fig, filename=None):
        if filename is None:
            fnam = os.path.join('plots', 'kcwi_drp_plot.png')
        else:
            fnam = os.path.join('plots', filename)

        options = webdriver.FirefoxOptions()

        options.add_argument("--headless")
        options.add_argument("--hide-scrollbars")
        options.add_argument("--force-device-scale-factor=1")
        options.add_argument("--force-color-profile=srgb")
        driver = webdriver.Firefox(firefox_binary="[YOUR/PATH/HERE]",
                                    firefox_options=options)
        export_png(fig, filename=fnam, webdriver=driver)
        driver.close()
        logger.info(">>> Saving to %s" % fnam)

   where :code:`[YOUR/PATH/HERE]` is replaced by the path found in the
   previous step
#. Navigate to the :code:`KCWI_DRP` directory, and run::

        python setup.py install

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