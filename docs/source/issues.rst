============
Known Issues
============

This page contains a list of known issues related to the pipeline that are not
currently being worked on. Only issues which cannot be fixed with a new release
will be added to this list. Wherever possible, a link to relevent GitHub issues
will be provided.

Firefox/geckodriver cannot be found
===================================

Error message:

.. code-block:: bash

    RuntimeError: Neither firefox and geckodriver nor a variant of chromium browser and chromedriver are available on system PATH. You can install the former with 'conda install -c conda-forge firefox geckodriver'.

This issue has been submitted by a small number of users. This error is thrown
by our plotting library, bokeh, when it can't automatically find the path to a
driver needed to generate plots. The following describes a one-off workaround.

These instructions assume you are using :code:`conda` to manage your environment.


#. If you installed the pipeline with :code:`pip`, uninstall with 
   :code:`pip uninstall kcwidrp`
#. Install the pipeline following the Installing for Development instructions.
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


