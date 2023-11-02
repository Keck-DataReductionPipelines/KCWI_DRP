===================
Installing KCWI_DRP
===================

This document describes how to install KCWI_DRP, for both users and developers.

If you encounter any issues beyond what is described on this page, see the
:ref:`Issues <issues:FAQ / Known Issues>` page. If nothing there helps you, see
:ref:`help <help:Support>`.

Environment and External Dependencies
=====================================

This section describes how to set up your Python environment for the KCWI DRP,
and covers a few external dependencies. 

We highly recommend that you use `Anaconda <https://www.anaconda.com/>`_ for this
installation. Using a ``conda`` environment ensures that the DRP's requirements
will not conflict with other software and packages you have installed on your
machine, among other benefits.


Creating the Conda environment
------------------------------

.. warning::
    If you have previously installed the DRP from source, and especially if it
    was a previous version, we strongly advise you to delete the ``kcwidrp`` 
    environment and create a new one. Installation conflicts are likely to occur
    if a stale environment is used.

    To remove a ``conda`` environment, use the following command:

    .. code-block:: bash

        conda remove --name kcwidrp --all

The following will create a :code:`conda` environment called :code:`kcwidrp`,
and activate it.

    .. code-block:: bash

        conda create --name kcwidrp python=3.7
        conda activate kcwidrp


Firefox
-------

The KCWI DRP uses `bokeh <http://bokeh.org/>`_ for plotting, which requires access
to a web browser. We recommend installing `Firefox <https://www.mozilla.org/en-US/firefox/new/>`_
for this purpose, as it is the browser used for testing and validating the DRP.
Other browsers have been reported to work, but we haven't tested them ourselves.

If you are running the DRP on a headless system and run into issues with driver
support, see :ref:`firefox/geckodriver issues <issues:Firefox/geckodriver cannot be found>`.


Getting and Installing the DRP
==============================

The following sections describe how to install the pipeline into your new 
environment. We recommend most users follow the ``pip`` section. If you
plan on editing the pipeline, follow the Installing Manually section.

Installing with ``pip`` (recommended)
-------------------------------------

Next, we install the pipeline using ``pip``. This step will automatically install
all the requirements needed for the pipeline.

    .. code-block:: bash

        pip install kcwidrp


Installing From Source
----------------------

To install directly from the source code, see the :ref:`developer install<development:Installing for Development>`
instructions.
