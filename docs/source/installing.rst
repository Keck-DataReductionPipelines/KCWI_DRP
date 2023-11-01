===================
Installing KCWI_DRP
===================

This document describes how to install KCWI_DRP, for both users and developers.

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


Installing Manually
-------------------

In some cases, it may be preferable to install the code manually, such as when
you are actively working on the pipeline. You can install the pipeline directly
from source by doing the following:

    .. code-block:: bash

        git clone https://github.com/Keck-DataReductionPipelines/KCWI_DRP.git
        cd KCWI_DRP
        python setup.py develop

Using :code:`python setup.py develop` makes the installation editable. This means
that any changes you make to the code will take effect immediately. Note that
installing using :code:`python setup.py` directly can lead to conflicts if you
later install using ``pip``. If you need to install using ``pip`` later, it is
highly recommended that you create a new conda environment first.

If you do not plan on editing the code, instead use :code:`pip install .` from
inside the source code directory.
