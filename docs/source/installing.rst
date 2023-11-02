===================
Installing KCWI_DRP
===================

This document describes how to install KCWI_DRP, for both users and developers.

If you encounter any issues beyond what is described on this page, see the
:ref:`Issues <issues:FAQ / Known Issues>` page. If nothing there helps you, see
:ref:`help <help:Support>`.

Environment and External Dependencies
=====================================

Installing Dependencies
=======================

We highly recommend that you use `Anaconda <https://www.anaconda.com/>`_ for the
majority of these installations.

Creating a Conda environment
----------------------------

The following will create a :code:`conda` environment called :code:`kcwidrp`,
and activate it.

    .. code-block:: bash

        conda create --name kcwidrp python=3.7
        conda activate kcwidrp

Getting and Installing the Code
-------------------------------

Now, we clone the repository from GitHub and use pip to install the pipeline.
The :code:`pip` step will install all the requirements needed for the pipeline.

    .. code-block:: bash

        git clone https://github.com/Keck-DataReductionPipelines/KCWI_DRP.git
        cd KCWI_DRP
        git checkout kcrm_merge
        pip install .


Getting and Installing the DRP
==============================

    We highly recommend that you use Anaconda for the majority
    of these installations. 

    Detailed installation instructions are presented below:

    Installing with environment.yml
    -------------------------------
    An environment.yml file is provided
    `here <https://github.com/Keck-DataReductionPipelines/KCWI_DRP/blob/kcrm_merge/environment.yml>`_
    which contains the majority of the required dependencies. To create the conda
    environment, download the environment file and run

    .. code-block:: bash

        conda env create -f environment.yml
        conda activate kcwidrp
        pip install kcwidrp

    This creates an environment called kcwidrp that contains most of the required 
    dependencies. 

Installing From Source
----------------------

To install directly from the source code, see the :ref:`developer install<development:Installing for Development>`
instructions.
