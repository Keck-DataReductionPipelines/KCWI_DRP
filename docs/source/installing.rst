===================
Installing KCWI_DRP
===================

This document describes how to install KCWI_DRP, for both users and developers.

Since version 1.1 is still not officially released, users MUST install the DRP
directly from GitHub. These instructions will change once the DRP is released to
PyPI, which will allow users to install the DRP with :code:`pip`.

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

..
    Installing Dependencies
    =======================

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


    Installing Manually
    -------------------

    This pipeline currently runs on `python <http://www.python.org>`_ 3.7.
    Instructions for installing the other dependencies are below:

    .. code-block:: bash

    conda install bokeh
    conda install -c conda-forge selenium geckodriver firefox phantomjs
    conda install -c astropy ccdproc pyregion
    conda install psutil
    conda install requests
    conda install pytest
    conda install cython
    conda install pandas
    pip install ref_index
    pip install keckdrpframework
    pip install kcwidrp


    Installing for Development
    --------------------------

    If you want to alter the pipeline, you can install it directly from source by
    skipping :code:`pip install kcwidrp` during the requirements section above, and
    instead running:

    .. code-block:: bash

        git clone https://github.com/Keck-DataReductionPipelines/KCWI_DRP.git
        cd KCWI_DRP
        python setup.py develop
