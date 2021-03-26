===================
Installing KCWI_DRP
===================

This document describes how to install KCWI_DRP.
Both for users and developers.

Installing Dependencies
=======================

We highly recommend that you use Anaconda for the majority
of these installations. 

Detailed installation instructions are presented below:

Installing with environment.yml
-------------------------------
An environment.yml file is provided
`here <https://github.com/MNBrod/KCWI_DRP/blob/deploy_pip/environment.yml>`_
which contains the majority of the required dependencies. To create the conda
environment, download the environment file and run

.. code-block:: bash

    conda env create -f environment.yml
    conda activate kcwidrp

This creates an environment called kcwidrp that contains most of the required 
dependencies. 

To finalize the installation, gather `LACosmicX <https://github.com/cmccully/lacosmicx>`_, which
provides cosmic ray rejection:

.. code-block:: bash

   git clone https://github.com/cmccully/lacosmicx.git
   cd lacosmicx
   python setup.py install

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

* `LACosmicX <https://github.com/cmccully/lacosmicx>`_

.. code-block:: bash

   git clone https://github.com/cmccully/lacosmicx.git
   cd lacosmicx
   python setup.py install

Installing for Development
--------------------------

If you want to alter the pipeline, you can install it directly from source by
skipping :code:`pip install kcwidrp` during the requirements section above, and
instead running:

.. code-block:: bash

    git clone https://github.com/Keck-DataReductionPipelines/KCWI_DRP.git
    cd KCWI_DRP
    python setup.py develop


