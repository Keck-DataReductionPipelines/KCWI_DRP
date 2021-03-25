===================
Installing KCWI_DRP
===================

This document describes how to install KCWI_DRP.
Both for users and developers.

Installing Dependencies
=======================

We highly recommend that you use Anaconda for the majority
of these installations. An environment.yml file is provided
`here <https://github.com/MNBrod/KCWI_DRP/blob/deploy_pip/environment.yml>`_
which contains all but one of the required packages. To create the conda
environment, download the environment file and run

.. code-block:: bash

    conda env create -f environment.yml

This will create a conda environment called :code:`kcwidrptest`

Detailed installation instructions are presented below:

Python Dependencies
-------------------

* `python <http://www.python.org>`_ version 3.6 or later
* Python packages

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

* `LACosmicX <https://github.com/cmccully/lacosmicx>`_

.. code-block:: bash

   git clone https://github.com/cmccully/lacosmicx.git
   cd lacosmicx
   python setup.py install

* `Keck DRP Framework <https://github.com/Keck-DataReductionPipelines/KeckDRPFramework>`_

.. code-block:: bash

   git clone https://github.com/Keck-DataReductionPipelines/KeckDRPFramework.git
   cd KeckDRPFramework
   python setup.py install (or develop)

* `KCWI DRP <https://github.com/Keck-DataReductionPipelines/KCWI_DRP>`_

.. code-block:: bash

   git clone https://github.com/Keck-DataReductionPipelines/KCWI_DRP.git
   cd KCWI_DRP
   python setup.py install (or develop)


