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

Python Dependencies
-------------------

* `python <http://www.python.org>`_ version 3.6 or later
* Python packages

.. code-block:: bash

   conda install bokeh
   conda install selenium
   conda install -c conda-forge phantomjs
   conda install -c astropy ccdproc
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



