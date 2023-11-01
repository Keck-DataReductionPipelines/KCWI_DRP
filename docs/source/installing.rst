===================
Installing KCWI_DRP
===================

This document describes how to install KCWI_DRP, for both users and developers.

Creating a Conda environment
----------------------------

We highly recommend that you use `Anaconda <https://www.anaconda.com/>`_ for this
installation. Using a ``conda`` environment ensures that the DRP's requirements
will not conflict with other software and packages you have installed on your
machine, among other benefits.

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

Getting and Installing the Code
-------------------------------

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

If you do not plan on editing the code, instead use :code:`pip install .`
