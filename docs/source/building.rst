====================
Building For Release
====================

These instructions summarize the steps required to build a new
release of kcwidrp for conda or pip.

Pip
---

In order to upload to pip, you will need access to a PyPI account with
owndership status for the kcwidrp project. For access to the KeckDRPS account,
ask `Luca <lrizzi@keck.hawaii.edu>`_ or `Max <mbrodheim@keck.hawaii.edu>`_.

After your pull request is merged into master, download the changes:

.. code-block:: bash

    git checkout master
    git pull origin master

Then, after ensuring the dist directory is empty:

.. code-block:: bash

    # Make sure you have the most recent version of twine installed
    pip install twine --upgrade

    # Construct the pip distribution
    python setup.py sdist bdist_wheel

    # Test the upload by uploading to TestPyPI
    twine upload --repository-url https://test.pypi.org/legacy/ dist/*

    #If all went well with the test, upload to the permanent PyPI
    twine upload dist/*

Conda
-----

.. code-block:: bash

    conda update conda
    conda install conda-build anaconda-client

    conda-build conda_build_files
    conda build conda_build_files --output

    anaconda login
    anaconda upload PATH-FROM-OUTPUT