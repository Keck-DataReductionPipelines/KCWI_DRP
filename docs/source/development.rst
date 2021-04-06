=====================
Development Procedure
=====================

The following document describes the procedure to follow for making changes to
KCWI_DRP.

Managing Git
------------

There are three persistent branches in the KCWI_DRP repo:

- ``master`` (default): This branch contains the current release version of the code, and provides the source used to build the package for ``conda`` and  ``pip``. 
- ``develop``: This branch contains new features and non-essential fixes. It should be stable at all times. A new release is triggered by this branch being merged into ``master``. Any new code should be added to this branch.
- ``documentation``: This branch contains the documentation. Any new features or changes should be recorded here. For instructions on updating the docs, see the Updating Documentation page.

When adding new code, please do your development in a branch off of ``develop``:

.. code-block:: bash

    git checkout develop
    git pull
    git checkout -b new_feature

While developing, it is important to stay up-to-date with any changes in the
``develop`` branch:

.. code-block:: bash

    git checkout develop
    git pull
    git checkout new_feature
    git merge --no-ff develop

Alternatively, you can do your own development in your own fork of the repo.

Pull Requests
^^^^^^^^^^^^^

Once your new feature is ready for merging into the develop branch, you need to
issue a pull request. The pull request must contain the following:

1. A short but thorough explanation of the feature and its purpose 
2. A pass from the automatic TravisCI tests
3. Evidence of successfully processing a night of KCWI data

Additionally, the code must meet the following:

- All code must be documented. This includes updated docstrings for altered code and new ones for new functions or classes
- Comments explaining any unusually complicated or unintuitive code blocks, in addition to the documentation above
- Remove all print statements. If you need to output information to terminal, use the logging interface
- At least two developers must approve of the request for it to be accepted

If your pull request meets all of the above, it will be merged into ``develop``
during the next developer meeting.

Merging ``develop`` into ``master``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At the developer's discretion, ``develop`` will be merged into ``master``, which
will signify a new sub-release. The version number should be incremented in
``develop``, and a summary of changes should be added to CHANGES.rst. Then,
issue and accept a pull request from ``develop`` into ``master``.

Building for Release
--------------------

These instructions summarize the steps required to build a new release of 
kcwidrp for conda or pip. These steps should be followed every time ``master``
is updated.

Pip 
^^^

In order to upload to pip, you will need access to a PyPI account with
owndership status for the kcwidrp project. For access to the KeckDRPs account,
ask `Max <mbrodheim@keck.hawaii.edu>`_ or `Luca <lrizzi@keck.hawaii.edu>`_.

After your pull request is merged into master, download the changes:

.. code-block:: bash

    git checkout master
    git pull

Then, after ensuring the ``dist`` directory is empty (ensure no previous
versions exist; these will conflict with the upload to PyPI):

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
^^^^^

Eventually, these steps will be rendered obsolete by the use of conda-forge. In
the meantime, the following instructions will build a conda package from the pip
package. This should be run whenever a new pip version is created.

.. code-block:: bash

    conda update conda
    conda install conda-build anaconda-client

    conda-build conda_build_files
    conda build conda_build_files --output

    anaconda login
    anaconda upload PATH-FROM-OUTPUT