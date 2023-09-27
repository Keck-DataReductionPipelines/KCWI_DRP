======================
Updating documentation
======================

The documentation is generated automatically by ReadTheDocs, and triggerd by
pushing commits to the ``main`` branch of the
`KCWI_DRP-Docs <https://github.com/Keck-DataReductionPipelines/KCWI_DRP-Docs/tree/main>`_ repository.

To update the documentation, first clone the documentation repo:

.. code-block:: bash

    git clone https://github.com/Keck-DataReductionPipelines/KCWI_DRP-Docs.git
    
If you already have a local copy of the repo, instead make sure you have the latest version:

.. code-block:: bash

   cd KCWI_DRP-Docs
   git checkout main
   git pull



Updating text
=============

The documentation is contained in ``.rst`` files in the ``documentation``
branch. New documents can be added to the source directory, and linked to the
correct section (usually in ``index.rst``).

Any non-``.rst`` files should be placed in a sub-directory. For example, ``.csv``
files that describe the data outputs are inside a ``outputs`` directory. Screenshots
or plots should be placed inside the ``_static`` directory.

Once the new document (or a modified document) is ready, follow this procedure:


.. code-block:: bash

   git checkout main
   git pull
   git add new_document.rst
   git commit -m "Description of the new document or change"
   git push

Legacy Documentation
====================

It may be the case that specific versions of the documentation should be kept 
around after the codebase has changed. For example, a copy of the blue-side-only
documentation is permanently available.

To preserve a locked version of the documentation, the documentation source code
should be put into a separate branch and tagged on GitHub. Then, the `DSI Team <dsi-team@keck.hawaii.edi>`_
can configure ReadTheDocs to display that branch as a version option.

