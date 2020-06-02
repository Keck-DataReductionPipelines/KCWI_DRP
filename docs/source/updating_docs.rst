======================
Updating documentation
======================

The documenation is generated automatically by ReadTheDocs, and triggerd by
pushing commits to the ``documentation`` branch.

Updating text
=============

The documentation is contained in ``.rst`` files in the ``documentation`` branch.
The ``docs`` subdirectory is only present in this branch. New documents can be added
to the source directory, and linked to the correct section (usually in ``index.rst``).

Once the new document (or a modified document) is ready, follow this procedure:

.. code-block:: bash

   git pull
   git add new_document.rst
   git commit -m "I have added a new document"
   git push

Keeping the documentation branch up to date
===========================================

Because the documentation uses autodoc methods (such as ``automodule``) to extract
doc strings from the functions and classes, it is important to keep the ``documenation``
branch in sync with the master branch (or with any other active branch).
To to this, follow this procedure:

.. code-block:: bash

   git checkout documentation
   git pull
   git rebase master (or develop or other active branches)
   git push --force

