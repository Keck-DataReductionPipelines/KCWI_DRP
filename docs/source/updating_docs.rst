======================
Updating documentation
======================

The documentation lives in the KCWI_DRP repository, next to the code it
describes:

- ``docs/source/``: the user documentation, in reStructuredText (``.rst``),
  built by Sphinx and published on Read the Docs.
- ``docs/dev/``: developer and planning documents, in Markdown (``.md``). Some
  of them, such as the :doc:`testing guide <testing>`, are included in the
  site.

Until version 3 is released, Read the Docs still builds the released 1.3.x
documentation from the former
`KCWI_DRP-Docs <https://github.com/Keck-DataReductionPipelines/KCWI_DRP-Docs>`_
repository. From the version 3 release on, it builds from KCWI_DRP: ``latest``
from ``develop`` and ``stable`` from the newest release tag.

Changing the documentation
==========================

Documentation changes go in the same pull request as the code they describe,
on the same branch (see :doc:`development`). A pull request that changes
user-visible behavior, a configuration key, a console script or an output file
updates ``docs/source/`` too.

New pages are added to ``docs/source/`` and linked from a table of contents
(usually in ``index.rst``). Non-``.rst`` files go in a sub-directory: for
example, the ``.csv`` tables that describe the data outputs are in
``outputs/``, and screenshots or plots go in ``_static/``.

Building the documentation
==========================

Build the documentation locally before opening the pull request, the same way
Read the Docs does:

.. code-block:: bash

   pip install -e ".[docs]"
   sphinx-build -W --keep-going docs/source docs/_build

Open ``docs/_build/index.html`` to read the result. The build fails on any
warning. Besides Sphinx's own warnings, this catches documentation that has
drifted from the code:

- a module that the API pages cannot import;
- a key in ``kcwidrp/configs/kcwi.cfg`` that ``configuration.rst`` does not
  mention;
- a console script with no section in ``scripts.rst``.

Paste the last line of the build into the pull request description.

Legacy documentation
====================

The documentation of older versions (1.0.1 to 1.3.x) will stay readable after
the switch to this repository: at the version 3 release it moves to a separate
Read the Docs project built from the archived KCWI_DRP-Docs repository, linked
from :doc:`versions`.
