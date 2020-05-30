============================
The KCWI DRP: basic concepts
============================

While most of the basic algorithms are inherited from the IDL KCWI pipeline (KDERP),
the architecture of this pipeline is completely different.

Most of the underlying architecture is a consequence of using the Keck DRP Framework.
This event-based framework implements a set of rules that connect events to processing code.

To illustrate this concept, we can use a simple event such as a new file on disk.
The KCWI_DRP associates this event to a Python class called ``ingest_file``. This means that
the namespace contains a class or a function with this name.
The processing code contained in the class or the function will then be applied to the file.

In more detail, if the code is a class, the ``_perform`` method of the class will be executed.

In our case, the ``ingest_file`` class looks at the header of the file and ingests it in its
entirety. It then trigger further processing steps based on the image type.

Recipes and primitives
----------------------

The original concept upon which this DRP is based is the separation of code in
primitives and recipes. Primitives are simple pieces of code usually dedicated to a single
operation. An example of a primitive would be ``subtract_bias``, which takes a science frame and
a master bias and performs the subtraction. Recipes are sequences of primitives, and they are
usually associated with a specific image type.
An example of a recipe is ``process_science`` which would contain calls to primitives such as
``subtract_overscan``, ``apply_flat_field`` or ``apply_wavelength_calibration`` and so on.

The Keck DRP Framework doesn't natively implement the concept of recipes, but it
allows to specify a *next* event: the *next* event is an event that should automatically
be triggered when the current event is finished. Using this, and linking individual
events via the *next* event, we can easily build sequences of events which, for all intent
and purposes, are the same as recipes.

These recipes are specified in a pipeline definition file, which can be found in ``kcwidrp/pipelines/kcwi_pipeline.py".
This file contains an ``event_table`` which provides the link between events and processing
code. For each entry, the third field is the *next* event.


