.. _versions:

=================
Previous versions
=================

Version 1.3.3 (2026)
--------------------

* Fixed ``StackFlats``: from 1.2.0 to 1.3.2 each stacked flat was made from N
  copies of the first input flat. Stacked flats now combine every input flat,
  and the stacked uncertainty is rescaled with the data.

Version 1.3.2 (2026)
--------------------

* When ``plot_firefox_compat = True``, Firefox and geckodriver installed
  outside snap are preferred over snap-installed ones.

Version 1.3.1 (2026)
--------------------

* Missing Bokeh sessions are handled better.
* New configuration parameters ``terminate_on_failed_bokeh_start`` and
  ``plot_firefox_compat``, for systems where Bokeh cannot start or cannot find
  Firefox (for example snap-installed Firefox on Ubuntu). See :doc:`issues`.

Version 1.3.0 (2026)
--------------------

* Replaced ``setup.py`` based builds with Poetry and ``pyproject.toml``;
  package resources are found with ``importlib`` instead of ``pkg_resources``.
* Added spectrophotometric standards (CALSPEC standards, GD 153); removed
  HD 19445, HD 60753, HZ 29 and HZ 43.
* Fixed the Tukey taper being applied repeatedly to the continuum bars during
  the dispersion calculation, and a mix-up between the taper fraction and the
  Tukey alpha.
* The ``NOSKYSUB`` extension is written as 32-bit floats and carries the WCS.
* ``-v`` shows debug messages in the log.
* Biases are no longer required to share a ``GROUPID`` to be combined.

Version 1.2.0 (2024)
--------------------

* Uses Python 3.12, which allows reductions on Apple Silicon (M-series) Macs as
  well as Intel and AMD machines.
* Added the ``check_cals`` script, which checks whether you have the
  calibration files needed before starting a full reduction.
* Input files can have arbitrary names (they no longer need to keep the names
  they were written with at the telescope).
* Two instances of the pipeline can run at the same time, for example one per
  channel.
* Added standard star spectra for flux calibration.


Version 1.1.0 (2023)
--------------------
Added support for red-side data reduction.

Version 1.0 (2021)
------------------

* Simplified installation via pip and conda environment
* Vacuum to air and heliocentric or barycentric correction (the algorithms used here are courtesy of Yuguang Chen at Caltech)
* Ability of using KOA file names or original file names
* Better provenance and traceability of DRP versions and execution steps in the headers
* Versatile sky subtraction modes including using external sky frames and ability of masking regions


Version 0.1 (2019)
------------------

* First end-to-end version including all the reduction step of the IDL pipeline
* Three execution modes

   * Reduce all files in a directory in the order in which they appear
   * Reduce all files after grouping them by file type and in the correct order
   * Monitor a directory for new files and reduce them as they appear

* Multi-threading for CPU intensive tasks such as wavelength calibration
* Multi-processing for large datasets

