"""Decide from the configuration whether the live plot viewer is wanted.

This module deliberately imports neither Matplotlib nor Bokeh, so the launch
scripts can use it before anything heavy is loaded.

Rules (bokeh plan 5.1):

* ``enable_plot_viewer`` present: it decides.  ``enable_bokeh`` is then
  ignored, but still gets the deprecation message.
* only ``enable_bokeh`` present: honoured, with the deprecation message.
* neither present: viewer off.

Neither key has any effect on saving PNG files; that follows ``plot_level``.
"""

import logging

_DEPRECATION = ("Config key 'enable_bokeh' is deprecated and will be removed; "
                "use 'enable_plot_viewer' (it controls the live plot viewer "
                "only; saved PNGs follow 'plot_level').")

_deprecation_logged = False

_TRUE = ("true", "yes", "on", "1")


def _present(value):
    """A key is present when it was given a value (not missing / None)."""
    return value is not None


def _as_bool(value):
    """Interpret a parsed config value.

    keckdrpframework's ConfigClass turns the text ``True``/``False`` into
    bools; anything else (``true``, ``yes``, ``1``) arrives as a string or
    number.  Strings are accepted case-insensitively.
    """
    if isinstance(value, str):
        return value.strip().lower() in _TRUE
    return bool(value)


def viewer_requested(config_instrument, logger=None):
    """Return True when the configuration asks for the live plot viewer.

    ``config_instrument`` is ``framework.config.instrument``.  Missing keys
    read as ``None`` there (ConfigClass) and via ``getattr`` default here.
    The ``enable_bokeh`` deprecation message is logged at warning level once
    per process.
    """
    global _deprecation_logged

    new = getattr(config_instrument, 'enable_plot_viewer', None)
    old = getattr(config_instrument, 'enable_bokeh', None)

    if _present(old) and not _deprecation_logged:
        _deprecation_logged = True
        (logger or logging.getLogger('KCWI')).warning(_DEPRECATION)

    if _present(new):
        return _as_bool(new)
    if _present(old):
        return _as_bool(old)
    return False
