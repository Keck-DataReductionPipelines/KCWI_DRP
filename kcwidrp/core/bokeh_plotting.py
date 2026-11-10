"""
Hand-off from the primitives to the live plot viewer.

The names are kept from the Bokeh era so primitives do not change: the
viewer is a :class:`kcwidrp.core.plot_viewer.PlotViewer`, still stored as
``context.bokeh_session`` (None when no viewer runs).

Created on Aug 8, 2019

@author: lrizzi
"""
import logging

import psutil

logger = logging.getLogger('KCWI')


def bokeh_plot(plot, session):
    """Show a recorded plot in the viewer ``session``.

    Returns at once; does nothing when ``session`` is None; never raises.
    """
    if session is None:
        return
    try:
        session.show(plot)
    except Exception as e:
        logger.warning("Could not send the plot to the plot viewer: %s", e)
        logger.debug("bokeh_plot failure", exc_info=True)


def bokeh_clear(session):
    """Show an empty figure in the viewer ``session``.

    Returns at once; does nothing when ``session`` is None; never raises.
    """
    if session is None:
        return
    try:
        session.clear()
    except Exception as e:
        logger.warning("Could not clear the plot viewer: %s", e)
        logger.debug("bokeh_clear failure", exc_info=True)


def init_plot_session(context):
    """Ensure context.bokeh_session exists (None if no viewer is running).

    Primitives read context.bokeh_session unconditionally, but
    StartPlotViewer only sets it when the viewer is enabled (B9).
    """
    if not hasattr(context, 'bokeh_session'):
        context.bokeh_session = None


def check_running_process(process=None):
    """
    Check if there is any running process that contains
    the given name processName.
    """
    # Iterate over the all the running process
    for proc in psutil.process_iter(attrs=["name", "cmdline"]):
        try:
            # Check if process name contains the given name string.
            for command in proc.cmdline():
                if process in command:
                    return True
        except (psutil.NoSuchProcess, psutil.AccessDenied,
                psutil.ZombieProcess):
            pass
    return False
