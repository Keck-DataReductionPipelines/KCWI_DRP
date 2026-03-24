import numpy as np
from bokeh.io import export_png

import atexit
import os
import shutil
import logging

from selenium import webdriver

logger = logging.getLogger('KCWI')

_firefox_driver = None


def configure_plot_driver(firefox_compat=False, prewarm=False):
    """Optionally pre-warm the Firefox WebDriver at pipeline startup.

    Called by StartBokeh when both plot_firefox_compat and plot_prewarm_firefox
    are True in kcwi.cfg.  Once the driver is created here, save_plot will
    automatically use the Firefox-compatible export path for all subsequent
    calls.
    """
    if firefox_compat and prewarm:
        logger.info("Pre-warming Firefox WebDriver for plot export")
        _get_driver()


def _get_driver():
    """Return a cached Firefox WebDriver, creating it on first call."""
    global _firefox_driver
    if _firefox_driver is not None:
        try:
            # Navigate to about:blank to cancel any pending navigation left by
            # Bokeh's reset step (it navigates to http://0.0.0.1/ after each
            # screenshot; on some networks that address times out rather than
            # refusing immediately, blocking the next export call).
            _firefox_driver.get("about:blank")
            return _firefox_driver
        except Exception:
            _firefox_driver = None

    from selenium.webdriver.firefox.service import Service as FirefoxService
    options = webdriver.FirefoxOptions()
    options.add_argument("--headless")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--hide-scrollbars")
    options.add_argument("--force-device-scale-factor=1")
    options.add_argument("--force-color-profile=srgb")
    firefox_bin = _find_firefox_binary()
    if firefox_bin:
        options.binary_location = firefox_bin
    geckodriver = _find_geckodriver()
    serv = FirefoxService(executable_path=geckodriver) if geckodriver else FirefoxService()
    _firefox_driver = webdriver.Firefox(options=options, service=serv)
    atexit.register(_close_driver)
    return _firefox_driver


def _close_driver():
    global _firefox_driver
    if _firefox_driver is not None:
        try:
            _firefox_driver.quit()
        except Exception:
            pass
        _firefox_driver = None


def get_plot_lims(data, padding=0.05, clip=True):
    """Get plot limits using data range plus padding fraction"""
    dmin = np.nanmin(data)
    dmax = np.nanmax(data)
    ddel = np.abs(dmax - dmin)
    if ddel < 0.0001 and clip:
        ddel = 1.0
    return dmin - padding * ddel, dmax + padding * ddel


def oplot_slices(fig, yrange):
    """Overplot slices vertically on plot"""
    for ix in range(1, 24):
        sx = ix * 5 - 0.5
        fig.line([sx, sx], yrange, color='black', line_dash='dashdot')


def set_plot_lims(fig, xlim=None, ylim=None):
    """Set bokeh figure plot ranges"""
    if xlim:
        fig.x_range.start = xlim[0]
        fig.x_range.end = xlim[1]
    if ylim:
        fig.y_range.start = ylim[0]
        fig.y_range.end = ylim[1]

def _find_firefox_binary():
    """Return path to the real Firefox binary, handling snap wrapper installs."""
    # On Ubuntu 24.04, /usr/bin/firefox is a shell script wrapper for the snap.
    # Selenium requires the actual ELF binary, not a wrapper script.
    candidates = [
        '/snap/firefox/current/usr/lib/firefox/firefox',
        '/usr/lib/firefox/firefox',
    ]
    for path in candidates:
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path
    # Fall back to whatever is on PATH (works on most non-snap installs)
    return shutil.which('firefox')


def _find_geckodriver():
    """Return path to geckodriver, searching common locations."""
    candidates = [
        '/snap/bin/geckodriver',
        '/usr/bin/geckodriver',
        '/usr/local/bin/geckodriver',
    ]
    for path in candidates:
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path
    return shutil.which('geckodriver')


def save_plot(fig, filename=None):
    """Save a Bokeh figure to a PNG file.

    Dispatches to the Firefox-compatible path or the simple path depending on
    whether configure_plot_driver() was called with firefox_compat=True.
    """
    if _firefox_driver is not None:
        _save_plot_firefox(fig, filename)
    else:
        _save_plot_simple(fig, filename)


def _save_plot_firefox(fig, filename=None):
    """Firefox-compatible export for systems with snap-confined Firefox."""
    if filename is None:
        fnam = os.path.join('plots', 'kcwi_drp_plot.png')
    else:
        fnam = os.path.join('plots', filename)

    # Resolve to absolute path before chdir so the PNG lands in the right place.
    fnam = os.path.abspath(fnam)
    os.makedirs(os.path.dirname(fnam), exist_ok=True)

    driver = _get_driver()
    # Snap-confined Firefox on Ubuntu 24.04 uses a private /tmp mount namespace,
    # so file:///tmp/... URLs are invisible to it.  $HOME is accessible via the
    # snap 'home' interface.  Bokeh writes its scratch HTML relative to the
    # process CWD, so chdir to $HOME to land the temp file somewhere Firefox
    # can actually reach.  fnam is already absolute so the PNG output goes to
    # the right place regardless.
    old_cwd = os.getcwd()
    try:
        os.chdir(os.path.expanduser('~'))
        export_png(fig, filename=fnam, webdriver=driver)
    finally:
        os.chdir(old_cwd)
    logger.info(">>> Saving to %s" % fnam)


def _save_plot_simple(fig, filename=None):
    """Simple export using Bokeh's built-in browser detection."""
    if filename is None:
        fnam = os.path.join('plots', 'kcwi_drp_plot.png')
    else:
        fnam = os.path.join('plots', filename)
    export_png(fig, filename=fnam)
    logger.info(">>> Saving to %s" % fnam)
