# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import configparser
import os
import pathlib
import re
import subprocess
import sys
from importlib.metadata import entry_points, version
from importlib.resources import files

from sphinx.util import logging

# -- Project information -----------------------------------------------------

project = 'Keck Cosmic Web Imager DRP'
copyright = '2020-2026, Don Neill, KCWI Team, Luca Rizzi, Max Brodheim'
author = 'Don Neill, KCWI Team, Luca Rizzi, Max Brodheim'

# The full version, including alpha/beta/rc tags, taken from the installed
# kcwidrp package (the build installs this repository with its "docs"
# extra, so the docs always describe the code in the same commit)
release = version('kcwidrp')
version = '.'.join(release.split('.')[:2])


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.viewcode',
              'sphinx.ext.autosectionlabel',
              'myst_parser'
              ]
# User pages are reStructuredText; Markdown (myst-parser) is for the
# developer pages that live in docs/dev/ (docs plan, decision 3)
source_suffix = {'.rst': 'restructuredtext', '.md': 'markdown'}
napoleon_numpy_docstring = True
autosectionlabel_prefix_document = True

# Class documentation should contain *both* the class docstring and
# the __init__ docstring
autoclass_content = "class"


master_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'logo_only': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'style_nav_header_background': 'white',
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_logo = "_static/Keck_Logo_Tag.jpg"


# -- Generated pages and consistency checks ----------------------------------

SOURCE_DIR = pathlib.Path(__file__).parent
logger = logging.getLogger(__name__)

# Scripts whose --help output is included in scripts.rst. The other scripts
# (wb, wr, kcwi_masksky_ds9) do not use argparse, and smart_reduce_kcwi does
# not import in kcwidrp 1.3.x.
HELP_SCRIPTS = ['reduce_kcwi', 'check_cals', 'modhead', 'start_kcwi_rti']


def kcwidrp_scripts():
    """Return {name: 'module:function'} for kcwidrp's console scripts."""
    return {ep.name: ep.value for ep in entry_points(group='console_scripts')
            if ep.value.startswith('kcwidrp.')}


def write_script_help(app):
    """Write each script's --help output to _generated/<name>.txt."""
    out_dir = SOURCE_DIR / '_generated'
    out_dir.mkdir(exist_ok=True)
    scripts = kcwidrp_scripts()
    # Fix the width so the output does not depend on the terminal
    env = dict(os.environ, COLUMNS='80')
    for name in HELP_SCRIPTS:
        module, function = scripts[name].split(':')
        code = (f"import sys; sys.argv = [{name!r}, '--help']; "
                f"from {module} import {function}; {function}()")
        result = subprocess.run([sys.executable, '-c', code], env=env,
                                capture_output=True, text=True)
        if result.returncode != 0 or not result.stdout:
            logger.warning(f"{name} --help failed: {result.stderr}")
        (out_dir / f'{name}.txt').write_text(result.stdout)


def check_scripts_documented(app):
    """Warn about console scripts that scripts.rst does not describe."""
    text = (SOURCE_DIR / 'scripts.rst').read_text()
    for name in kcwidrp_scripts():
        if not re.search(rf'^{re.escape(name)}\n-+$', text, re.MULTILINE):
            logger.warning(f"script {name} has no section in scripts.rst")


def check_config_documented(app):
    """Warn about kcwi.cfg keys that configuration.rst does not mention."""
    config = configparser.ConfigParser()
    config.optionxform = str  # keep the case of keys such as TAPERFRAC
    config.read_string(files('kcwidrp').joinpath('configs/kcwi.cfg').read_text())
    keys = {key for section in config.sections() for key in config[section]}
    text = (SOURCE_DIR / 'configuration.rst').read_text()
    for key in sorted(keys):
        if not re.search(rf'\b{re.escape(key)}\b', text):
            logger.warning(f"kcwi.cfg key {key} is not in configuration.rst")


def setup(app):
    app.connect('builder-inited', write_script_help)
    app.connect('builder-inited', check_scripts_documented)
    app.connect('builder-inited', check_config_documented)
