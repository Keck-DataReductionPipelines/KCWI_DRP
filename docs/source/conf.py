# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

from importlib.metadata import version

# -- Project information -----------------------------------------------------

project = 'Keck Cosmic Web Imager DRP'
copyright = '2020-2026, Don Neill, KCWI Team, Luca Rizzi, Max Brodheim'
author = 'Don Neill, KCWI Team, Luca Rizzi, Max Brodheim'

# The full version, including alpha/beta/rc tags, taken from the installed
# kcwidrp package (see docs/requirements.txt)
release = version('kcwidrp')
version = '.'.join(release.split('.')[:2])


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.napoleon',
              'sphinx.ext.viewcode',
              'sphinx.ext.autosectionlabel'
              ]
napoleon_numpy_docstring = True
autosectionlabel_prefix_document = True

# Class documentation should contain *both* the class docstring and
# the __init__ docstring
autoclass_content = "class"

# kcwi_masksky_ds9 imports regions, which kcwidrp 1.3.x does not list as a
# dependency (it lists pyregion instead), so autodoc could not import it
autodoc_mock_imports = ['regions']


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
