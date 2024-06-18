# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))  # Adjust as needed to locate your package


# -- Project information -----------------------------------------------------

project = 'Project Name'
author = 'Author Name'
release = '0.1.0'  # The full version, including alpha/beta/rc tags

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.todo',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.githubpages',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- Extension configuration -------------------------------------------------

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True

# Intersphinx mapping
intersphinx_mapping = {'https://docs.python.org/3/': None}

# Todo extension settings
todo_include_todos = True

# -- Options for autodoc -----------------------------------------------------
autodoc_member_order = 'bysource'

# -- Options for LaTeX output ------------------------------------------------
latex_elements = {
    'papersize': 'a4paper',
    'pointsize': '10pt',
    'preamble': '',
}

# -- Custom configurations ---------------------------------------------------
# Add any other configurations you need here
