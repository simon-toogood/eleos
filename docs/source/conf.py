# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import datetime
import os
import sys
from pathlib import Path

import sphinx_rtd_theme

print()

sys.path.insert(0, Path("../../eleos").resolve())  # Adjust if needed

import planetmapper
from planetmapper.common import __version__

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'eleos'
copyright = '2024, Simon Toogood'
author = 'Simon Toogood'
release = 'Alpha'
version = __version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx_rtd_theme',
    'sphinx.ext.autosummary',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

default_role = 'code'
master_doc = 'index'

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'matplotlib': ('https://matplotlib.org/stable', None),
    'numpy': ('https://numpy.org/doc/stable', None),
}

# Autodoc
autodoc_member_order = 'bysource'
# autoclass_content = 'both'
# autodoc_typehints = 'both'
# autodoc_typehints_description_target = 'documented_params'
autodoc_inherit_docstrings = False


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
# html_extra_path = ['google065a4d650d8ee82d.html']
# html_logo = 'images/logo_rtd_transparent.png'
html_theme_options = {
    'logo_only': True,
}