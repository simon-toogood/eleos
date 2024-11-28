# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sys
import os
import sphinx_rtd_theme

project = 'eleos'
copyright = '2024, Simon Toogood'
author = 'Simon Toogood'
release = 'Alpha'

sys.path.insert(0, os.path.abspath('../..'))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',	        # To generate autodocs
    'sphinx.ext.mathjax',           # autodoc with maths
    'sphinx.ext.napoleon',          # For auto-doc configuration
    'sphinx.ext.coverage',
    'sphinx.ext.autosummary',
    "sphinx_rtd_theme"  
]

autosummary_generate = True  # Turn on sphinx.ext.autosummary

templates_path = ['_templates']
exclude_patterns = []

napoleon_google_docstring = True   # Turn on googledoc strings
napoleon_numpy_docstring = False     # Turn off numpydoc strings
napoleon_use_ivar = True 	     # For maths symbology

autoclass_content = 'both'


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

