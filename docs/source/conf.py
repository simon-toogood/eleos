import os
import sys

# Add project root to path
sys.path.insert(0, os.path.abspath(".."))

project = "Eleos"
author = "Your Name"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",      # Auto API docs from docstrings
    "sphinx.ext.napoleon",     # Support for Google/NumPy style docstrings
    "sphinx.ext.viewcode",     # Add links to source code
]

templates_path = ["_templates"]
exclude_patterns = []

html_theme = "sphinx_rtd_theme"

def skip(app, what, name, obj, would_skip, options):
    if name == "__init__":
        return False
    return would_skip

def setup(app):
    app.connect("autodoc-skip-member", skip)