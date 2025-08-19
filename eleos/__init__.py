"""
Eleos: A utility library for working with planetary spectra and gneerating NEMESIS radiative transfer cores.

Modules available directly under `eleos`:
- constants
- cores
- parsers
- profiles
- results
- shapes
- spectra
- spx
- utils
"""

# Explicitly declare the public API
__all__ = [
    "constants",
    "cores",
    "parsers",
    "profiles",
    "results",
    "shapes",
    "spectra",
    "spx",
    "utils",
]

# Import submodules into package namespace
from . import (
    constants,
    cores,
    parsers,
    profiles,
    results,
    shapes,
    spectra,
    spx,
    utils,
)