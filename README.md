# Eleos - A Python interface to NEMESIS

Please note, this is a WIP project and features may be added/removed at any time. It is absolutely *not* guarunteed to be backward compatible until it has matured significantly

## Documentation

Documentation can be found on ReadTheDocs [here](https://eleos.readthedocs.io/en/latest/). This is still WIP, full documentation will be added in due course.

## Installation

Eleos is available on PyPI so it can be installed as any other Python package. On Unix-like systems:

`pip install nemesis_eleos`

And on Windows:

`py -m pip install nemesis_eleos`


While `eleos` only creates file for and reads files created by NEMESIS, it is mandatory to have NEMESIS installed for this package to work. This is primarily due to required the utilities `Makephase` and `Normxsc` on the PATH. See the NEMESIS GitHub page for full instructions on [how to download the software](https://github.com/nemesiscode/radtrancode/blob/master/README.md) and [how to compile it](https://github.com/nemesiscode/radtrancode/blob/master/AACOMPILE.txt).


This library was built with the intention of running on the University of Leicester’s HPC system ALICE3. Therefore, functions like `cores.generate_alice_job` and `cores.run_alice_job` are only guaranteed to work on ALICE3, which uses the SLURM job scheduler. Other HPC facilities will require their own template submission files in `data/statics` and functions in `cores`.

## Working Directory Structure

This library will function best with the following structure in your working directory:

```
parent_directory
|-- core_1
|-- core_2
|   ...
generate.py
analyse.py
```

where `parent_directory/` contains a set of cores for a retrieval, `core_N/` is the core folder that NEMESIS is run from, `generate.py` is the code that generates these cores using `eleos.cores`, and `analyse.py` is the code that analyses the results using `eleos.results`. This is only a recommendation, and you can structure your working directories however you like.

## Core Generation Example (generate.py)

This code generates 4 cores, each with a different forward modelling error factor. Scattering is on by default, so these consider both thermal emission from the planet and reflected sunlight.

```python

from eleos import shapes, profiles, cores

# ------ Gas Profile definitions ------ #

ph3 = profiles.GasProfile(gas_name="PH3", 
                          shape=shapes.Shape2(
                                scale_factor=1, scale_factor_error=0.5))

nh3 = profiles.GasProfile(gas_name="NH3", 
                          shape=shapes.Shape2(
                                scale_factor=1, scale_factor_error=0.5))

# ------ Tropospheric clouds ------ #

aero1 = profiles.AerosolProfile(label="Troposphere",
                                shape=shapes.Shape48(
                                    base_pressure=1, base_pressure_error=1,
                                    top_pressure=0.5, top_pressure_error=0.5,
                                    opacity=5.0,       opacity_error=5.0,
                                    fsh=0.8,           fsh_error=0.8))
                
n1 = profiles.ImagRefractiveIndexProfile(label="Troposphere",
                                         shape=shapes.Shape444(
                                             radius=2,                   radius_error=1,
                                             variance=0.1,               variance_error=0.1,
                                             refractive_index=1.3+1e-3j, refractive_index_error=1e-3))

# ------ Stratospheric haze ------ #

aero2 = profiles.AerosolProfile(label="Stratosphere",
                                shape=shapes.Shape48(
                                    base_pressure=0.1, base_pressure_error=0.1,
                                    top_pressure=0.0001, top_pressure_error=0.0001,
                                    opacity=2.0,       opacity_error=2.0,
                                    fsh=0.7,           fsh_error=0.7))

# ------ Core creation ------ #

# This is the directory all the cores will live in (ie. it will generate nemesis/example/core_1/nemesis.apr,...
cd = "nemesis/example/"

# It is recommended to clear the parent directory before running
cores.clear_parent_directory(cd)

# For a full list of options to the NemesisCore constructor, see the documentation
core = cores.NemesisCore(parent_directory=cd,
                         spx_file="data/zonal_spectra/sparse_55.0degS.spx",
                         profiles=[ph3, nh3],
                         fmerror_factor=10)

# Multiple aerosol modes can be added, with either a 444 profile or by specifying parameters (if they are not to be fitted).
core.add_aerosol_mode(aero1, n1)
core.add_aerosol_mode(aero2, radius=0.1, variance=0.1, refractive_index=1.3+1e-3j)

# You can set the error in regions of the spectrum to be very small so NEMESIS will always fit the model there
core.fix_peak(4.5, 0.2)

# Create the files required for NEMESIS to run
core.generate_core()

# Generate a SLURM submission script and submit it to the scheduler
cores.generate_alice_job(cores=core, username="none", hours=1)
cores.run_alice_job(cd)
```

## Result Analysis Example (analyse.py)

This code takes the result of running NEMESIS on the output of the above example and plots the retrieved spectrum and chi squared values then saves it to a file.

```python
import matplotlib.pyplot as plt
from eleos import results

# Read in the core after NEMESIS has been run successfully
res = results.NemesisResult("nemesis/example/core_1/")

# Create a new Figure object with two Axes
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))

# Plot the model and measured spectrum on one axis
res.plot_spectrum(ax=ax1)
res.plot_chisq(ax=ax2)

# Save the figure
plt.tight_layout()
fig.savefig("example.png", dpi=500)
```


## Command Line

Eleos can be run at the command line in order to quickly generate summary plots and print a human-readable representation of the.mre file for a core directory using the command

`python -m eleos --make-summary path/to/core/directory`

At the moment, this only works for cores created by Eleos, as it requires the core.pkl file to be present which contains a serialisation of the NemesisCore object used to create the core. This allows very easy access to all the profiles, shapes, attributes etc… In the future, this restriction might be lifted.


## Limitations and future work

Currently, this library only supports Jupiter (although expansion should be fairly easy when providing .ref files), preset k-tables, and is not as flexible generally as using Nemesis directly. In the future, these will be relaxed once the basics of the library have been debugged and tested for both forward and retrieval modes. The advantage of Eleos is the ablilty to create hundreds of consistent cores through a standard interface, without relying on bespoke code to interface with NEMESIS.