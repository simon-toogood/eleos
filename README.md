# Eleos - A Python interface to NEMESIS

Please note, this is a WIP project and features may be added/removed at any time. It is absolutely *not* guarunteed to be backward compatible until it has matured significantly.

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

```python
from eleos import cores, shapes, profiles


# Create a profile to retrieve ammonia - this is parameterised as model 20 in NEMESIS
nh3 = profiles.GasProfile(gas_name="NH3", 
                          shape=shapes.Shape20(
                          knee_pressure=1.0, 
                          tropopause_pressure=0.1,
                          deep_vmr=2e-4, deep_vmr_error=10e-6,
                          fsh=0.2,           fsh_error=0.1))

# Create a profile to retrieve phosphine - this is a simple scaling of the prior distribution
ph3 = profiles.GasProfile(gas_name="PH3", 
                          shape=shapes.Shape2(
                          scale_factor=2, scale_factor_error=0.1))


# Create a Gaussian-shaped cloud layer at 1.2 bar with set optical properties
deep = profiles.AerosolProfile(label="Deep Cloud",
                               retrieve_optical=False,
                               shape=shapes.Shape47(
                               central_pressure=1.23, central_pressure_error=0.04,
                               pressure_width=0.1,    pressure_width_error=0.1,
                               opacity=1.75,          opacity_error=0.5),
                               radius=3.02,           
                               variance=0.5,          
                               real_n=1.3,
                               imag_n=1e-3)

# Create another Gaussian cloud at 0.6 bar and retrieve the optical properties as well (particle radius, variance and imag. refractive index)
main = profiles.AerosolProfile(label="Main Cloud",
                               retrieve_optical=True,
                               shape=shapes.Shape47(
                               central_pressure=0.6,  central_pressure_error=0.01,
                               pressure_width=0.2,    pressure_width_error=0.2,
                               opacity=1.0,           opacity_error=0.5),
                               radius=3.02,           radius_error=0.1,
                               variance=0.1,          variance_error=0.1,
                               real_n=1.3, 
                               imag_n=1e-3,           imag_n_error=1e-3)

# Create an aerosol layer that represents a uniformally distributed haze between 0.5bar and 0.1bar and retrieve the optical properties
haze = profiles.AerosolProfile(label="Haze",
                               retrieve_optical=True,
                               shape=shapes.Shape37(
                               bottom_pressure=0.5, 
                               top_pressure=0.1,
                               opacity=0.25,       opacity_error=0.1),
                               radius=0.34,        radius_error=0.02,
                               variance=0.3,       variance_error=0.3,
                               real_n=1.3,
                               imag_n=1e-3,        imag_n_error=1e-3)

# Set the directory of the core and clear it of any previous retrievals
cd = "example/"
cores.clear_parent_directory(cd)

# Create the NemesisCore object with the profiles we defined
core = cores.NemesisCore(cd,
                         planet="jupiter",
                         instrument_ktables="NIRSPEC",
                         spx_file="nearnadir.spx",
                         profiles=[ph3, nh3, deep, main, haze],
                         fmerror_factor=5,
                         num_iterations=2,
                         scattering=True,
                         reference_wavelength=4)

# Set some pressure limits - we only have sensitivity between approx 1mbar and 10bar
core.set_pressure_limits(min_pressure=1e-3, max_pressure=10)

# If there is a specific feature at a given wavelength that we want NEMESIS to always fit, we can use the fix_peak method
core.fix_peak(central_wavelength=4.07, width=0.05)

# Create all the files necessary for NEMESIS to run
core.generate_core()

# Run the job on ALICE - this will also run Eleos from the command line to make some summary plots in example/core_1/plots/
cores.generate_alice_job(cd, python_env_name="pythonmain", username="scat2", hours=3)
cores.run_alice_job(cd)
```

## Result Analysis Example (analyse.py)

```python
import matplotlib.pyplot as plt
from eleos import results


# Load the retrieved core as a NemesisResult object
res = results.NemesisResult("example/core_1")

# Get some information about the run
res.print_summary()

# Create two plots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))

# Plot the retrieved spectrum and aerosol optical thickness as a function of pressure 
res.plot_spectrum(ax=ax1)
res.plot_aerosol_profiles(ax=ax2)

# Save the figure
fig.tight_layout()
fig.savefig("example.png", dpi=500)
```


## Command Line

Eleos can be run at the command line in order to quickly generate summary plots and print a human-readable representation of the.mre file for a core directory using the command

`python -m eleos --make-summary path/to/core/directory`

At the moment, this only works for cores created by Eleos, as it requires the core.pkl file to be present which contains a serialisation of the NemesisCore object used to create the core. This allows very easy access to all the profiles, shapes, attributes etc… In the future, this restriction might be lifted.


## Limitations and future work

Currently, this library only supports Jupiter (although expansion should be fairly easy when providing .ref files), preset k-tables, and is obviously not as flexible generally as using Nemesis directly. In the future, these will be relaxed once the basics of the library have been debugged and tested for both forward and retrieval modes. The advantage of Eleos is the ablilty to create hundreds of consistent cores through a standard interface, without relying on bespoke code to interface with NEMESIS.