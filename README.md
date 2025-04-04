# Eleos - A Python interface to NEMESIS

Please note, this is a WIP project and features may be added/removed at any time. It is absolutely *not* guarunteed to be backward compatible until it has matured significantly.

## Documentation

Documentation can be found on ReadTheDocs [here](https://eleos.readthedocs.io/en/latest/). This is still WIP, full documentation will be added in due course.

## Installation

Eleos is available on PyPI so it can be installed as any other Python package. On Unix-like systems:

`pip install nemesis_eleos`

And on Windows:

`py -m pip install nemesis_eleos`


While `eleos` only creates files for, and reads files created by, NEMESIS, it is mandatory to have NEMESIS installed for this package to work. This is primarily due to requiring the utilities `Makephase` and `Normxsc` on the PATH. See the NEMESIS GitHub page for full instructions on [how to download the software](https://github.com/nemesiscode/radtrancode/blob/master/README.md) and [how to compile it](https://github.com/nemesiscode/radtrancode/blob/master/AACOMPILE.txt).


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

**More examples can be found in the** `examples/` **directory.**


The `eleos.cores` module is designed to be used to generate cores ready for NEMESIS to run.  This code generates a core with NH3 and PH3 profiles and three aerosol layers for a JWST NIRSPEC observation of Jupiter. It uses the built-in .ref files, ktables and other miscellaneous NEMESIS input files (available in `eleos/data/`). Then it generates a submission script for ALICE and submits the job to the scheduler. After the core has successfully run, there will be a selection of summary plots in the `parent/core/plots` directory.

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
                         spx_file="example.spx",
                         profiles=[ph3, nh3, deep, main, haze],
                         fmerror_factor=5,
                         num_iterations=20,
                         scattering=True,
                         reference_wavelength=4,
                         min_pressure=1e-3, 
                         max_pressure=10)

# If there is a specific feature at a given wavelength that we want NEMESIS to always fit, we can use the fix_peak method
core.fix_peak(central_wavelength=4.07, width=0.05)

# Create all the files necessary for NEMESIS to run
core.generate_core()

# Run the job on ALICE - this will also run Eleos from the command line to make some summary plots in example/core_1/plots/
cores.generate_alice_job(cd, python_env_name="pythonmain", username="scat2", hours=3)
cores.run_alice_job(cd)
```

## Result Analysis Example (analyse.py)

Once NEMESIS has been run on the core, the `eleos.results` module is used to analyse the output. This script takes the output of the previous and generates a set of summary tables containing the prior and retrieved values for each parameter, a plot of the retrieved spectrum, and a plot of the aerosol densities.


```python
import matplotlib.pyplot as plt
from eleos import results


# Load the retrieved core as a NemesisResult object
res = results.NemesisResult("example/core_1")

# Get some information about the run
res.print_summary()

# Create two plots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))

# Plot the retrieved spectrum, and aerosol density as a function of pressure 
res.plot_spectrum(ax=ax1)
res.plot_aerosol_profiles(ax=ax2)

# Save the figure
fig.tight_layout()
fig.savefig("example.png", dpi=500)
```


## Sensitivity Analysis

To determine the effect of each parameter on the spectrum, a sensitivity analysis can be run for a given core. This varies each parameter in each profile by a series of factors (by default 80%, 90%, 95%, 105%, 110%, and 120%) and looks at the relative change in the spectrum compared to not changing anything. To generate the analysis cores we can use the `load_from_previous` function in `eleos.cores` to load the core we want to analyse  and `create_sensitivity_analysis` to generate the forward model cores with the tweaked parameters.


```python
from eleos import cores


# Define and clear the new parent directory
cd = "sensitivity/"
# cores.clear_parent_directory(cd)

# Load the previously retrieved core
core = cores.load_from_previous("example/core_1/", cd)

# Create a sensitivity analysis and run it
cores.create_sensitivity_analysis(core)
cores.generate_alice_job(cd, python_env_name="pythonmain", username="scat2", hours=2, memory=1)
cores.run_alice_job(cd)
```


After all the cores have run, we can use the `eleos.results` module again to analyse the results.


```python
from eleos import results

cd = "sensitivity/"

sens = results.SensitivityAnalysis(cd)
sens.make_parameters_plot()
sens.savefig("sensitivity.png", dpi=400)
```


or, alternatively, we can run the equivalent script from the command-line using

`python -m eleos --make-sensitivity-summary sensitivity/`

## Command Line

Eleos can be run at the command line in order to quickly generate summary plots and print a human-readable representation of the.mre file for a core directory using the command

`python -m eleos --make-summary path/to/core/directory`

At the moment, this only works for cores created by Eleos, as it requires the core.pkl file to be present which contains a serialisation of the NemesisCore object used to create the core. This allows very easy access to all the profiles, shapes, attributes etc… In the future, this restriction might be lifted.

## Limitations and future work

Eleos is not and will never be as flexible as NEMESIS. It is intended to provide a framework for automating the boilerplate code used in generating hundreds of consistent cores through a standard interface, without relying on bespoke code to interface with NEMESIS.

Currently, this library only supports Jupiter (although expansion should be fairly easy when providing .ref files), preset k-tables. In the future, these will be relaxed once the basics of the library have been debugged and tested for both forward and retrieval modes.

A big limitation is the lack of TemperatureProfiles. While these have been partially implemented, for my work in the NIR range we keep the temeprature fixed so adding this is not a priority for me. In the future it will be added though!