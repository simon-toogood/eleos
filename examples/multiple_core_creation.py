from eleos import cores, shapes, profiles
import numpy as np


# Define a parameterised phosphine profile
ph3 = profiles.GasProfile(gas_name="PH3", 
                          shape=shapes.Shape20(
                          knee_pressure=1.0, 
                          tropopause_pressure=0.1,
                          deep_vmr=1.86e-6,          deep_vmr_error=0.2e-6,
                          fsh=0.2,                   fsh_error=0.1))
        

# And a thick cloud layer
main = profiles.AerosolProfile(label="Main Cloud",
                               retrieve_optical=True,
                               shape=shapes.Shape48(
                               base_pressure=2,      base_pressure_error=0.5,
                               top_pressure=0.1,     top_pressure_error=0.2,
                               opacity=1,            opacity_error=0.5,
                               fsh=0.3,              fsh_error=0.05),
                               radius=3.0,           radius_error=0.1,
                               variance=0.5,         variance_error=0.2,
                               real_n=1.3,
                               imag_n=1e-3,          imag_n_error=1e-3)

# And a thin haze layer
haze = profiles.AerosolProfile(label="Haze",
                               retrieve_optical=True,
                               shape=shapes.Shape37(
                               bottom_pressure=0.1, 
                               top_pressure=0.01,
                               opacity=0.25,         opacity_error=0.1),
                               radius=0.34,          radius_error=0.02,
                               variance=0.3,         variance_error=0.3,
                               real_n=1.3,  
                               imag_n=1e-3,          imag_n_error=1e-3)

# Clear the target directory
cd = "multiple_example/"
cores.clear_parent_directory(cd)

# Loop over the desired haze opacities
for haze_opacity in np.arange(1e-20, 2.25, 0.25):

    # Set the haze opacity
    haze.opacity = haze_opacity

    # Create and generate the core as normal
    core = cores.NemesisCore(cd,
                            spx_file="data/nearnadir.spx",
                            profiles=[ph3, main, haze],
                            fmerror_factor=4,
                            reference_wavelength=4)

    core.set_pressure_limits(min_pressure=1e-3, max_pressure=10)

    core.generate_core()

# Run the cores an an array job on ALICE
cores.generate_alice_job(cd, python_env_name="pythonmain", username="scat2", hours=24)
cores.run_alice_job(cd)