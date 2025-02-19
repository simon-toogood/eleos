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