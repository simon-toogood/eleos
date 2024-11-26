import cores
import profiles
import shapes


nh3_shape = shapes.Shape1(knee_pressure=0.1, 
                          deep_vmr=1e-4, 
                          deep_vmr_error=1e-4, 
                          fsh=0.3, 
                          fsh_error=0.3)
aero_shape = shapes.Shape32(base_pressure=0.8, 
                            base_pressure_error=0.5, 
                            opacity=1, 
                            opacity_error=0.3,
                            fsh=0.4,
                            fsh_error=0.2)


nh3_profile = profiles.GasProfile(gas_name="NH3", 
                                  isotope_id=0, 
                                  shape=nh3_shape)
aero_profile = profiles.AerosolProfile(aerosol_id=1, 
                                       shape=aero_shape)
temp_profile = profiles.TemperatureProfile(filepath="./data/jupiter/tempapr.dat")


core_list = []
for n in range(1, 5):
    core = cores.NemesisCore(parent_directory=f"cores/",
                             spx_file="/home/s/scat2/JWST/2022_JupSouthPole/zonal_spectra/sparse_55.0degS.spx",
                             ref_file="data/jupiter/jupiter.ref",
                             profiles=[temp_profile, nh3_profile, aero_profile],
                             fmerror_factor=n)
    core_list.append(core)


cores.generate_alice_job(cores=core_list, username="scat2")
cores.run_alice_job("cores/")
