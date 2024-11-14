import cores
import profiles
import shapes


nh3_shape = shapes.Shape1(knee_pressure=1, 
                          deep_vmr=2, 
                          deep_vmr_error=3, 
                          fsh=4, 
                          fsh_error=5)
nh3_profile = profiles.GasProfile(gas_name="NH3", 
                                  isotope_id=0, 
                                  shape=nh3_shape)

aero_shape = shapes.Shape32(base_pressure=1e2, 
                            base_pressure_error=1e1, 
                            opacity=1, 
                            opacity_error=0.3,
                            fsh=0.4,
                            fsh_error=0.2)
aero_profile = profiles.AerosolProfile(aerosol_id=1, shape=aero_shape)

temp_profile = profiles.TemperatureProfile(filepath="tempapr.dat")


core = cores.NemesisCore(directory=f"core_test/",
                        spx_file="/home/s/scat2/NEMESIS/2022_JupSouthPole/jwst_mrs.spx",
                        ref_file="data/jupiter.ref",
                        profiles=[temp_profile, nh3_profile, aero_profile])

core._generate_default_core()

# import results
# res = results.NemesisResult("nemesis", "/home/s/scat2/NEMESIS/2022_JupSouthPole/core_0/")
# fig, ax = res.plot_spectrum()
# fig.savefig("spectrum.png", dpi=500)