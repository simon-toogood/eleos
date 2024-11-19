# Eleos - A Python interface to NEMESIS

This is a WIP project and features may be added/removed at any time. Full documentation will be added in due course.

## Core Generation Example

This code generates 4 cores, each with a different forward modelling error factor. It retrieves the temeprature profile using a prior from a file (`tempapr.dat`), the ammonia profile represented as a knee pressure (model 1 in NEMESIS), and an aerosol layer represented as model 32. It then generates a submission script to run NEMESIS using those cores on ALICE.

```python

from eleos import shapes, profiles, cores

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
```

## Result Analysis Example

This code takes the result of running NEMESIS on the output of the above example and plots the retrieved spectrum and temperature profile, then saves it to a file. 

```python
res = NemesisResult("cores/core_1/")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))
res.plot_spectrum(ax=ax1)
res.plot_temperature(ax=ax2)
plt.tight_layout()
fig.savefig("nosync/temp.png", dpi=500)
```


