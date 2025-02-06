from eleos import shapes, profiles, cores


# ------ Gas Profile definitions ------ #

ph3 = profiles.GasProfile(gas_name="PH3", 
                          shape=shapes.Shape4(
                                knee_pressure=0.8, knee_pressure_error=0.8,
                                deep_vmr=0.0005,   deep_vmr_error=0.0005,
                                fsh=0.2,           fsh_error=0.2))

nh3 = profiles.GasProfile(gas_name="NH3", 
                          shape=shapes.Shape2(
                                scale_factor=1.0, scale_factor_error=0.1))


# ------ Tropospheric cloud layer ------ #

aero1 = profiles.AerosolProfile(shape=shapes.Shape32(
                                    base_pressure=0.8, base_pressure_error=0.5,
                                    opacity=5.0,       opacity_error=2.0,
                                    fsh=0.7,           fsh_error=0.9))
                

# ------ Stratospheric haze ------ #

aero2 = profiles.AerosolProfile(shape=shapes.Shape32(
                                    base_pressure=0.1, base_pressure_error=0.5,
                                    opacity=2.0,       opacity_error=2.0,
                                    fsh=0.9,           fsh_error=0.9))
                
n2 = profiles.ImagRefractiveIndexProfile(shape=shapes.Shape444(
                                             radius=3,                   radius_error=1,
                                             variance=0.5,               variance_error=0.1,
                                             refractive_index=1.2+2e-3j, refractive_index_error=1e-3))


# ------ Core creation ------ #

cd = "example/"
cores.clear_parent_directory(cd)

core = cores.NemesisCore(parent_directory=cd,
                         spx_file="your_spx_file_here.spx",
                         profiles=[ph3, nh3],
                         fmerror_factor=10)

core.add_aerosol_mode(aero1, radius=2, varaince=0.1, refractive_index=1.3+1e-3j)
core.add_aerosol_mode(aero2, n2)

core.generate_core()

cores.generate_alice_job(cd, username="none", hours=2)
cores.run_alice_job(cd)
