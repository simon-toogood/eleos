import os
import subprocess


class NemesisCore:
    def __init__(self, name, directory):
        self.name = name
        self.directory = directory
        os.mkdir(self.directory)
        os.mkdir("tmp")

    def generate_fmerror(self):
        pass

    def generate_xsc(self, num_aerosol_modes, start_wl, end_wl, delta_wl, radius, variance, real_refactive_index, imag_refractive_index):
        """Run Makephase and Normxsc with the given inputs to generate the .xsc file containing aerosol refractive indicies.
        Currently only supports mode 1 (Mie scattering, gammma dist.) and constant refractive index over range"""

        # Generate the makephase.inp file
        with open("tmp/makephase.inp", mode="w+") as file:
            file.write(f"{num_aerosol_modes}\n1\n{start_wl} {end_wl} {delta_wl}\n{self.name}.xsc\ny\n1\n{radius} {variance}\n2\n1\n{real_refactive_index} {imag_refractive_index}")

        # Generate the normxsc.inp file
        with open("tmp/normxsc.inp", mode="w+") as file:
            file.write(f"{self.name}.xsc\n1 1")

        # Run Makephase and Normxsc
        cwd = os.getcwd()
        os.chdir(self.directory)
        os.system("Makephase < ../tmp/makephase.inp")
        os.system("Normxsc < ../tmp/normxsc.inp")
        os.chdir(cwd)

    def generate_cloudf(self, wavelengths, imag_refractive_index, imag_refractive_index_err):
        """Generate the cloudf1.dat file. Contains some header info (???) and then 3 columns:
           wavelength, imag refractive index and the error on it"""
        
        header = f"    0.01    0.1\n    0.01    0.001\n{len(wavelengths)}   -1    !NWAVE,CLEN\n2.73  1.55\n2.73        !V_OD_NORM"
        with open("fcloud1.dat", mode="w+") as file:
            file.write(header + "\n")
            for wl in wavelengths:
                file.write(f"{wl} {imag_refractive_index} {imag_refractive_index_err}\n")


core = NemesisCore("test", "core_1")
