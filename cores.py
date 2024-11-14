import os
import shutil
import pandas as pd
from jwstools import spx

import constants

# Files:
#   Thermal:
#     * aerosol.ref (aerosol profiles as function of height)
#     * parah2.ref  (don't know, just copy it)
#     * fmerror.dat (additional error for each wavelength)
#     * nemesis.fla (various flags)
#     * nemesis.cia (pointer to collision-induced opacity file)
#     * nemesis.kls (list of k-tables)
#     - nemesis.ref (reference heights, pressures, temps and composition)
#     * nemesis.xsc (aerosol cross-sections)
#     * nemesis.set (scattering angles, layer info)
#     * nemesis.apr (a priori file containing vars to retrieve)
#     * nemesis.spx (spectrum)
#   **  .abo, .nam
#   

class NemesisCore:
    def __init__(self, directory, spx_file, ref_file, profiles=list(), planet="jupiter", scattering=False):
        self.directory = directory
        self.spx_file = spx_file
        self.ref_file = ref_file
        self.profiles = profiles
        self.planet = planet
        self.scattering = scattering
        if os.path.exists(self.directory):
            shutil.rmtree(self.directory)
        os.mkdir(self.directory)
        os.mkdir(self.directory + "tmp")

    def __str__(self):
        return f"<NemesisCore: {self.directory}>"

    def parse_ref_file(self):
        """Read in the .ref file provided and return a DataFrame
        
        Args:
            None
            
        Returns:
            pandas.DataFrame with columns for height (in km), pressure (in atm), temperature (in K), and the VMRs for all specified gases"""
        with open(self.ref_file) as file:
            for i, line in enumerate(file):
                if "height" in line:
                    skip_to = i
                if i == 2:
                    tokens = line.split()
                    n_gas = int(tokens[-1])
        df = pd.read_table(self.ref_file, skiprows=skip_to+1, sep="\s+", names=["height", "pressure", "temp"] + [f"VMR gas {n+1}" for n in range(n_gas)])
        return df

    def copy_template_files(self):
        template_core = "./data/template_core/"
        shutil.copy(template_core+"nemesis.cia", self.directory)
        shutil.copy(template_core+"nemesis.abo", self.directory)
        shutil.copy(template_core+"nemesis.nam", self.directory)
        shutil.copy(template_core+"parah2.ref" , self.directory)

    def generate_apr(self):
        """Generate the nemesis.apr file from the profile list

        Args:
            None
        
        Returns:
            None
            
        Creates:
            nemesis.apr"""
        out = f"*******Apriori File*******\n           {len(self.profiles)}\n"
        for profile in self.profiles:
            out += profile.generate_apr_data() + "\n"
        with open(self.directory + "nemesis.apr", mode="w") as file:
            file.write(out)

    def generate_set(self, num_layers=120, bottom_layer_height=-80):
        """Generate the settings file for NEMESIS
        
        Args:
            num_layers: The number of atmospheric layers
            bottom_layer_height: The altitude of the bottom layer of the atmosphere
            
        Returns:
            None
            
        Creates:
            nemesis.set"""
        
        with open("./data/template_core/nemesis.set", mode="r") as file:
            out = file.read()
        out = out.replace("<DISTANCE>", f"{constants.DISTANCES[self.planet]:.3f}")
        out = out.replace("<SUNLIGHT>", f"{int(self.scattering)}")
        out = out.replace("<BOUNDARY>", f"{int(self.scattering)}")
        out = out.replace("<BASE>", f"{bottom_layer_height:.2f}")
        out = out.replace("<NUM_LAYERS>", f"{num_layers}")
        with open(self.directory+"nemesis.set", mode="w+") as file:
            file.write(out)

    def generate_flags(self, inormal=0, iray=1, ih2o=0, ich4=0, io3=0, inh3=0, iptf=0, imie=0, iuv=0):
        """Generate the flags file. As a general rule, leave all this at the defaults. The descriptions are copied directly fron the NEMESIS manual.
        Args:
            inormal: 0 or 1 depending on whether the ortho/para-H2 ratio is in equilibrium (0) or normal 3:1 (1).
            iray: Sets the Rayleigh optical depth calculation:
                    0 = Rayleigh scattering optical depth is not included.
                    1 = Rayleigh optical depths suitable for gas giant atmosphere
                    2 = Rayleigh optical depths suitable for CO2 dominated atmosphere
                    3 = Rayleigh optical depths suitable for a N2-O2 atmosphere.
                    4 = New Rayleigh optical depths suitable for a gas giant atmosphere.
                    5 = New Rayleigh optical depths suitable for a gas giant atmosphere PLUS Raman scattering.
                    6 = New Rayleigh optical depths suitable for a gas giant atmosphere PLUS Sromovsky Raman scattering, PLUS Sromovsky Polarisation Correction.
            ih2o: Turns additional H2O continuum off (0) or on (1)
            ich4: Turns additional CH4 continuum off (0) or on (1)
            io3:  Turns additional O3 UV continuum off (0) or on (1)
            inh3: Turns additional NH3 continuum off (0) or on (1)
            iptf: Used in only a few routines to switch between normal partition function calculation (0) or the high-temperature partition function for CH4 for Hot Jupiters.
            imie: Only relevant for scattering calculations. If set to 0, the phase function is computed from the associated Henyey-Greenstein hgphase*.dat files. However, if set to 1, the phase function is computed from the MieTheory calculated PHASEN.DAT file.
        
        Returns:
            None
            
        Creates:
            nemesis.fla"""
        
        with open(self.directory+"nemesis.fla", mode="w+") as file:
            file.write('       {aa}	! Inormal (0=eqm, 1=normal)\n       {bb}	! Iray (0=off, 1=on)\n       {cc}	! IH2O (1 = turn extra continuum on)\n       {dd}	! ICH4 (1 = turn extra continuum on)\n       {ee}	! IO3 (1 = turn extra continuum on)\n       {ff}	! INH3 (1 = turn extra continuum on)\n       {gg}	! Iptf (0=default, 1=CH4 High-T)\n       {hh}	! IMie\n       {ii}	! UV Cross-sections\n       0	! INLTE (0=LTE)'.format(aa=inormal, bb=iray, cc=ih2o, dd=ich4, ee=io3, ff=inh3, gg=iptf, hh=imie, ii=iuv))
            
    def generate_aerosol_ref(self, num_aerosol_modes):
        """Generate a number of aerosol reference profile that is 0 at all altitudes.
        
        Args:
            num_aerosol_modes: Number of aerosol modes to generate
            
        Returns:
            None
            
        Creates:
            aerosol.ref"""
        heights = self.parse_ref_file().height
        print(heights)
        with open(self.directory+"aerosol.ref", mode="w+") as file:
            file.write(f"# Generated by Eleos\n")
            file.write(f"{len(heights)} {num_aerosol_modes}\n")
            for h in heights:
                file.write(f"{h:>12.5f} ")
                for i in range(num_aerosol_modes):
                    file.write("0.000000E+00 ")
                file.write("\n")

    def generate_kls(self, instrument="NIRSPEC"):
        """Copy the ktables from the template core for the given instrument.
        Args:
            instrument: Either 'NIRSPEC' or 'MIRI'
            
        Returns:
            None
            
        Creates:
            nemesis.kls
            
        TODO:
            Add way to include/exclude different elements"""
        if instrument == "NIRSPEC":
            shutil.copy("./data/template_core/nirspec.kls", self.directory+"nemesis.kls")
        elif instrument == "MIRI":
            shutil.copy("./data/template_core/miri.kls", self.directory+"nemesis.kls")

    def generate_fmerror(self, factor=3):
        """For each wavelength in the spx file, multiply the error by factor. only supports 1 spx geom"""
        spx_data = spx.read(self.spx_file).geometries[0]
        num_entries = len(spx_data.wavelengths)
        if num_entries > 2048:
            raise IndexError(f"spx file has too many wavelengths! ({num_entries}/2048)")
        with open(self.directory + "fmerror.dat", mode="w+") as file:
            file.write(f"{num_entries}\n")
            for wl, err in zip(spx_data.wavelengths, spx_data.error):
                file.write(f"{wl:.6e}  {err*factor:.6e}\n")

    def generate_xsc(self, num_aerosol_modes, start_wl, end_wl, delta_wl, radius, variance, real_refactive_index, imag_refractive_index):
        """Run Makephase and Normxsc with the given inputs to generate the .xsc file containing aerosol refractive indicies.
        Currently only supports mode 1 (Mie scattering, gammma dist.) and constant refractive index over range.
        
        Args:
            num_aerosol_modes: The number of aerosol modes to create. CURRENTLY ONLY SUPPORTS 1
            start_wl: The shortest wavelength (in microns) to compute the refractive index for
            end_wl: The longest wavelength (in microns) to computethe refractive index for
            delta_wl: Wavelength resolution (in microns)
            radius: The mean radius of the particles (in microns)
            variance: The varaince of the particle radius (in microns)
            real_refractive_index: The real component of the refractive index of the particles
            imag_refractive_index: The imaginary component of frefractive index for the particles
            
        Returns:
            None
        
        Creates:
            tmp/makephase.inp
            tmp/normxsc.inp
            Whatever files Makephase makes
            
        TODO:
            Add suppport for multiple aerosol modes
            Add support for other aerosol scattering properties
            Swap real and imap parts for a standard Python complex number?
            Read start_wl and end_wl from spx file"""

        # Generate the makephase.inp file
        with open(self.directory+"tmp/makephase.inp", mode="w+") as file:
            file.write(f"{num_aerosol_modes}\n1\n{start_wl} {end_wl} {delta_wl}\nnemesis.xsc\ny\n1\n{radius} {variance}\n2\n1\n{real_refactive_index} {imag_refractive_index}")

        # Generate the normxsc.inp file
        with open(self.directory+"tmp/normxsc.inp", mode="w+") as file:
            file.write(f"nemesis.xsc\n1 1")

        # Run Makephase and Normxsc
        cwd = os.getcwd()
        os.chdir(self.directory)
        os.system("Makephase < tmp/makephase.inp")
        os.system("Normxsc < tmp/normxsc.inp")
        os.chdir(cwd)

    def _generate_default_core(self):
        self.copy_template_files()
        self.generate_apr()
        self.generate_set()
        self.generate_flags()
        self.generate_aerosol_ref(1)
        self.generate_kls()
        self.generate_fmerror()
        self.generate_xsc(1, 1.8, 5.3, 0.1, 1, 0.1, 1.3, 1e-3)


    def generate_cloudf(self, wavelengths, imag_refractive_index, imag_refractive_index_err):
        """Generate the cloudf1.dat file. Contains some header info (???) and then 3 columns:
        wavelength, imag refractive index and the error on it.
           
        Args:
            wavelengths: List of wavelengths 
            imag_refractive_index: The imaginary part of refractive index for the clouds
            imag_refractive_index_err: The estimated uncertainty on the imaginary refractive index

        Returns:
            None

        Creates:
            fcloud1.dat
        """
        
        header = f"    0.01    0.1\n    0.01    0.001\n{len(wavelengths)}   -1    !NWAVE,CLEN\n2.73  1.55\n2.73        !V_OD_NORM"
        with open("fcloud1.dat", mode="w+") as file:
            file.write(header + "\n")
            for wl in wavelengths:
                file.write(f"{wl} {imag_refractive_index} {imag_refractive_index_err}\n")


