import os
import shutil
import pandas as pd
import pickle
import time
from pathlib import Path
import numpy as np
import sys

import warnings
warnings.formatwarning = lambda msg, *_: f"Warning: {msg}\n"

from . import constants
from . import spx
from . import parsers
from . import utils
from . import profiles as profiles_ # to avoid namespace shadowing by NemesisCore.profiles


CORE_ID_COUNTER = 0


class NemesisCore:
    """The class that constructs a core directory and all the required files for NEMESIS to run."""
    def __init__(self, 
                 parent_directory,
                 spx_file, 
                 ref_file=None, 
                 profiles=list(), 
                 planet="jupiter", 
                 scattering=True, 
                 forward=False,
                 prompt_if_exists=True, 
                 num_iterations=30,
                 num_layers=120, 
                 bottom_layer_height=None, 
                 instrument_ktables="NIRSPEC", 
                 fmerror_factor=0,
                 fmerror_pct=None,
                 fmerror_value=None,
                 cloud_cover=1.0,
                 reference_wavelength=None):
        
        """Create a NEMESIS core directory with a given set of profiles to retrieve
        
        Args:
            parent_directory (str):       The directory in which to create the core folder 
            spx_file (str):               Path to the spectrum file to fit
            ref_file (str):               Path to the ref file to use (if left blank it will use the default for that planet)
            profiles List[str]:           List of profiles.Profile objects to retrieve
            planet (str):                 Name of the planet being observed. Must be one of 'jupiter', 'saturn', 'uranus', 'neptune' or 'titan'.
            scattering (bool):            Whether to run a scattering retrieval or not
            forward (bool):               Whether of not to run a forward model (ie. set number of iterations = 0)
            prompt_if_exists (bool):      If the core directory already exists, then ask before continuing
            num_iterations (int):         Number of iterations to run in the retrieval (if forward is set this has no effect)
            num_layers (int):             The number of atmospheric layers to simulate
            bottom_layer_height (int):    The height in km of the bottom of the atmosphere (by defauylt use the lowest height in the .ref file)
            instrument_ktables (str):     Either 'NIRSPEC' or 'MIRI'; determines which set of ktables to use
            fmerror_factor (float):       The factor by which to multiply the error on the spectrum (see also, fmerror_pct and fmerror_value)
            fmerror_pct (float):          If given, instead of using fmerror_factor or fmerror_value, use a flat percentage of the brightness (eg. 0.1 = 10%) (see also, fmerror_factor and fmerror_value)
            fmerror_value (float):        If given, instead of using fmerror_factor or fmerror_pct, use a flat value in W/cm2/sr/um (see also, fmerror_factor and fmerror_pct)
            cloud_cover (bool):           If scattering mode is on, then this is the fractional cloud cover between 0 and 1 (usually doesn't need to be changed)
            reference_wavelength (float): If scattering mode is on, then normalise the cross-sections at this wavelength.
        """
        # Increment the global core counter and store local version
        global CORE_ID_COUNTER
        CORE_ID_COUNTER += 1
        self.id_ = CORE_ID_COUNTER

        # Assign attributes passed in
        self.spx_file = Path(spx_file).resolve()
        self.profiles = profiles
        self.planet = planet.lower()
        self.scattering = scattering
        self.forward = forward
        self.num_iterations = num_iterations
        self.num_layers = num_layers
        self.instrument_ktables = instrument_ktables
        self.fmerror_factor = fmerror_factor
        self.fmerror_pct = fmerror_pct
        self.fmerror_value = fmerror_value
        self.cloud_cover = cloud_cover
        self.reference_wavelength = reference_wavelength

        # Set the directories of the parent folder and own core
        self.parent_directory = Path(parent_directory).resolve()
        self.directory = self.parent_directory / f"core_{self.id_}"

        # Create the directory tree if it doesn't already exist and clear it if it does
        os.makedirs(self.parent_directory, exist_ok=True)
        if os.path.exists(self.directory) and prompt_if_exists:
            if os.path.exists(self.directory / "nemesis.mre"):
                msg = "There is already a core that has already been run in"
            elif os.path.exists(self.directory / "nemesis.ref"):
                msg = "There is already a core that has not been run yet in"
            else:
                msg = "There is already data in"

            x = input(f"{msg} {self.directory} - erase and continue? Y/N")
            if x.lower() == "n":
                return
            shutil.rmtree(self.directory)
        os.makedirs(self.directory)
        os.mkdir(self.directory / "plots")

        # Check if we are perfoming a scattering run with more than 39 layers (NEMESIS hates this...)
        if num_layers > 39 and self.scattering:
            warnings.warn(f"Too many atmospheric layers specified for a scattering run ({num_layers} vs. 39). Automatically reducing to 39")
            self.num_layers = 39

        # Set ref file if not specified:
        if ref_file is None:
            self.ref_file = constants.PATH / f"data/{planet}/{planet}.ref"
            warnings.warn(f"No ref file specified. Using the default in {self.ref_file}")
        else:
            self.ref_file = ref_file
        
        # Copy in the boilerplate files
        self._copy_input_files()
        self._copy_template_files()

        # Parse the ref file:
        self.ref = parsers.NemesisRef(self.directory / "nemesis.ref")

        # Set bottom layer height if not defined
        if bottom_layer_height is None:
            self.bottom_layer_height = self.ref.data.iloc[0].height
        else:
            self.bottom_layer_height = bottom_layer_height

        # Set layer type (by default this is 1 - equal log pressure grid)
        self.layer_type = 1
        self.min_pressure = self.ref.data.iloc[-1].pressure
        self.max_pressure = self.ref.data.iloc[0].pressure

        # Set number of aerosol modes (incremented by add_aerosol_mode)
        self.num_aerosol_modes = 0

        # Add a list to hold any FixedPeak classes the user wants
        self.fixed_peaks = []

        # Add a reference to self in each Profile and check that there are no Aerosol profiles
        for profile in self.profiles:
            if isinstance(profile, (profiles_.AerosolProfile, profiles_.ImagRefractiveIndexProfile)):
                raise ValueError("Aerosol-related profile specified in constructor - please use core.add_aerosol_profile(...) instead!")
            profile.core = self

        # If in forward mode, set the number of iterations to 0
        if self.forward:
            self.num_iterations = 0

        # Raise an error if trying to use features not implemented yet
        if planet != "jupiter":
            raise Exception("Eleos does not support planets other than Jupiter yet! ")

    def __str__(self):
        return f"<NemesisCore: {self.directory}>"

    def _save_core(self):
        with open(self.directory / "core.pkl", mode="wb") as file:
            pickle.dump(self, file)

    def _add_profile(self, profile):
        profile.core = self
        self.profiles.append(profile)

    def _copy_input_files(self):
        """Copy the given .spx and .ref file into the core, as well as the eleos generation script itself
        
        Args:
            None
            
        Returns:
            None
            
        Creates:
            nemesis.spx
            nemesis.ref
            eleos_generation.py"""

        shutil.copy(self.ref_file, self.directory / "nemesis.ref")
        shutil.copy(self.spx_file, self.directory / "nemesis.spx")
        shutil.copy(sys.argv[0], self.directory / "eleos_generation.py")

    def _copy_template_files(self):
        shutil.copy(constants.PATH / "data/statics/nemesis.cia", self.directory)
        shutil.copy(constants.PATH / "data/statics/nemesis.abo", self.directory)
        shutil.copy(constants.PATH / "data/statics/nemesis.nam", self.directory)
        shutil.copy(constants.PATH / "data/statics/nemesis.sol", self.directory)
        shutil.copy(constants.PATH / "data/statics/makephase.inp", self.directory)

    def _generate_inp(self):
        """Generate the nemesis input file

        Args:
            num_iterations: Number of iterations to run (has no effect if NemesisCore.forward is True)
        
        Returns:
            None
            
        Creates:
            nemesis.inp"""
        with open(constants.PATH / "data/statics/nemesis.inp") as file:
            out = file.read()
            out = out.replace("<SCATTERING>", str(int(self.scattering)))
            out = out.replace("<N_ITERATIONS>", str(self.num_iterations))
        with open(self.directory / "nemesis.inp", mode="w+") as file:
            file.write(out)

    def _generate_apr(self):
        """Generate the nemesis.apr file from the profile list

        Args:
            None
        
        Returns:
            None
            
        Creates:
            nemesis.apr"""
        out = f"*******Apriori File*******\n           {len(self.profiles)}\n"
        for profile in self.profiles:
            profile.shape.create_required_files(self.directory)
            out += profile.generate_apr_data() + "\n"
        with open(self.directory / "nemesis.apr", mode="w") as file:
            file.write(out)

    def _generate_set(self):
        """Generate the settings file for NEMESIS
        
        Args:
            None

        Returns:
            None
            
        Creates:
            nemesis.set"""

        with open(constants.PATH / "data/statics/nemesis.set", mode="r") as file:
            out = file.read()
        out = out.replace("<DISTANCE>", f"{constants.DISTANCES[self.planet]:.3f}")
        out = out.replace("<SUNLIGHT>", f"{int(self.scattering)}")
        out = out.replace("<BOUNDARY>", f"{int(self.scattering)}")
        out = out.replace("<N_LAYERS>", f"{int(self.num_layers)}")
        out = out.replace("<BASE>", f"{self.bottom_layer_height:.2f}")
        out = out.replace("<LAYER_TYPE>", str(self.layer_type))
        with open(self.directory / "nemesis.set", mode="w+") as file:
            file.write(out)

    def _generate_flags(self, inormal=0, iray=1, ih2o=0, ich4=0, io3=0, inh3=0, iptf=0, imie=0, iuv=0):
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
        
        with open(self.directory / "nemesis.fla", mode="w+") as file:
            file.write('       {aa}	! Inormal (0=eqm, 1=normal)\n       {bb}	! Iray (0=off, 1=on)\n       {cc}	! IH2O (1 = turn extra continuum on)\n       {dd}	! ICH4 (1 = turn extra continuum on)\n       {ee}	! IO3 (1 = turn extra continuum on)\n       {ff}	! INH3 (1 = turn extra continuum on)\n       {gg}	! Iptf (0=default, 1=CH4 High-T)\n       {hh}	! IMie\n       {ii}	! UV Cross-sections\n       0	! INLTE (0=LTE)'.format(aa=inormal, bb=iray, cc=ih2o, dd=ich4, ee=io3, ff=inh3, gg=iptf, hh=imie, ii=iuv))

    def _generate_aerosol_ref(self):
        """Generate a number of aerosol reference profile that is 0 at all altitudes.
        
        Args:
            None

        Returns:
            None
            
        Creates:
            aerosol.ref"""
        
        # If there are no clouds defined return immediately
        if self.num_aerosol_modes == 0:
            return
        
        heights = self.ref.data.height
        with open(self.directory / "aerosol.ref", mode="w+") as file:
            file.write(f"# Generated by Eleos\n")
            file.write(f"{len(heights)} {self.num_aerosol_modes}\n")
            for h in heights:
                file.write(f"{h:>12.5f} ")
                for i in range(self.num_aerosol_modes):
                    file.write("0.000000E+00 ")
                file.write("\n")

    def _generate_kls(self):
        """Copy the ktables from the template core for the given instrument.
        Args:
            instrument: Either 'NIRSPEC' or 'MIRI'
            
        Returns:
            None
            
        Creates:
            nemesis.kls
            
        TODO:
            Add way to include/exclude different elements
            Add option to use ktables on ALICE rather than prepackaged"""
        if self.instrument_ktables == "NIRSPEC":
            shutil.copy(constants.PATH / "data/jupiter/nirspec.kls", self.directory / "nemesis.kls")
        elif self.instrument_ktables == "MIRI":
            shutil.copy(constants.PATH / "data/jupiter/miri.kls", self.directory / "nemesis.kls")

    def _generate_fmerror(self):
        """For each wavelength in the spx file, adjust the error by a factor (either fmerror_factor, _pct or _value) and write to the fmerror file. 
        Currently only supports 1 spx geometry
        Args:
            None
            
        Returns:
            None
        
        Creates:
            fmerror.dat """
        spx_data = spx.read(self.spx_file).geometries[0]
        num_entries = len(spx_data.wavelengths)

        if num_entries > 2048:
            raise IndexError(f"spx file has too many wavelengths! ({num_entries}/2048)")
        
        with open(self.directory / "fmerror.dat", mode="w+") as file:
            # Header with number of lines 
            file.write(f"{num_entries+2}\n")

            # Catch any cases outside the wavelength range (lower)
            file.write(f"{0.1:.6e}  {1e-8:.6e}\n")

            # Scale the spx error by either fmerror_factor, fmerror_pct or fmerror_value
            for wl, err, spc in zip(spx_data.wavelengths, spx_data.error, spx_data.spectrum):
                # Check to see if this peak is fixed
                flag = False
                for fxp in self.fixed_peaks:
                    if fxp.isin(wl):
                        file.write(f"{wl:.6e}  {fxp.error:.6e}\n")
                        flag = True
                # Otherwise, compute the new error
                if self.fmerror_factor is not None and not flag:
                    file.write(f"{wl:.6e}  {err*self.fmerror_factor:.6e}\n")
                elif self.fmerror_pct is not None and not flag:
                    file.write(f"{wl:.6e}  {spc*self.fmerror_pct:.6e}\n")
                elif self.fmerror_value is not None and not flag:
                    file.write(f"{wl:.6e}  {self.fmerror_value:.6e}\n")

            # Catch any cases outside the wavelength range (upper)
            file.write(f"{100:.6e}  {1e-8:.6e}\n")

    def _generate_xsc(self):
        """Run Makephase and Normxsc with the given inputs to generate the .xsc and phase files.
        
        Args:
            None
            
        Returns:
            None
        
        Creates:
            makephase.inp
            normxsc.inp
            nemesis.xsc
            PHASE{N}.DAT
            hgphase{n}.dat"""
        
        # Check to see if any aerosols have been added - if not then return immediately
        if self.num_aerosol_modes == 0:
            return
        
        # Read in wavelengths bounds from spx file
        wls = spx.read(self.spx_file).geometries[0].wavelengths
        start_wl = min(wls) - 0.1
        end_wl = max(wls) + 0.1
        if self.reference_wavelength is not None:
            ref_wl_idx, _ = utils.find_nearest(wls, self.reference_wavelength)
        else:
            ref_wl_idx = 0
            self.reference_wavelength = start_wl
            warnings.warn(f"No reference wavelength for aerosol cross-sections specified - using the shortest wavelength ({start_wl:4f}um)")

        # Replace the number of aerosol modes and start/end/delta wavelengths
        with open(self.directory / "makephase.inp", mode="r+") as file:
            lines = file.read().split("\n")
            lines[0] = str(self.num_aerosol_modes)
            lines[2] = f"{start_wl} {end_wl} 0.1"
            file.seek(0)
            file.write("\n".join(lines))
            file.truncate()

        # Generate the normxsc.inp file
        with open(self.directory / "normxsc.inp", mode="w+") as file:
            file.write(f"nemesis.xsc\n{ref_wl_idx + 1} 1")

        # Run Makephase and Normxsc
        cwd = os.getcwd()
        os.chdir(self.directory)
        os.system("Makephase < makephase.inp > makephase.out")
        os.system("Normxsc < normxsc.inp > normxsc.out")
        os.chdir(cwd)

    def _generate_fcloud_ref(self):
        """Generate the fcloud.ref file
        
        Args:
            None
            
        Returns:
            None
            
        Creates:
            fcloud.ref"""
        
        # Return immediately if there are no cloud layers defined
        if self.num_aerosol_modes == 0:
            return 
        
        with open(self.directory / "fcloud.ref", mode="w+") as file:
            file.write(f"{len(self.ref.data)}    {self.num_aerosol_modes}\n")
            for height in self.ref.data.height:
                file.write(f"      {height:> 9.4f}      {int(self.cloud_cover)}")
                for n in range(self.num_aerosol_modes):
                    file.write("      1")
                file.write("\n")

    def _generate_parah2_ref(self):
        """Generate the parah2.ref file from a profile in the eleos/data/[planet] directory, with any
        pressure limits imposed by NemesisCore.set_pressure_limits()
        
        Args:
            None
            
        Returns:
            None
            
        Creates:
            parah2.ref"""
        df = pd.read_table(constants.PATH / f"data/{self.planet}/parah2.ref", skiprows=1, header=None, sep="\s+")
        df.columns = ["height", "parah2"]
        minh, maxh = self.get_height_limits()
        df = df[(df.height >= minh) & (df.height <= maxh)]
        with open(self.directory / "parah2.ref", mode="w+") as file:
            file.write(str(len(df)) + "\n")
            file.write(df.to_string(header=False, index=False))

    def _generate_summary(self):
        """Dump the information used to create the NEMESIS core to a human-readable text file
        
        Args:
            None
            
        Returns:
            None
            
        Creates:
            eleos_inputs.txt"""
        out = ""
        for k, v in self.__dict__.items():
            if k == "profiles":
                out += "profiles:"
                for p in v:
                    out += utils.indent(p.__repr__(), level=1)
                    out += utils.indent(p.compact_str(), level=2)
            else:
                out += f"{k}: {v}"
            out += "\n"

        with open(self.directory / "eleos_inputs.txt", mode="w+") as file:
            file.write(out)

    def get_aerosol_mode(self, id=None, label=None):
        """Given an aerosol ID (positive, as used by NEMESIS), or an eleos label, return the corresponding AerosolProfile object.
        Note that as there is no uniqueness requirement on eleos labels, this function will return the Profile with the lowest aerosol ID
        in the case of multiple similar labels. This may be updated in the future to return a list of all matching profiles.
        
        Args:
            id (int): The aerosol ID
            label (str): The label associated with the AerosolProfile object
            
        Returns:
            AerosolProfile: The corresponding aerosol profile"""

        if id > self.num_aerosol_modes:
            raise IndexError(f"Aerosol index is too large! ({id} vs {self.num_aerosol_modes} max.)")
        for profile in self.profiles:
            if isinstance(profile, profiles_.AerosolProfile):
                if id is not None and profile.aerosol_id == id:
                    return profile
                if label is not None and profile.label == label:
                    return profile

    def add_aerosol_mode(self, aerosol_profile=None, refractive_index_profile=None, radius=None, variance=None, refractive_index=None):
        """Add an aerosol mode. This uses the Mie scattering option with constant refractive index over the wavelength range, which is set to be the same as the input .spx files.
        Can be called with either of the following signatures:

        core.add_aerosol_mode(AerosolProfile(...), radius, variance, refractive_index)
        or
        core.add_aerosol_mode(AerosolProfile(...), ImagRefractiveIndexProfile(...))

        to pull parameters from an ImagRefractiveIndexProfile or specified manually respectively
        
        Args:
            refractive_index_profile (ImagRefractiveIndexProfile): If given, use the aerosol radius distribution and refractive index from the Profile specified
            radius (float): The mean radius of the aerosol particles
            variance (float): The variance of the particle size distribution
            refractive_index (complex): The aerosol's refractive index as a complex number (eg. 1.2 + 1e-3j)
            
        Returns:
            None"""

        # Increment aerosol mode counter        
        self.num_aerosol_modes += 1

        # Add the profile to the core
        self._add_profile(aerosol_profile)

        # Assign aerosol IDs
        aerosol_profile.aerosol_id = self.num_aerosol_modes

        # Pull parameters from the refractive index profile and add to the core if given
        if refractive_index_profile is not None:
            refractive_index_profile.aerosol_id = self.num_aerosol_modes
            refractive_index_profile.shape.aerosol_id = self.num_aerosol_modes
            self._add_profile(refractive_index_profile)

            radius = refractive_index_profile.shape.radius
            variance = refractive_index_profile.shape.variance
            refractive_index = refractive_index_profile.shape.refractive_index

        # Add the mode to the bottom of the Makephase input file
        with open(self.directory / "makephase.inp", mode="a") as file:
            file.write(f"1\n{radius} {variance}\n2\n1\n{refractive_index.real} {refractive_index.imag}\n")

    def generate_core(self):
        print(f"Generating core {self.id_}")

        # Consistency check for aerosol mode ID - sometimes this screws up and I'm not sure why
        ids = []
        for p in self.profiles:
            if isinstance(p, profiles_.AerosolProfile):
                ids.append(p.aerosol_id)
        if len(set(ids)) != len(ids):
            raise ValueError("Aerosol IDs are not consistent!")

        self._generate_inp()
        self._generate_set()
        self._generate_flags()
        self._generate_aerosol_ref()
        self._generate_parah2_ref()
        self._generate_kls()
        self._generate_fmerror()
        self._generate_fcloud_ref()
        self._generate_xsc()
        self._generate_apr()
        self._generate_summary()
        self._save_core()        

    def fix_peak(self, central_wavelength, width, error=1e-15):
        self.fixed_peaks.append(FixedPeak(central_wavelength, width, error))

    def get_height_limits(self):
        """Get the heights of the top and bottom layers of the atmosphere in km.
        
        Args:
            None
            
        Returns:
            top_height: Height of the top of the atmosphere
            bottom_height: Height at the bottom of the atmosphere"""
        
        heights = self.ref.data.height
        return min(heights), max(heights)

    def set_pressure_limits(self, min_pressure=None, max_pressure=None):
        self.layer_type = 4
        self.min_pressure = min_pressure
        self.max_pressure = max_pressure
        pressures = np.logspace(np.log10(max_pressure), 
                                np.log10(min_pressure), 
                                self.num_layers)
        with open(self.directory / "pressure.lay", mode="w+") as file:
            file.write("Created by Eleos\n")
            file.write(str(self.num_layers) + "\n")
            for p in pressures:
                file.write(str(p) + "\n")


class FixedPeak:
    """Used internally to specify if any spectral regions should be fixed so that NEMESIS always fits it there. Don't instantiate, instead use NemesisCore.fix_peak"""
    def __init__(self, central_wavelength, width, error):
        self.central_wavelength = central_wavelength
        self.width = width
        self.error = error

    def __str__(self):
        return f"<FixedPeak at {self.central_wavelength}±{self.width/2:.4f}um>"

    def isin(self, wl):
        return (wl > self.central_wavelength - self.width/2) & (wl < self.central_wavelength + self.width/2)
    

def load_core(core_directory):
    """Load a NemesisCore object saved using NemesisCore._save_core
    
    Args:
        None
        
    Returns:
        NemesisCore: The unpickled core"""
    with open(Path(core_directory) / "core.pkl", 'rb') as f:
        core = pickle.load(f)

    # Refresh the directory attributes in the NemesisCore object in case the folder has been moved
    core.parent_directory = (core_directory / "..").resolve()
    core.directory = core_directory.resolve()

    return core


def reset_core_numbering():
    """Reset the automatic core numbering. Useful for creating multiple sets of cores in one program.
    
    Args:
        None
        
    Returns:
        None"""
    global CORE_ID_COUNTER
    CORE_ID_COUNTER = 0


def clear_parent_directory(parent_directory, prompt_if_exists=True):
    """Attempt to delete all the child files/direcotries from the specified folder. It is recommended to call this function at the start
    of every core generation script as it ensures that the script is fully stateless. If a script generates N cores when first run
    and is subsequently modified to produce N/2 cores, then the remaining N/2 core directories will remain in the parent directory.
    This should not matter as the SLURM script will only run the correct subset of cores, but it can get confusing. If the parent directory
    specified does not exist then it will print a warning and return.
    
    Args:
        parent_directory (str): The directory to clear
        prompt_if_exists (bool): Whether to confirm with the user before deleting the directory
        
    Returns:
        None"""
    
    parent_directory = Path(parent_directory)
    if os.path.exists(parent_directory):
        if prompt_if_exists:
            x = input(f"{parent_directory} already exists - are you sure you want to erase and continue? Y/N")
            if x.lower() == "n":
                print("Quitting")
                exit()
            else:
                shutil.rmtree(parent_directory)
        else:
            shutil.rmtree(parent_directory)


def generate_alice_job(parent_directory, python_env_name, username=None, memory=16, hours=24):
    """Generate an sbatch submission script for use on ALICE. The job is an array job over all the specified cores.
    After running NEMESIS it will run Eleos to create some summary plots in the parent/core_N/plots directory
    
    Args:
        parent_directory (str): The directory containing the NEMESIS core directories to be run
        python_env_name (str): The name of a conda environment which has Eleos installed
        username: The username of the user running the job (eg, scat2, lnf2)
        memory: The amount of memory to use (in GB)
        hours: The number of hours to schedule the job for
        
    Returns:
        None
        
    Creates:
        submitjob.sh in the parent directory of the cores"""

    parent_directory = Path(parent_directory)
    script_path = parent_directory / "submitjob.sh"

    ncores = len([f for f in os.listdir(parent_directory) if not os.path.isfile(parent_directory / f)])

    # Read the submission script and replace template fields
    with open(constants.PATH / "data/statics/template.job", mode="r") as file:
        out = file.read()
        out = out.replace("<MEMORY>", str(memory))
        out = out.replace("<HOURS>", f"{hours:02}")
        out = out.replace("<N_CORES>", str(ncores))
        out = out.replace("<CORE_DIR>", str(parent_directory.resolve()))
        out = out.replace("<USERNAME>", str(username))
        out = out.replace("<PYTHON_ENV>", python_env_name)

    # Write the filled template 
    with open(script_path, mode="w+") as file:
        file.write(out)


def run_alice_job(parent_directory, print_queue_delay=2):
    """Run the script created by generate_alice_job
    
    Args:
        parent_directory (str): Path to the directory containing all the cores
        print_queue_delay (float): Duration to wait (in seconds) before printing the current job queue

    Returns:
        None
    """
    os.system(f"sbatch {parent_directory}submitjob.sh")
    time.sleep(print_queue_delay)
    os.system("squeue --me")
