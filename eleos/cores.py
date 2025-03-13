import math
import pickle
import time
import sys
import signal
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import numpy as np

import warnings
warnings.formatwarning = lambda msg, *_: f"Warning: {msg}\n"

from . import constants
from . import spx
from . import parsers
from . import utils
from . import profiles as profiles_ # to avoid namespace shadowing by NemesisCore.profiles
from . import results

CORE_ID_COUNTER = 0


class NemesisCore:
    """The class that constructs a core directory and all the required files for NEMESIS to run.
    
    Attributes:
        parent_directory (Path):      The directory in which the core directory is created
        directory (Path):             The directory in which the core files are stored
        id_ (int):                    The unique ID of the core
        profiles (dict):              A dictionary of profiles.Profile objects to retrieve, where the keys are the profile labels and the values are the Profile objects themselves
        spx_file (Path):              The path to the spectrum file to fit
        ref_file (Path):              The path to the ref file to use
        ref (parsers.NemesisRef):     The parsed ref file
        planet (str):                 The name of the planet being observed
        scattering (bool):            Whether to run a scattering retrieval or not
        num_aerosol_modes (int):      The number of aerosol profiles in the retrieval
        forward (bool):               Whether of not to run a forward model (ie. set number of iterations = 0)
        num_iterations (int):         Number of iterations to run in the retrieval (if forward is set this has no effect)
        num_layers (int):             The number of atmospheric layers to simulate
        layer_type (int):             The type of layering to use 
        min_pressure (float):         The minimum pressure in the atmosphere
        max_pressure (float):         The maximum pressure in the atmosphere
        bottom_layer_height (int):    The height in km of the bottom of the atmosphere (by defauylt use the lowest height in the .ref file)
        fixed_peaks (list):           A list of FixedPeak objects that specify regions of the spectrum to fix
        instrument_ktables (str):     Either 'NIRSPEC' or 'MIRI'; determines which set of ktables to use.
        fmerror_factor (float):       The factor by which to multiply the error on the spectrum (see also, fmerror_pct and fmerror_value)
        fmerror_pct (float):          If given, instead of using fmerror_factor or fmerror_value, use a flat percentage of the brightness (eg. 0.1 = 10%) (see also, fmerror_factor and fmerror_value)
        fmerror_value (float):        If given, instead of using fmerror_factor or fmerror_pct, use a flat value in W/cm2/sr/um (see also, fmerror_factor and fmerror_pct)
        cloud_cover (bool):           If scattering mode is on, then this is the fractional cloud cover between 0 and 1 (usually doesn't need to be changed)
        reference_wavelength (float): If scattering mode is on, then normalise the cross-sections at the closest wavelength to this value in the .xsc file
        """
    def __init__(self, 
                 parent_directory,
                 spx_file, 
                 ref_file=None, 
                 profiles=list(), 
                 planet="jupiter", 
                 scattering=True, 
                 forward=False,
                 confirm=True, 
                 num_iterations=30,
                 num_layers=120, 
                 bottom_layer_height=None, 
                 min_pressure=None,
                 max_pressure=None,
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
            confirm (bool):               If the core directory already exists, then ask before continuing
            num_iterations (int):         Number of iterations to run in the retrieval (if forward is set this has no effect)
            num_layers (int):             The number of atmospheric layers to simulate
            bottom_layer_height (int):    The height in km of the bottom of the atmosphere (by default it uses the lowest height in the .ref file)
            min_pressure (float):         The pressure at the top of the model atmosphere
            max_pressure (float):         The pressure at the bottom of the model atmosphere
            instrument_ktables (str):     Either 'NIRSPEC' or 'MIRI'; determines which set of ktables to use.
            fmerror_factor (float):       The factor by which to multiply the error on the spectrum (see also, fmerror_pct and fmerror_value)
            fmerror_pct (float):          If given, instead of using fmerror_factor or fmerror_value, use a flat percentage of the brightness (eg. 0.1 = 10%) (see also, fmerror_factor and fmerror_value)
            fmerror_value (float):        If given, instead of using fmerror_factor or fmerror_pct, use a flat value in W/cm2/sr/um (see also, fmerror_factor and fmerror_pct)
            cloud_cover (bool):           If scattering mode is on, then this is the fractional cloud cover between 0 and 1 (usually doesn't need to be changed)
            reference_wavelength (float): If scattering mode is on, then normalise the cross-sections at the closest wavelength to this value in the .xsc file. This parameter will be updated with the exact wavelength
        """
        # Increment the global core counter and store local version
        global CORE_ID_COUNTER
        CORE_ID_COUNTER += 1
        self.id_ = CORE_ID_COUNTER

        # Assign attributes passed in
        self.spx_file = Path(spx_file).resolve()
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
        self._create_directory_tree(confirm)

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

        # Parse the ref file:
        self.ref = parsers.NemesisRef(self.ref_file)

        # Set bottom layer height if not defined
        if bottom_layer_height is None:
            self.bottom_layer_height = self.ref.data.iloc[0].height
        else:
            self.bottom_layer_height = bottom_layer_height

        # Set min/max pressures
        if min_pressure is None:
            self.min_pressure = self.ref.data.pressure.min()
        if max_pressure is None:
            self.max_pressure = self.ref.data.pressure.max()

        # Set number of aerosol modes (incremented by add_profile)
        self.num_aerosol_modes = 0

        # Add a list to hold any FixedPeak classes the user wants
        self.fixed_peaks = []

        # Add a reference to self in each Profile and set up AerosolProfiles correctly
        self.profiles = dict()
        for profile in profiles:
            self.add_profile(profile)

        # If in forward mode, set the number of iterations to 0
        if self.forward:
            self.num_iterations = 0

        # Raise an error if trying to use features not implemented yet
        if planet != "jupiter":
            raise Exception("Eleos does not support planets other than Jupiter yet! ")

    def __str__(self):
        return f"<NemesisCore: {self.directory}>"

    def _create_directory_tree(self, prompt_if_exists):
        os.makedirs(self.parent_directory, exist_ok=True)
        if os.path.exists(self.directory) and prompt_if_exists:
            if os.path.exists(self.directory / "nemesis.mre"):
                msg = "There is already a core that has been run in"
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

    def _save_core(self):
        """Dump the core object to a pickle file in the core directory
        
        Args:
            None
        
        Returns:
            None
            
        Creates:
            core.pkl"""
        with open(self.directory / "core.pkl", mode="wb") as file:
            pickle.dump(self, file)

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
        """Copy some boilerplate files from the library data directory into the core directory
        
        Args:
            None
            
        Returns:
            None
            
        Creates:
            nemesis.cia
            nemesis.abo
            nemesis.nam
            nemesis.sol
            maskephase.inp"""
        shutil.copy(constants.PATH / "data/statics/nemesis.cia", self.directory)
        shutil.copy(constants.PATH / "data/statics/nemesis.abo", self.directory)
        shutil.copy(constants.PATH / "data/statics/nemesis.nam", self.directory)
        shutil.copy(constants.PATH / "data/statics/nemesis.sol", self.directory)

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
        
        out = f"*******Apriori File*******\n           <NUM_PROFILES>\n"

        num_profiles = 0
        for profile in self.profiles.values():

            profile.shape.create_required_files(self.directory)

            if isinstance(profile, profiles_.AerosolProfile) and profile.retrieve_optical:
                profile._generate_cloudfn_dat(self.directory)
                num_profiles += 1

            out += profile.generate_apr_data() + "\n"
            num_profiles += 1

        with open(self.directory / "nemesis.apr", mode="w") as file:
            out = out.replace("<NUM_PROFILES>", str(num_profiles))
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
        out = out.replace("<DISTANCE>", f"{constants.DISTANCE[self.planet]:.3f}")
        out = out.replace("<SUNLIGHT>", f"{int(self.scattering)}")
        out = out.replace("<BOUNDARY>", f"{int(self.scattering)}")
        out = out.replace("<N_LAYERS>", f"{int(self.num_layers)}")
        out = out.replace("<BASE>", f"{self.bottom_layer_height:.2f}")
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
            nemesis.xsc
            PHASE{N}.DAT
            hgphase{n}.dat"""
        
        # Check to see if any aerosols have been added - if not then return immediately
        if self.num_aerosol_modes == 0:
            return

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

    def _generate_makephase_inp(self):
        """Add the optical properties of any attached AerosolProfiles to the Makephase input file
        
        Args:
            None
            
        Returns:
            None
            
        Creates:
            makephase.inp
            normxsc.inp"""
        
        # Get wavelength bounds including the specified reference wavelength
        start, end, step, idx = self._get_xsc_wavelengths(step=0.1)
        
        # Generate the makephase.inp file
        with open(self.directory / "makephase.inp", mode="w+") as file:
            file.write(str(self.num_aerosol_modes) + "\n")
            file.write("1\n")
            utils.write_nums(file, start, end, step)
            file.write("nemesis\n")
            file.write("y\n")

            for name, profile in self.profiles.items():
                if isinstance(profile, profiles_.AerosolProfile):
                    file.write(f"1\n{profile.radius} {profile.variance}\n2\n")
                    if not profile.lookup:
                        file.write(f"1\n{profile.real_n} {profile.imag_n}\n")
                    else:
                        file.write(f"{constants.MAKEPHASE_GASES[profile.n_lookup]}\n")

        # Generate the normxsc.inp file
        with open(self.directory / "normxsc.inp", mode="w+") as file:
            file.write(f"nemesis.xsc\n{idx + 1} 1")

    def _generate_pressure_lay(self):
        """Generate the pressure.lay file that splits the .ref file up into layers for the model
        
        Args:
            min_pressure (float): The pressure at the top of the atmosphere
            max_pressure (float): The pressure at the bottom of the atmosphere
            
        Returns:
            None
            
        Creates:
            pressure.lay"""
        
        pressures = np.logspace(np.log10(self.max_pressure), 
                                np.log10(self.min_pressure), 
                                self.num_layers)
        
        with open(self.directory / "pressure.lay", mode="w+") as file:
            file.write("Created by Eleos\n")
            file.write(str(self.num_layers) + "\n")
            for p in pressures:
                file.write(str(p) + "\n")

    def _generate_summary(self):
        """Dump the information used to create the NEMESIS core to a human-readable text file
        
        Args:
            None
            
        Returns:
            None
            
        Creates:
            eleos_inputs.txt
            aerosol_names.txt"""
        
        out = ""
        for k, v in self.__dict__.items():
            if k == "profiles":
                out += "profiles:\n"
                for label, profile in v.items():
                    out += utils.indent(profile.__repr__(), level=1)
                    out += "\n"
                    out += utils.indent("\n".join([f"{k}: {v}" for k,v in profile.__dict__.items()]), level=2)
                    out += "\n"
            else:
                out += f"{k}: {v}"
            out += "\n"

        with open(self.directory / "summary.txt", mode="w+") as file:
            file.write(f"Generated: {time.asctime()}\n")
            file.write(out)

        with open(self.directory / "aerosol_names.txt", mode="w+") as file:
            for name, profile in self.profiles.items():
                if isinstance(profile, profiles_.AerosolProfile):
                    file.write(name + "\n")
            file.truncate(file.tell() - 1)

    def _get_xsc_wavelengths(self, step):
        """Get the wavelengths for use in Makephase, including the full data range and, if necessary, 
        expanding to include the reference wavelength specified. 
        
        Args:
            step: Wavelength step size (microns)
            
        Returns:
            float: Start wavelength
            float: End wavelength
            float: Wavelength step
            int: reference wavelength index
        """

        # Get wavelengths from spx file
        wls = spx.read(self.spx_file).geometries[0].wavelengths

        # Initial guesses for start and end wavelengths
        start_wl = math.floor(min(wls) / step) * step
        end_wl = math.ceil(max(wls) / step) * step

        # Expand the range to include the reference wavelength if necessary
        if self.reference_wavelength < start_wl:
            start_wl = math.floor(self.reference_wavelength / step) * step  # Round down to nearest step
        elif self.reference_wavelength > end_wl:
            end_wl = math.ceil(self.reference_wavelength / step) * step  # Round up to nearest step
        
        # Find the reference wavelength index and set the reference value to the closest value
        idx, ref = utils.find_nearest(np.arange(start_wl, end_wl+step, step), self.reference_wavelength)
        self.reference_wavelength = ref

        return start_wl, end_wl, step, idx

    def generate_core(self, verbosity=1):
        """Create all the files necessary for a NEMESIS retrieval in the directory specified by NemesisCore.directory.
        
        Args:
            verbosity (int): One of 0, 1, 2. If 0 then do not print any progress. If 1 then print core number. If 2 then print status for each file
            
        Returns:
            None

        Creates:
            See _generate_* methods
        """
        if verbosity == 1: print(f"Generating core {self.id_}")

        if verbosity == 2: print(f"Generating summary file")
        self._generate_summary()

        if verbosity == 2: print(f"Copying .spx and .ref files")
        self._copy_input_files()

        if verbosity == 2: print(f"Copying boilerplate files")
        self._copy_template_files()

        if verbosity == 2: print(f"Generating pressure.lay")
        self._generate_pressure_lay()

        if verbosity == 2: print(f"Generating nemesis.inp")
        self._generate_inp()
        
        if verbosity == 2: print(f"Generating nemesis.set")
        self._generate_set()

        if verbosity == 2: print(f"Generating nemesis.fla")
        self._generate_flags()

        if verbosity == 2: print(f"Generating aerosol.ref")
        self._generate_aerosol_ref()

        if verbosity == 2: print(f"Generating parah2.ref")
        self._generate_parah2_ref()

        if verbosity == 2: print(f"Generating nemesis.kls")
        self._generate_kls()

        if verbosity == 2: print(f"Generating fmerror.dat")
        self._generate_fmerror()

        if verbosity == 2: print(f"Generating fcloud.ref")
        self._generate_fcloud_ref()

        if verbosity == 2: print(f"Generating makephase.inp and normxsc.inp")
        self._generate_makephase_inp()

        if verbosity == 2: print(f"Generating nemesis.xsc and [hg]phase.dat files")
        self._generate_xsc()

        if verbosity == 2: print(f"Generating nemesis.apr")
        self._generate_apr()

        if verbosity == 2: print("Saving Eleos core object")
        self._save_core()        

    def run(self):
        """Run NEMESIS on the core. This method should only be used for short forward models as
        it does not schedule the jobs on the ALICE compute nodes (see run_alice_job() for that), 
        instead running it on the login node.
        
        Args:
            None
            
        Returns:
            NemeisResult: The results object for the core"""
        
        self.generate_core(print_progress=True)
        print("Running NEMESIS...")
        subprocess.run("Nemesis < nemesis.nam > nemesis.prc", 
                       shell=True, 
                       cwd=self.directory,
                       stdout=subprocess.PIPE)
        res = results.NemesisResult(self.directory)
        res.make_summary_plot()
        print("Finished!")
        return res

    def generate_prior_distributions(self):
        """Run NEMESIS briefly in a temporary directory to generate the prior gas and aerosol profiles.
        
        Args:
            None
            
        Returns:
            pd.DataFrame: Data from the .prf files generated. Columns are: 'height', 'pressure', 'temperature', the gas profiles, and the aerosol profiles
        """

        # Change core directory to eleos/tmp temporarily
        old_dir = self.directory
        self.directory = constants.PATH / "tmp"

        # Clear tmp directory
        for fp in self.directory.glob("*"):
            if fp.is_file():
                fp.unlink()

        # Copy any extra files from the mian directory
        for filename in old_dir.glob("*"):
            if filename.is_file():
                shutil.copy(filename, self.directory)

        # Generate the core in this directory
        self.generate_core(verbosity=0)

        # Run NEMESIS
        with open(self.directory / "nemesis.nam") as inp:
            proc = subprocess.Popen("Nemesis < nemesis.nam > nemesis.prc", 
                                    shell=True, 
                                    cwd=self.directory,
                                    preexec_fn=os.setsid)
            # Wait for the aerosol.prf file to be generated - potential race condition??
            while True:
                time.sleep(0.05)
                if (self.directory / "aerosol.prf").exists() and (self.directory / "nemesis.prf").exists():
                    break
        
        # Kill NEMESIS
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM)

        # Read in the .prf files
        aerosol_prf = parsers.AerosolPrf(self.directory / "aerosol.prf")
        nemesis_prf = parsers.NemesisPrf(self.directory / "nemesis.prf")

        # Reset the directory and clear tmp
        for fp in self.directory.glob("*"):
            if fp.is_file():
                fp.unlink()
        self.directory = old_dir

        nemesis_prf.data.height = np.round(nemesis_prf.data.height * 1000).astype(int) 
        aerosol_prf.data.height = np.round(aerosol_prf.data.height * 1000).astype(int) 
        merged = pd.merge(nemesis_prf.data, aerosol_prf.data, on="height")
        merged.height /= 1000

        return merged

    def get_aerosol_profile(self, id=None):
        """Given an aerosol ID (as used by NEMESIS) return the corresponding AerosolProfile object.
        
        Args:
            id (int): The aerosol ID
            label (str): The label associated with the AerosolProfile object
            
        Returns:
            AerosolProfile: The corresponding aerosol profile"""

        id = abs(id)

        if id > self.num_aerosol_modes:
            raise IndexError(f"Aerosol index is too large! ({id} vs {self.num_aerosol_modes} max.)")

        for profile in self.profiles.values():
            if isinstance(profile, profiles_.AerosolProfile):
                if profile.aerosol_id == id:
                    return profile

    def add_profile(self, profile):
        """Add a profile to retrieve. It's recommended to do this during instantiation using the profiles argument.

        Args:
            profile: profiles.Profile object to add

        Returns:
            None"""
        
        self.profiles[profile.label] = profile
        profile.core = self

        if isinstance(profile, profiles_.AerosolProfile):
            # Increment aerosol mode counter        
            self.num_aerosol_modes += 1

            # Assign aerosol IDs
            profile.aerosol_id = self.num_aerosol_modes

    def fix_peak(self, central_wavelength, width, error=1e-15):
        """Set a region of the spectra to be fixed (ie. NEMESIS will be forced to match it). This is done internally by setting the
        error on that region to be extrmemely small.
        
        Args:
            central_wavelength (float): The central wavelength of the peak
            width (float): The width of the spectral feature (this is a full width, so the region being fixed is centre-width/2 to centre+width/2)
            error (float): Optionally, give the error to be assigned to the region. This should be very small and there is not much need to change this
        
        Returns:
            None
        """
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

    def set_random_priors(self):
        """For each parameter in the associated profiles, set the values randomly based on a normal distribution with mean and standard deviation given by the pre-existing values. 
        
        For example, if there is a single attached GasProfile with the following parameters:
        GasProfile(gas_name="PH3", 
                   shape=shapes.Shape20(
                   knee_pressure=1.0, 
                   tropopause_pressure=0.1,
                   deep_vmr=1.86e-6,          deep_vmr_error=0.2e-6,
                   fsh=0.3,                   fsh_error=0.1))
        then calling this function will set deep_vmr to a value randomly sampled from a normal distribution with mean 1.86e-6 and standard deviation 0.2e-6, and 
        fsh to a value randomly sampled from a normal distribution with mean 0.3 and standard deviation 0.1. The core can then be modified / generated as normal.
        
        Args:
            None

        Returns:
            None
        """
        print("Generating random priors...")
        for label, profile in self.profiles.items():
            for name in profile.NAMES:
                try:
                    mean = getattr(profile, name)
                    sigma = getattr(profile, name + "_error")
                    n = -1
                    while n < 0:
                        n = np.random.normal(mean, sigma)
                    if sigma == 0:
                        setattr(profile, name+"_error", 1e-8)
                    setattr(profile, name, n)
                except AttributeError:
                    continue
        print("Generated priors!")


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
    
    core_directory = Path(core_directory)

    with open(core_directory / "core.pkl", 'rb') as f:
        core = pickle.load(f)

    # Refresh the directory attributes in the NemesisCore object in case the folder has been moved
    core.parent_directory = (core_directory / "..").resolve()
    core.directory = core_directory.resolve()

    return core


def load_from_previous(previous_directory, parent_directory, confirm=True):
    """Load a core from a previous retrieval to use as a template for creating a new core. This method loads the core
    using the core.pkl file and then calls from_previous_retrieval on all Profile objects stored, as well as resetting core ID,
    parent_directory, etc...
    
    Args:
        previous_directory (str): The path to the core directory to load
        parent_directory: THe new parent directory to create the new core under
        
    Returns:
        NemesisCore: The new core object"""
    
    core = load_core(previous_directory)
    res = results.NemesisResult(previous_directory)

    global CORE_ID_COUNTER
    CORE_ID_COUNTER += 1
    core.id_ = CORE_ID_COUNTER

    core.parent_directory = Path(parent_directory)
    core.directory = core.parent_directory / f"core_{core.id_}"
    core._create_directory_tree(confirm)

    for label, profile in core.profiles.items():
        core.profiles[label] = profile.from_previous_retrieval(res, label)

    return core


def reset_core_numbering(): 
    """Reset the automatic core numbering. Useful for creating multiple sets of cores in one program.
    
    Args:
        None
        
    Returns:
        None"""
    global CORE_ID_COUNTER
    CORE_ID_COUNTER = 0


def clear_parent_directory(parent_directory, confirm=True):
    """Attempt to delete all the child files/directories from the specified folder. It is recommended to call this function at the start
    of every core generation script as it ensures that the script is fully stateless. If a script generates N cores when first run
    and is subsequently modified to produce N/2 cores, then the remaining N/2 core directories will remain in the parent directory.
    This should not matter as the SLURM script will only run the correct subset of cores, but it can get confusing. If the parent directory
    specified does not exist then it will print a warning and return.
    
    Args:
        parent_directory (str): The directory to clear
        confirm (bool): Whether to confirm with the user before deleting the directory
        
    Returns:
        None"""
    
    parent_directory = Path(parent_directory)
    if os.path.exists(parent_directory):
        if confirm:
            x = input(f"{parent_directory} already exists - are you sure you want to erase and continue? Y/N")
            if x.lower() == "n":
                print("Quitting")
                exit()
            else:
                shutil.rmtree(parent_directory)
        else:
            shutil.rmtree(parent_directory)


def generate_alice_job(parent_directory, 
                       python_env_name, 
                       memory=16, 
                       hours=24,
                       username=None,
                       notify=("end", "fail")):
    """Generate an sbatch submission script for use on ALICE. The job is an array job over all the specified cores.
    After running NEMESIS it will run Eleos to create some summary plots in the parent/core_N/plots directory
    
    Args:
        parent_directory (str): The directory containing the NEMESIS core directories to be run
        python_env_name (str): The name of a conda environment which has Eleos installed
        memory (int): The amount of memory to use (in GB)
        hours (int): The number of hours to schedule the job for
        username (str): The username of the user running the job (eg, scat2, lnf2)
        notify (List[str] or str): What email notifications to send. Can be any combination of 'begin', 'end', 'fail'
        
    Returns:
        None
        
    Creates:
        submitjob.sh in the parent directory of the cores"""

    parent_directory = Path(parent_directory)
    script_path = parent_directory / "submitjob.sh"
    ncores = len([f for f in os.listdir(parent_directory) if not os.path.isfile(parent_directory / f)])

    # Coerce notify to an iterable and create the string
    try:
        len(notify)
    except:
        notify = [notify]
    notify_str = ",".join(notify).upper()   

    # Assert that hours and memoryt are ints
    assert isinstance(hours, int)
    assert isinstance(memory, int) 

    # Read the submission script and replace template fields
    with open(constants.PATH / "data/statics/template.job", mode="r") as file:
        out = file.read()
        out = out.replace("<MEMORY>", str(memory))
        out = out.replace("<HOURS>", f"{hours:02}")
        out = out.replace("<N_CORES>", str(ncores))
        out = out.replace("<CORE_DIR>", str(parent_directory.resolve()))
        out = out.replace("<USERNAME>", str(username))
        out = out.replace("<PYTHON_ENV>", python_env_name)
        out = out.replace("<NOTIFY>", notify_str)

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
    parent_directory = Path(parent_directory)
    os.system(f"sbatch {parent_directory / 'submitjob.sh'}")
    time.sleep(print_queue_delay)
    os.system("squeue --me")
