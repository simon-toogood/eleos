import os
import shutil
import pandas as pd
import pickle
import time
import warnings

from . import constants
from . import spx


# Files:
#   Thermal:
#     * aerosol.ref (aerosol profiles as function of height)
#     * parah2.ref  (don't know, just copy it)
#     * fmerror.dat (additional error for each wavelength)
#     * nemesis.fla (various flags)
#     * nemesis.cia (pointer to collision-induced opacity file)
#     * nemesis.kls (list of k-tables)
#     * nemesis.ref (reference heights, pressures, temps and composition)
#     * nemesis.xsc (aerosol cross-sections)
#     * nemesis.set (scattering angles, layer info)
#     * nemesis.apr (a priori file containing vars to retrieve)
#     * nemesis.spx (spectrum)
#     * nemesis.inp (input flags)
#    ** .abo, .nam
#   

# NOTES:
# aerosol radius for xsc can now be specified in two places; ImagrefractiveIndexProfile and here.
#  What happens if theyre differnet? Do i need to add a check for the presence of a 444 and overwrite that?


class NemesisCore:
    core_id = 0

    def __init__(self, 
                 parent_directory, 
                 spx_file, 
                 ref_file=None, 
                 profiles=list(), 
                 planet="jupiter", 
                 scattering=True, 
                 forward=False,
                 num_iterations=15,
                 num_layers=120, 
                 bottom_layer_height=-80, 
                 num_aerosol_modes=1, 
                 instrument="NIRSPEC", 
                 fmerror_factor=3,
                 aerosol_radius=2, 
                 aerosol_variance=0.1, 
                 aerosol_refractive_index=1.3+1e-3j,
                 cloud_cover=1.0):
        """Create a NEMESIS core directory with a given set of profiles to retrieve
        
        Args:
            parent_directory: The directory in which to create the core folder 
            spx_file: Path to the spectrum file to fit
            ref_file: Path to the ref file to use (if left blank it will use the default for that planet)
            profiles: List of profiles.Profile objects to retrieve
            planet: Name of the planet being observed. Must be one of 'jupiter', 'saturn', 'uranus', 'neptune' or 'titan'.
            scattering: Whether to run a scattering retrieval or not
            forward: Whether of not to run a forward model (ie. set number of iterations = 0)
            num_iterations: Number of iterations to run in the retrieval (if forward is set this has no effect)
            num_layers: The number of atmospheric layers to simulate
            bottom_layer_height: The height in km of the bottom layer of the atmosphere
            num_aerosol_modes: The number of aerosol modes to use
            instrument: Either 'NIRSPEC' or 'MIRI'; determines which set of ktables to use
            fmerror_factor: The factor by which to multiply the error on the spectrum
            aerosol_radius: The mean radius of the aerosol particles (this may change in future for many clouds)
            aerosol_variance: The variance of the aerosol size distributions (this may change in future for many clouds), 
            aerosol_refractive_index: A complex value for the refractive inde of the aersol (this may change in future for many clouds),
            cloud_cover: If scattering, then this is the fractional cloud cover between 0 and 1 (should not usually be changed)
        """
        # Increment the global core counter
        NemesisCore.core_id += 1

        # Set the directories of the parent folder and own cores
        self.parent_directory = parent_directory
        self.directory = parent_directory + f"core_{self.core_id}/"

        # Set ref file if not specified:
        if ref_file is None:
            self.ref_file = constants.PATH + f"data/{planet}/{planet}.ref"
            warnings.warn(f"No ref file specified. Using the default in {self.ref_file}")
        else:
            self.ref_file = ref_file

        # Assign attributes passed in
        self.spx_file = spx_file
        self.profiles = profiles
        self.planet = planet.lower()
        self.scattering = scattering
        self.forward = forward
        self.num_iterations = num_iterations
        self.num_layers = num_layers
        self.bottom_layer_height = bottom_layer_height
        self.num_aerosol_modes = num_aerosol_modes
        self.instrument = instrument
        self.fmerror_factor = fmerror_factor
        self.aerosol_radius = aerosol_radius
        self.aerosol_variance = aerosol_variance
        self.aerosol_refractive_index = aerosol_refractive_index
        self.cloud_cover = cloud_cover

        # Add a reference to self in each Profile
        for profile in profiles:
            profile.core = self

        # If in forward mode, set the number of iterations to 0
        if self.forward:
            self.num_iterations = 0

        # Parse the ref file:
        self.ref = parse_ref_file(self.ref_file)

        # Raise an error if trying to use features not implemented yet
        if planet != "jupiter":
            raise Exception("Eleos does not fully support planets other than Jupiter yet!")
        if num_aerosol_modes > 1:
            raise Exception("Eleos does not support multiple aerosol modes yet!")
        
        # Create the directory tree if it doesn't already exist and clear it if it does
        os.makedirs(self.parent_directory, exist_ok=True)
        if os.path.exists(self.directory):
            shutil.rmtree(self.directory)
        os.makedirs(self.directory)

        # Create the core!
        self.generate_core()

    def __str__(self):
        return f"<NemesisCore: {self.directory}>"

    def _save_core(self):
        with open(self.directory+"core.pkl", mode="wb") as file:
            pickle.dump(self, file)

    def _copy_input_files(self):
        """Copy the given .spx and .ref file into the core
        
        Args:
            None
            
        Returns:
            None
            
        Creates:
            nemesis.spx
            nemesis.ref"""

        shutil.copy(self.ref_file, self.directory+"nemesis.ref")
        shutil.copy(self.spx_file, self.directory+"nemesis.spx")

    def _copy_template_files(self):
        shutil.copy(constants.PATH + "data/statics/nemesis.cia", self.directory)
        shutil.copy(constants.PATH + "data/statics/nemesis.abo", self.directory)
        shutil.copy(constants.PATH + "data/statics/nemesis.nam", self.directory)
        shutil.copy(constants.PATH + "data/statics/nemesis.sol", self.directory)
        shutil.copy(constants.PATH + f"data/{self.planet}/parah2.ref" , self.directory)

    def _generate_inp(self):
        """Generate the nemesis input file

        Args:
            num_iterations: Number of iterations to run (has no effect if NemesisCore.forward is True)
        
        Returns:
            None
            
        Creates:
            nemesis.inp"""
        with open(constants.PATH+"data/statics/nemesis.inp") as file:
            out = file.read()
            out = out.replace("<SCATTERING>", str(int(self.scattering)))
            out = out.replace("<N_ITERATIONS>", str(self.num_iterations))
        with open(self.directory+"nemesis.inp", mode="w+") as file:
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
        with open(self.directory + "nemesis.apr", mode="w") as file:
            file.write(out)

    def _generate_set(self):
        """Generate the settings file for NEMESIS
        
        Args:
            num_layers: The number of atmospheric layers
            bottom_layer_height: The altitude of the bottom layer of the atmosphere
            
        Returns:
            None
            
        Creates:
            nemesis.set"""

        with open(constants.PATH+"data/statics/template.set", mode="r") as file:
            out = file.read()
        out = out.replace("<DISTANCE>", f"{constants.DISTANCES[self.planet]:.3f}")
        out = out.replace("<SUNLIGHT>", f"{int(self.scattering)}")
        out = out.replace("<BOUNDARY>", f"{int(self.scattering)}")
        out = out.replace("<BASE>", f"{self.bottom_layer_height:.2f}")
        with open(self.directory+"nemesis.set", mode="w+") as file:
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
        
        with open(self.directory+"nemesis.fla", mode="w+") as file:
            file.write('       {aa}	! Inormal (0=eqm, 1=normal)\n       {bb}	! Iray (0=off, 1=on)\n       {cc}	! IH2O (1 = turn extra continuum on)\n       {dd}	! ICH4 (1 = turn extra continuum on)\n       {ee}	! IO3 (1 = turn extra continuum on)\n       {ff}	! INH3 (1 = turn extra continuum on)\n       {gg}	! Iptf (0=default, 1=CH4 High-T)\n       {hh}	! IMie\n       {ii}	! UV Cross-sections\n       0	! INLTE (0=LTE)'.format(aa=inormal, bb=iray, cc=ih2o, dd=ich4, ee=io3, ff=inh3, gg=iptf, hh=imie, ii=iuv))

    def _generate_aerosol_ref(self):
        """Generate a number of aerosol reference profile that is 0 at all altitudes.
        
        Args:
            num_aerosol_modes: Number of aerosol modes to generate
            
        Returns:
            None
            
        Creates:
            aerosol.ref"""
        heights = self.ref.height
        with open(self.directory+"aerosol.ref", mode="w+") as file:
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
        if self.instrument == "NIRSPEC":
            shutil.copy(constants.PATH+"data/jupiter/nirspec.kls", self.directory+"nemesis.kls")
        elif self.instrument == "MIRI":
            shutil.copy(constants.PATH+"data/jupiter/miri.kls", self.directory+"nemesis.kls")

    def _generate_fmerror(self):
        """For each wavelength in the spx file, multiply the error by factor. Currently only supports 1 spx geometry
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
        with open(self.directory + "fmerror.dat", mode="w+") as file:
            # Header with number of lines 
            file.write(f"{num_entries+2}\n")
            # Catch any cases outside the wavelength range (lower)
            file.write(f"{0.1:.6e}  {1e-8:.6e}\n")
            # Scale the spx error
            for wl, err in zip(spx_data.wavelengths, spx_data.error):
                file.write(f"{wl:.6e}  {err*self.fmerror_factor:.6e}\n")
            # Catch any cases outside the wavelength range (upper)
            file.write(f"{100:.6e}  {1e-8:.6e}\n")

    def _generate_xsc(self):
        """Run Makephase and Normxsc with the given inputs to generate the .xsc file containing aerosol refractive indicies.
        Very limited in functionality currently!!
        Currently only supports mode 1 (Mie scattering, gammma dist.) and constant refractive index over range.
        
        Args:
            None
            
        Returns:
            None
        
        Creates:
            makephase.inp
            normxsc.inp
            nemesis.xsc
            PHASE{N}.DAT
            hgphase{n}.dat
            
        TODO:
            Add suppport for multiple aerosol modes
            Add support for other aerosol scattering properties"""

        # Read in wavelengths bounds from spx file
        wls = spx.read(self.spx_file).geometries[0].wavelengths
        start_wl = min(wls) - 0.1
        end_wl = max(wls) + 0.1

        # Split the refractive index into real and imag
        real_n = self.aerosol_refractive_index.real
        imag_n = self.aerosol_refractive_index.imag

        # Generate the makephase.inp file
        with open(self.directory+"makephase.inp", mode="w+") as file:
            file.write(f"{self.num_aerosol_modes}\n1\n{start_wl} {end_wl} 0.1\nnemesis.xsc\ny\n1\n{self.aerosol_radius} {self.aerosol_variance}\n2\n1\n{real_n} {imag_n}")

        # Generate the normxsc.inp file
        with open(self.directory+"normxsc.inp", mode="w+") as file:
            file.write(f"nemesis.xsc\n1 1")

        # Run Makephase and Normxsc
        cwd = os.getcwd()
        os.chdir(self.directory)
        os.system("Makephase < makephase.inp")
        os.system("Normxsc < normxsc.inp")
        os.chdir(cwd)

    def _generate_fcloud_ref(self):
        """Generate the fcloud.ref file
        
        Args:
            None
            
        Returns:
            None
            
        Creates:
            fcloud.ref"""
        with open(self.directory+"fcloud.ref", mode="w+") as file:
            file.write(f"{len(self.ref)}    {self.num_aerosol_modes}\n")
            for height in self.ref.height:
                # This will break with more than 1 aerosol mode!
                file.write(f"      {height:> 9.4f}      {int(self.cloud_cover)}      1\n")

    def generate_core(self):
        self._copy_input_files()
        self._copy_template_files()
        self._generate_inp()
        self._generate_set()
        self._generate_flags()
        self._generate_aerosol_ref()
        self._generate_kls()
        self._generate_fmerror()
        self._generate_xsc()
        self._generate_apr()
        self._save_core()
        if self.scattering:
            self._generate_fcloud_ref()


def parse_ref_file(ref_file):
    """Read in the .ref file provided and return a DataFrame
    
    Args:
        ref_file: Path to the .ref file
        
    Returns:
        pandas.DataFrame with columns for height (in km), pressure (in atm), temperature (in K), and the VMRs for all specified gases"""
    with open(ref_file) as file:
        for i, line in enumerate(file):
            if "height" in line:
                skip_to = i
            if i == 2:
                tokens = line.split()
                n_gas = int(tokens[-1])
    df = pd.read_table(ref_file, skiprows=skip_to+1, sep="\s+", names=["height", "pressure", "temp"] + [f"VMR gas {n+1}" for n in range(n_gas)])
    return df


def load_core(core_directory):
    """Load a NemesisCore object saved using NemesisCore._save_core
    
    Args:
        None
        
    Returns:
        NemesisCore: The unpickled core"""
    with open(core_directory+"core.pkl", 'rb') as f:
        return pickle.load(f)


def reset_core_numbering():
    """Reset the automatic core numbering. Useful for creating multiple sets of cores in one program.
    
    Args:
        None
        
    Returns:
        None"""
    NemesisCore.core_id = 0


def clear_parent_directory(parent_directory):
    """Attempt to delete all the child files/direcotries from the specified folder. It is recommended to call this function at the start
    of every core generation script as it ensures that the script is fully stateless. If a script generates N cores when first run
    and is subsequently modified to produce N/2 cores, then the remaining N/2 core directories will remain in the parent directory.
    This should not matter as the SLURM script will only run the correct subset of cores, but it can get confusing. If the parent directory
    specified does not exist then it will print a warning and return.
    
    Args:
        parent_directory (str): The directory to clear
        
    Returns:
        None"""
    try:
        shutil.rmtree(parent_directory)
    except FileNotFoundError:
        warnings.warn(f"Attempted to clear directory {parent_directory} but it does not exist.")
        return


def generate_alice_job(cores, username, memory=16, hours=24):
    """Generate an sbatch submission script for use on ALICE. The job is an array job over all the specified cores.
    
    Args:
        cores: Either a single core or a list of cores
        username: The username of the user running the job (eg, scat2, lnf2)
        memory: The amount of memory to use (in GB)
        hours: The number of hours to schedule the job for
        
    Returns:
        None
        
    Creates:
        submitjob.sh in the parent directory of the cores"""
    
    # If only one core is passed, make it into a 1-length list for consistency
    if isinstance(cores, NemesisCore):
        cores = [cores,]

    script_path = cores[0].parent_directory + "submitjob.sh"

    # Read the submission script and replace template fields
    with open(constants.PATH+"data/statics/template.job", mode="r") as file:
        out = file.read()
        out = out.replace("<MEMORY>", str(memory))
        out = out.replace("<HOURS>", f"{hours:02}")
        out = out.replace("<N_CORES>", str(len(cores)))
        out = out.replace("<CORE_DIR>", os.path.abspath(cores[0].parent_directory))
        out = out.replace("<USERNAME>", username)

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