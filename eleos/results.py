import io
import itertools as it
import pandas as pd
import matplotlib.pyplot as plt
import glob
import re
from functools import wraps

from . import profiles
from . import utils
from . import cores


def plotting(func):
    @wraps(func)
    def wrapper(self, ax=None, *args, **kwargs):
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        else:
            fig = ax.get_figure()
        func(self, ax, *args, **kwargs)
        return fig, ax
    return wrapper


class NemesisResult:
    """Class for storing the results of a NEMESIS retrieval.
    
    Attributes:
        core_directory (str): The directory of the core being analysed
        ref (pandas.DataFrame): A DataFrame containing the data in the .ref file
        core (eleos.cores.NemesisCore): The NemesisCore object that generated the core directory
        profiles (list[eleos.porfiles.Profile]): A list of all the retrieved Profile objects from the run
        num_retrievals (int): I don't know what this param is - first field in the .mre file
        latitude (float): Latitude of the observed spectrum
        longitude (float): Longitude of the observed spectrum
        # Read in the fitted spectrum as a DataFrame
        self.fitted_spectrum = pd.read_table(mre_file, 
                                             names=["wavelength", "measured", "error", "pct_error", "model", "pct_diff"],
                                             index_col=0, sep="\s+", skiprows=5, nrows=blocks[0]-7)
"""
    def __init__(self, core_directory):
        """Inititalise a NemesisResult class
        
        Args:
            core_directory: The directory of a single core"""
        self.core_directory = core_directory
        self.ref = cores.parse_ref_file(self.core_directory+"nemesis.ref")
        self.core = cores.load_core(self.core_directory)
        self.profiles = self.core.profiles
        self._read_mre()
        self.chi_sq = self.get_chi_sq()

    def _parse_header_line(self, line, num_fields=-1, cast_to=float):
        """Take in a header line from the .mre file and split each token and cast to a dtype.
        Eg. "10     2     4     7.2" -> [10.0, 2.0, 4.0, 7.2]
        
        Args:
            line (str): The line to tokenise
            num_fields: The number of fields to include, starting from the left-hand side. Default is to include all fields
            cast_to: The type to cast the fields to.
        
        Returns:
            list[cast_to]: List of fields"""
        fields = [cast_to(x) for x in line.split()[:num_fields]]
        if num_fields == 1:
            return fields[0]
        else:
            return fields
                    
    def _read_mre(self):
        """Read in the nemesis.mre file"""
        mre_file = self.core_directory + "nemesis.mre"
        with open(mre_file) as file:
            mre_data = file.read().split("\n")

        header = []
        blocks = []
        for i, line in enumerate(mre_data):
            # Read the first 3 lines to the header array
            if i < 3:
                header.append(line)
            # Find the boundaries between result blocks
            elif "Variable" in line:
                blocks.append(i)
        blocks.append(i)

        # Set some attributes from the header info    
        self.num_retrievals = self._parse_header_line(header[0], num_fields=1, cast_to=int)
        self.ispec, self.ngeom, self.ny1, self.nx, self.ny2 = self._parse_header_line(header[1], num_fields=5, cast_to=int)
        self.latitude, self.longitude = self._parse_header_line(header[2], num_fields=2, cast_to=float)

        # Read in the fitted spectrum as a DataFrame
        self.fitted_spectrum = pd.read_table(mre_file, 
                                             names=["wavelength", "measured", "error", "pct_error", "model", "pct_diff"],
                                             index_col=0, sep="\s+", skiprows=5, nrows=blocks[0]-7)

        # Read in each retrieved parameter and add the result to the Shape object
        with open(mre_file) as file:
            for profile, (start, end) in zip(self.profiles, it.pairwise(blocks)):
                data = utils.read_between_lines(file, start, end)
                df = pd.read_table(io.StringIO(data), skiprows=4, sep="\s+", names=["i", "ix", "prior", "prior_error", "retrieved", "retrieved_error"])
                df.drop(["i", "ix"], axis=1, inplace=True)
                profile._add_result(df)

    def _format_time(self, decimal_hours):
        hours = int(decimal_hours)
        minutes = (decimal_hours*60) % 60
        seconds = (decimal_hours*3600) % 60
        return "%dh %02dm %02ds" % (hours, minutes, seconds)

    def get_chi_sq(self, all_iterations=False):
        """Get the chi squared values from the .prc file. If all_iterations is True, return a list containing chi squared values for all iterations"""
        if self.core.forward:
            pattern = "chisq/ny is equal to :    "
        else:
            pattern = "chisq/ny =    "
        with open(self.core_directory+"nemesis.prc") as file:
            lines = [line for line in file if pattern in line]
            values = [float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]) for line in lines]  # Extract all floats
        if all_iterations:
            return values
        else:
            return values[-1]

    def print_summary(self):
        print(f"Summary of retrieval in {self.core_directory}: ")
        print(f"Time taken: {self._format_time(self.elapsed_time)}")
        print(f"Chi squared value: {self.chi_sq}")
        
        for p in self.profiles:
           print(p)

    @property
    def elapsed_time(self):
        """Return the time taken for the retrieval in hours"""
        with open(self.core_directory+"nemesis.prc") as file:
            lines = file.read().splitlines()
            time = float(re.findall(r"[-+]?\d*\.\d+|\d+", lines[-1])[0])
        return time / 3600

    @plotting
    def plot_chi_sq(self, ax):
        print(self.core.num_iterations)
        ax.plot( self.get_chi_sq(all_iterations=True))
        ax.axhline(y=1, ls="dashed")
        ax.set_xlabel("Iteration Number")
        ax.set_ylabel("$\chi^2$")

    @plotting
    def plot_temperature(self, ax):
        """Plot the prior and retrieved temperature profile on a matplotlib Axes.
        
        Args:
            ax: The matplotlib.Axes object to plot to. If omitted then create a new Figure and Axes
            
        Returns:
            matplotlib.Figure: The Figure object to which the Axes belong
            matplotlib.Axes: The Axes object onto which the data was plotted"""
        # Find retrieved temperature profile
        for p in self.profiles:
            if isinstance(p, profiles.TemperatureProfile):
                temp_profile = p
                break
        else:
            raise AttributeError("No retrieved temperature profile found!")

        ax.plot(temp_profile.shape.data.retrieved, self.ref.height, c="k", lw=0.5, label="Retrieved")
        ax.plot(temp_profile.shape.data.prior, self.ref.height, c="r", lw=0.5, label="Prior")
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Height (km)")
        ax.legend()

    @plotting
    def plot_spectrum(self, ax):
        """Plot the measured and model spectrum on a matplotlib Axes.
        
        Args:
            ax: The matplotlib.Axes object to plot to. If omitted then create a new Figure and Axes
            
        Returns:
            matplotlib.Figure: The Figure object to which the Axes belong
            matplotlib.Axes: The Axes object onto which the data was plotted"""
        ax.plot(self.fitted_spectrum.wavelength, self.fitted_spectrum.measured, lw=0.5, label="Measured")
        ax.plot(self.fitted_spectrum.wavelength, self.fitted_spectrum.model, c="r", lw=0.5, label="Model")
        ax.set_yscale("log")
        ax.set_xlabel("Wavelength (μm)")
        ax.set_ylabel("Radiance\n(μW cm$^{-2}$ sr$^{-1}$ μm$^{-1}$)")
        ax.legend()
 

def load_multiple_cores(parent_directory):
    """Read in all the cores in a given directory and return a list of NemesisResult objects.
    
    Args:
        parent_directory (str): The directory containing all the individual core directories
        
    Returns:
        list[NemesisResult]: A list containing the result object for each core"""
    
    out = []
    for core in sorted(glob.glob(parent_directory+"core_*/")):
        out.append(NemesisResult(core))
    return out
