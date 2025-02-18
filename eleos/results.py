import pandas as pd
import matplotlib.pyplot as plt
import shutil
import re
import itertools
import numpy as np
from functools import wraps
from pathlib import Path

from . import profiles
from . import utils
from . import cores
from . import constants
from . import parsers


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


def plotting_altitude(func):
    @wraps(func)
    def wrapper(self, ax=None, pressure=True, *args, **kwargs):
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        else:
            fig = ax.get_figure()
        if pressure:
            ax.set_yscale("log")
            ax.set_ylabel("Pressure (bar)")
            ax.set_ylim(self.core.min_pressure, self.core.max_pressure)
            ax.invert_yaxis()
        else:
            ax.set_ylabel("Height (km)")
        func(self, ax, pressure, *args, **kwargs)
        return fig, ax
    return wrapper


class NemesisResult:
    """Class for storing and using the results of a NEMESIS retrieval.
    
    Attributes:
        core_directory (str):                   The directory of the core being analysed
        core (NemesisCore):                     The NemesisCore object that generated the core directory
        profiles (list[Profile]):               A list of all the retrieved Profile objects from the run
        latitude (float):                       Latitude of the observed spectrum
        longitude (float):                      Longitude of the observed spectrum
        chi_sq (float):                         The chi-squared value of the retrieval
        retrieved_spectrum (pandas.DataFrame):  A DataFrame containing the measured and modelled spectrum
        retrieved_aerosols (pandas.DataFrame):  A DataFrame containing the retrieved aerosol profiles
        retrieved_gases (pandas.DataFrame):     A DataFrame containing the retrieved chemical profiles
"""
    def __init__(self, core_directory):
        """Inititalise a NemesisResult class
        
        Args:
            core_directory: The directory of a single core"""
        self.core_directory = Path(core_directory)
        self.core = cores.load_core(self.core_directory)
        self.profiles = self.core.profiles
        self._read_mre()
        self.retrieved_aerosols = self._read_aerosol_prf()
        self.retrieved_gases = self._read_nemesis_prf()
        self.chi_sq = self.get_chi_sq()

    def _read_mre(self):
        mre = parsers.NemesisMre(self.core_directory / "nemesis.mre")
        self.__dict__ |= mre.__dict__
        for profile, df in zip(self.profiles, mre.retrieved_parameters):
            profile._add_result(df)

    def _read_aerosol_prf(self):
        header = ["height"] + [f"aerosol_{x}" for x in range(1, self.core.num_aerosol_modes+1)]
        data = pd.read_table(self.core.directory / "aerosol.prf", sep="\s+", skiprows=2, names=header)
        data.insert(1, "pressure", self.core.ref.data.pressure)
        return data

    def _read_itr(self):
        """Read in the nemesis.itr file. WARNING: This will almost certainly fail as it has file number offset hardcoded. It's on the todo list!
        
        Args:
            None
            
        Returns:
            pandas.DataFrame: A DataFrame containing the data from the .itr file with columns for each parameter"""
        
        # Reading the mre file to get the order of the parameters in the state vector
        with open(self.core_directory / "nemesis.mre") as file:
            lines = file.read().split("\n")
            toggle = False
            good = []
            for line in lines:
                if "i, ix" in line:
                    toggle = True
                elif "Variable" in line:
                    toggle = False
                if toggle:
                    good.append(line)
            great = []
            for line in good:
                if "i, ix" not in line:
                    great.append([float(x) for x in line.split()])
            del great[-1]
            vec = pd.DataFrame(great)
            vec.columns = ["i", "ix", "xa", "sa_err", "xn", "xn_err"]

        # Reading the itr file and extracting state vector for each iteration and the prior vector
        with open(self.core_directory / "nemesis.itr") as file:
            # Get the prior state vector
            prior = [float(x) for x in file.readlines()[5].split()]
            exps = []
            for p, v in zip(prior, vec.xa):
                exps.append(np.isclose(p, v))
            
            # Re-read the file and get the state vector for each iteration (which is the 3rd line after each blank line)
            file.seek(0)
            data = []
            count = -1
            for i, line in enumerate(file.read().split("\n")):
                if line == " ":
                    count = 3
                else:
                    count -= 1
                if count == 0:
                    d = []
                    for i, v in enumerate(line.split()):
                        d.append(float(v) if exps[i] else np.exp(float(v)))
                    data.append(d)

            # Assign the results to a DataFrame
            data = pd.DataFrame(data)
            state_vector_names = [f"{p.get_name()} {name.replace('_', ' ')}" for p in self.profiles for name in p.shape.NAMES]
            data.columns = state_vector_names

        return data

    def _read_nemesis_prf(self):
        with open(self.core_directory / "nemesis.prf") as file:
            lines = file.read().split("\n")
        
        names = []
        for i, line in enumerate(lines):
            # Skip first two lines and any blank lines
            if i in (0,1) or line == "":
                continue

            # Get the index of the starting line for the profiles
            if "height" in line:
                break
            
            # Otherwise it's a gas id
            else:
                gas_id, isotope_id = [int(x) for x in re.findall(r'\d+', line)]
                gas_name = constants.GASES[constants.GASES.radtrans_id == gas_id].name.iloc[0]
                names.append(f"{gas_name} {isotope_id}")

        out = pd.read_table(self.core_directory / "nemesis.prf", skiprows=i+1, sep="\s+")
        out.columns = ["height", "pressure", "temp"] + names

        return out

    def get_chi_sq(self, all_iterations=False):
        """Get the chi squared values from the .prc file. If all_iterations is True, return a list containing chi squared values for all iterations"""
        if self.core.forward:
            pattern = "chisq/ny is equal to :    "
        else:
            pattern = "chisq/ny =    "
        with open(self.core_directory / "nemesis.prc") as file:
            lines = [line for line in file if pattern in line]
            values = [float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]) for line in lines]  # Extract all floats
        if all_iterations:
            return values
        else:
            return values[-1]

    def print_summary(self):
        print(f"Summary of retrieval in {self.core_directory}")
        print(f"Time taken: {utils.format_decimal_hours(self.elapsed_time)}")
        print(f"Chi squared value: {self.chi_sq}")
        
        for p in self.profiles:
           print(p)

    @property
    def elapsed_time(self):
        """Return the time taken for the retrieval in hours"""
        with open(self.core_directory / "nemesis.prc") as file:
            lines = file.read().splitlines()
            time = float(re.findall(r"[-+]?\d*\.\d+|\d+", lines[-1])[0])
        return time / 3600

    @plotting
    def plot_chisq(self, ax):
        ax.plot( self.get_chi_sq(all_iterations=True))
        ax.axhline(y=1, ls="dashed")
        ax.set_xlabel("Iteration Number")
        ax.set_ylabel("$\chi^2$")

    @plotting
    def plot_spectrum(self, ax):
        """Plot the measured and model spectrum on a matplotlib Axes.
        
        Args:
            ax: The matplotlib.Axes object to plot to. If omitted then create a new Figure and Axes
            
        Returns:
            matplotlib.Figure: The Figure object to which the Axes belong
            matplotlib.Axes: The Axes object onto which the data was plotted"""
    
        ax.set_yscale("log")
        ax.plot(self.retrieved_spectrum.wavelength, self.retrieved_spectrum.measured, lw=0.5, label="Measured")
        ax.fill_between(self.retrieved_spectrum.wavelength, self.retrieved_spectrum.measured-self.retrieved_spectrum.error, self.retrieved_spectrum.measured+self.retrieved_spectrum.error, alpha=0.5)

        ax.plot(self.retrieved_spectrum.wavelength, self.retrieved_spectrum.model, c="r", lw=0.5, label="Model")

        plt.text(0.95, 0.05, f"$\chi^2 = ${self.chi_sq:.3f}",
            horizontalalignment='right',
            verticalalignment='bottom',
            transform = ax.transAxes)

        ax.set_xlabel("Wavelength (μm)")
        ax.set_ylabel("Radiance\n(μW cm$^{-2}$ sr$^{-1}$ μm$^{-1}$)")
        ax.legend()

    @plotting
    def plot_spectrum_residuals(self, ax):
        #ax.plot(self.retrieved_spectrum.wavelength, self.retrieved_spectrum.measured, lw=0.5, label="Measured")
        residuals = (self.retrieved_spectrum.model - self.retrieved_spectrum.measured) / self.retrieved_spectrum.measured

        ax.plot(self.retrieved_spectrum.wavelength, residuals, label="Residuals")
        ax.fill_between(self.retrieved_spectrum.wavelength, -self.retrieved_spectrum.error / self.retrieved_spectrum.measured, self.retrieved_spectrum.error / self.retrieved_spectrum.measured, alpha=0.5, label="Error")

        ax.set_xlabel("Wavelength (μm)")
        ax.set_ylabel("Residuals / Measurement")
        ax.legend()

    @plotting_altitude
    def plot_temperature(self, ax, pressure):
        """Plot the prior and retrieved temperature profile on a matplotlib Axes.
        
        Args:
            ax: The matplotlib.Axes object to plot to. If omitted then create a new Figure and Axes
            pressure: Whether to plot the temperature profile against pressure (if True) or height (if False)

        Returns:
            matplotlib.Figure: The Figure object to which the Axes belong
            matplotlib.Axes: The Axes object onto which the data was plotted"""

        print(self.core.ref)

        # Get the appropriate y axis data
        y = self.core.ref.pressure if pressure else self.core.ref.height
        
        # Find retrieved temperature profile
        temp_profile = None
        for p in self.profiles:
            if isinstance(p, profiles.TemperatureProfile):
                temp_profile = p
                ax.plot(temp_profile.shape.data.retrieved, y, c="k", lw=0.5, label="Retrieved")
                ax.plot(temp_profile.shape.data.prior, y, c="r", lw=0.5, label="Prior")
                ax.legend()
                break

        # If temperature profile not retrieved, plot the temperature profile in the .ref file
        if temp_profile is None:
            ax.plot(self.core.ref.temp, y)

        ax.set_xlabel("Temperature (K)")

    @plotting_altitude
    def plot_aerosol_profiles(self, ax, pressure, unit="tau/bar"):
        """Plot the retrieved aerosol profiles either in units of particles/gram of atmosphere or in units of optical thickness/bar
        against either height or pressure
        
        Args:
            ax: The matplotlib.Axes object to plot to. If omitted then create a new Figure and Axes
            unit: The unit to convert the aerosol profiles to. Valid values are 'tau/bar' for optical thickness/bar, 'cm2/g' for the native prf units
            pressure: Whether to plot the aerosol profiles against pressure (if True) or height (if False)
            
        Returns:
            matplotlib.Figure: The Figure object to which the Axes belong
            matplotlib.Axes: The Axes object onto which the data was plotted"""
        
        # Iterate over every retrieved aerosol
        for name in self.retrieved_aerosols.columns:
            if name in ("height", "pressure"):
                continue

            if unit == "tau/bar":
                # only god himself knows whats going on with these units...
                self.retrieved_aerosols[name] *= 1e5 / 10 / utils.get_planet_gravity(self.core.planet)
                unit_label = "Optical thickness / bar"
            elif unit == "cm2/g":
                unit_label = "Aerosol cross-section (cm$^2$ / g)"
            else:
                raise ValueError("Invalid unit! - Must be one of 'tau/bar', 'cm2/g'")

            # Get a label to use as the legend label - either the custom label or aerosol_<ID>
            id = int(name.removeprefix("aerosol_"))
            profile = self.core.get_aerosol_mode(id=id)
            if profile.label is not None:
                leg_label = profile.label
            else:
                leg_label = name

            y = self.retrieved_aerosols.pressure if pressure else self.retrieved_aerosols.height
            ax.plot(self.retrieved_aerosols[name], y, label=leg_label)

        ax.set_xlabel(unit_label)
        ax.legend()

    @plotting_altitude
    def plot_gas_profiles(self, ax, pressure, gas_names=None, include_priors=False):
        """Plot gas profiles from the .prf file.

        Args:
            ax (matplotlib.Axes): The matplotlib.Axes object to plot to. If omitted then create a new Figure and Axes.
            pressure (bool): Whether to plot the gas profiles against pressure (if True) or height (if False).
            gas_names (list[str], optional): List of gas names to plot. If None, plot all gases.
            include_priors (bool): Whether to plot the prior as well as the retrieved gas profiles
        
        Returns:
            matplotlib.Figure: The Figure object to which the Axes belong
            matplotlib.Axes: The Axes object onto which the data was plotted"""
        
        # Get the appropriate y axes
        y = self.retrieved_gases["pressure"] if pressure else self.retrieved_gases["height"]
        y2 = self.core.ref.data["pressure"] if pressure else self.core.ref.data["height"]
        
        # If no gas names specififed then get every gas profile
        if gas_names is None:
            gas_names = [x for x in self.retrieved_gases.columns if x not in ("height", "pressure", "temp")]
        
        # Plot the profiles
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for i, gas_name in enumerate(gas_names):
            ax.plot(self.retrieved_gases[gas_name], y, label=gas_name, c=colors[i])
            if include_priors:
                ax.plot(self.core.ref.data[gas_name], y2, label=f"{gas_name} Prior", c=colors[i])
        
        # Set a limit on lowest VMR
        x1, x2 = ax.get_xlim()
        if x1 < 1e-20:
            ax.set_xlim(1e-20, x2)

        ax.set_xscale("log")
        ax.set_xlabel("Volume Mixing Ratio")
        ax.legend()

    def make_summary_plot(self, figsize=(11, 10)):
        """Make a summary plot with prior and retrieved spectra, error on the spectra, aerosol and chemical profiles, and chi-squared values.
        
        Args:
            figsize (int, int): matplotlib figure size
            
        Returns:
            matplotlib.Figure: The produced figure
            dict(str: matplotlib.Axes): The axes of the produced figure"""
        names = []
        for profile in self.profiles:
            if isinstance(profile, profiles.GasProfile):
                names.append(f"{profile.gas_name} {profile.isotope_id}")

        fig, axs = plt.subplot_mosaic("AAA\nBBB\nCDE", 
                              gridspec_kw={"hspace": 0.25, "wspace": 0.35},
                              figsize=figsize)

        self.plot_spectrum(ax=axs["A"])
        self.plot_spectrum_residuals(ax=axs["B"])
        self.plot_chisq(ax=axs["C"])
        self.plot_aerosol_profiles(ax=axs["D"])
        self.plot_gas_profiles(ax=axs["E"], gas_names=names)

        fig.savefig(self.core_directory / "plots/summary.png", bbox_inches="tight", dpi=400)

        return fig, axs

    def make_iterations_plot(self, figsize=(14, 10)):
        """Plot the state vector for each iteration of the retrieval. Each retreived parameter is plotted on a separate axis.
        
        Args:
            None
            
        Returns:
            matplotlib.Figure: The Figure object to which the Axes belong
            matplotlib.Axes: The Axes object onto which the data was plotted"""
                
        data = self._read_itr()
            
        nrows = int(np.ceil(np.sqrt(len(data.columns))))
        ncols = int(np.ceil(len(data.columns) / nrows))

        fig, axs = plt.subplots(nrows, ncols, figsize=figsize, sharex=True)
        for ax, param in itertools.zip_longest(axs.reshape(-1), data):
            if param is None:
                ax.axis("off")
                continue
            ax.plot(data[param])
            ax.set_ylabel(param)
        
        fig.supxlabel('Iteration Number')
        fig.tight_layout()
        fig.savefig(self.core_directory / "plots/iterations.png", bbox_inches="tight", dpi=400)
        return fig, axs

    def savefig(self, name, fig=None, **kwargs):
        """
        Save a matplotlib figure to a file in the core's `plots` directory.

        Args:
            name (str): The name of the file to save the figure as.
            fig (matplotlib.Figure, optional): The figure to save. If None, the current figure will be saved. Default is None.
            **kwargs: Additional keyword arguments to pass to `savefig`.

        Returns:
            None
        """
        x = plt if fig is None else fig
        x.savefig(self.core_directory / "plots/" + name, bbox_inches="tight", **kwargs)

    def delete(self, confirm=True):
        """Delete the NemesisResult object AND delete the corresponding core directory. This action is irreviersible!
        
        Args:
            confirm (bool): Whether to prompt for confirmation
            
        Returns:
            None"""
        
        if confirm:
            answer = input(f"Are you sure you want to delete this core ({self.core_directory})? y/n ").lower()
            if answer != "y":
                return
            
        shutil.rmtree(self.core_directory)
        del self


def load_multiple_cores(parent_directory, raise_errors=True):
    """Read in all the cores in a given directory and return a list of NemesisResult objects.
    
    Args:
        parent_directory (str): The directory containing all the individual core directories
        raise_errors (bool): Whether to raise an error if a retieval failed (True) or to silently skip it (False) 
        
    Returns:
        list[NemesisResult]: A list containing the result object for each core"""
    
    parent_directory = Path(parent_directory)

    out = []
    for core in sorted(parent_directory.glob("core_*")):
        try:
            out.append(NemesisResult(core))
        except Exception as e:
            if raise_errors:
                raise e
    return out


