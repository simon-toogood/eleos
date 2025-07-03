import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import shutil
import re
import itertools
import numpy as np
from functools import wraps
from pathlib import Path

import warnings
warnings.formatwarning = lambda msg, *_: f"Warning: {msg}\n"

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
        profiles (dict[Profile]):               A dictionary of all the retrieved Profile objects from the run. The keys are the labels given to the Profiles on creation (eg. GasProfiles have the form "<gas_name> <isotope_id>" such as "PH3 0")
        latitude (float):                       Latitude of the observed spectrum
        longitude (float):                      Longitude of the observed spectrum
        chi_sq (float):                         The chi-squared value of the retrieval
        retrieved_spectrum (pandas.DataFrame):  A DataFrame containing the measured and modelled spectrum
        retrieved_aerosols (pandas.DataFrame):  A DataFrame containing the retrieved aerosol profiles
        retrieved_gases (pandas.DataFrame):     A DataFrame containing the retrieved chemical profiles
    """

    def __init__(self, core_directory): 
        """Constructor for NemesisResult
        
        Args:
            core_directory (str): The directory of the core
        """   
        # Load core directory
        self.core_directory = Path(core_directory)
        self.core = cores.load_core(self.core_directory)
        self.profiles = self.core.profiles

        # Parse some files
        self.mre = parsers.NemesisMre(self.core_directory / "nemesis.mre")
        self.aerosol_prf = parsers.AerosolPrf(self.core_directory / "aerosol.prf")
        self.nemesis_prf = parsers.NemesisPrf(self.core_directory / "nemesis.prf")
        self.aerosol_prf.data["pressure"] = self.nemesis_prf.data["pressure"]

        # Parse the iterations file for retrievals
        if not self.core.forward:
            self.itr = parsers.NemesisItr(self.core_directory / "nemesis.itr")
            self.itr.add_column_names(self.profiles)

        # Add the mre object attributes to the NemesisResult object for convenience
        self.__dict__ |= self.mre.__dict__
        self.retrieved_aerosols = self.aerosol_prf.data
        self.retrieved_gases = self.nemesis_prf.data
        self.chi_sqs = self._get_chi_squareds()
        self.chi_sq = self.chi_sqs[-1]

        # Add results to the profiles
        self._add_results_to_profiles()

        # Set some misc params
        self._time = None

    def _add_results_to_profiles(self):
        for label, profile in self.profiles.items():
            if isinstance(profile, profiles.AerosolProfile) and profile.retrieve_optical:
                profile._add_result(self.mre.retrieved_parameters.pop(0),
                                    self.mre.retrieved_parameters.pop(0))
            else:
                profile._add_result(self.mre.retrieved_parameters.pop(0))

    def _get_chi_squareds(self):
        """Read the .prc file and cache the chi squared values to a new .chi file if not already done.
        Otherwise, read the .chi file and return the chi squared values
        
        Args:
            None
        
        Returns:
            None
            
        Creates:
            nemesis.chi"""

        if (self.core_directory / "nemesis.chi").is_file():
            with open(self.core_directory / "nemesis.chi") as file:
                return [float(line) for line in file]
        else:
            # Pattern to match to find values in the .prc file
            if self.core.forward:
                pattern = "chisq/ny is equal to :"
            else:
                pattern = "chisq/ny = "
            
            # Parse the .prc file
            with open(self.core_directory / "nemesis.prc") as file:
                lines = [line for line in file if pattern in line]
                try:
                    values = [float(re.findall(r"[-+]?\d*\.\d+|\d+", line)[0]) for line in lines]  # Extract all floats
                except IndexError:
                    values = np.array([np.nan for _ in lines])
            # Write results to .chi file
            with open(self.core_directory / "nemesis.chi", "w+") as file:
                file.write("\n".join(map(str, values)))

            if len(values) == 0:
                values = [np.nan]

            return values

    def print_summary(self, colors=False):
        print(f"Summary of retrieval in {self.core_directory}")
        print(f"Time taken: {utils.format_decimal_hours(self.elapsed_time)}")
        print(f"Chi squared value: {self.chi_sq}")
        
        for name, profile in self.profiles.items():
           profile.print_table(colors=colors)

    @property
    def elapsed_time(self):
        """Return the time taken for the retrieval in hours"""
        if self._time is None:
            with open(self.core_directory / "nemesis.prc") as file:
                lines = file.read().splitlines()
                time = float(re.findall(r"[-+]?\d*\.\d+|\d+", lines[-1])[0])
            return time / 3600
        else:
            return self._time

    @plotting
    def plot_chisq(self, ax):
        if self.core.forward:
            ax.scatter(0, self.chi_sq)
        else:
            ax.plot(self.chi_sqs)
        ax.axhline(y=1, ls="dashed")
        ax.set_xlabel("Iteration Number")
        ax.set_ylabel("$\chi^2$")

    @plotting
    def plot_spectrum(self, ax, show_chisq=True, legend=True):
        """Plot the measured and model spectrum on a matplotlib Axes.
        
        Args:
            ax: The matplotlib.Axes object to plot to. If omitted then create a new Figure and Axes
            show_chisq (bool): Whether to display the chi-squared value of the fit
            legend (bool): Whether to draw the legend

        Returns:
            matplotlib.Figure: The Figure object to which the Axes belong
            matplotlib.Axes: The Axes object onto which the data was plotted"""
    
        ax.set_yscale("log")
        ax.plot(self.retrieved_spectrum.wavelength, self.retrieved_spectrum.measured, lw=0.5, label="Measured" if legend else None)
        ax.fill_between(self.retrieved_spectrum.wavelength, self.retrieved_spectrum.measured-self.retrieved_spectrum.error, self.retrieved_spectrum.measured+self.retrieved_spectrum.error, alpha=0.5)

        ax.plot(self.retrieved_spectrum.wavelength, self.retrieved_spectrum.model, c="r", lw=0.5, label="Model" if legend else None)

        if show_chisq:
            plt.text(0.95, 0.05, f"$\chi^2 = ${self.chi_sq:.3f}",
                horizontalalignment='right',
                verticalalignment='bottom',
                transform = ax.transAxes)

        ax.set_xlabel("Wavelength (μm)")
        ax.set_ylabel("Radiance\n(μW cm$^{-2}$ sr$^{-1}$ μm$^{-1}$)")
        if legend:
            ax.legend()

    @plotting
    def plot_spectrum_residuals(self, ax):
        residuals = np.log(self.retrieved_spectrum.model) - np.log(self.retrieved_spectrum.measured)

        ax.plot(self.retrieved_spectrum.wavelength, residuals, label="Residuals", lw=1)
        ax.fill_between(self.retrieved_spectrum.wavelength,
                         -self.retrieved_spectrum.error / self.retrieved_spectrum.measured, 
                         self.retrieved_spectrum.error / self.retrieved_spectrum.measured, 
                         alpha=0.5, label="Error")
        ax.axhline(y=0, zorder=-1, c="k", ls="dashed")
        ax.set_xlabel("Wavelength (μm)")
        ax.set_ylabel("Log Residuals")
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
        max_value = -1

        if self.core.num_aerosol_modes == 0:
            return
        
        for label in self.retrieved_aerosols.columns:
            if label in ("height", "pressure"):
                continue

            x = self.retrieved_aerosols[label]
            y = self.retrieved_aerosols.pressure if pressure else self.retrieved_aerosols.height

            if unit == "tau/bar":
                # only god himself knows whats going on with these units...
                x *= 1e5 / 10 / constants.GRAVITY[self.core.planet]
                unit_label = f"Optical depth / bar at {self.core.reference_wavelength:.1f}µm"
            elif unit == "particles/g":
                unit_label = f"Aerosol specific density (particles / gram) at {self.core.reference_wavelength:.1f}µm)"
            else:
                raise ValueError(f"Invalid unit! - Must be one of 'tau/bar', 'particles/g' - not {unit}")
            
            if x.max() > max_value:
                max_value = x.max()

            ax.plot(x, y, label=label)

        ax.set_xlabel(unit_label)
        ax.set_xscale("log")
        ax.set_xlim(1e-6, max_value*2)
        ax.legend()

    @plotting_altitude
    def plot_gas_profiles(self, ax, 
                          pressure, 
                          unit="", 
                          gas_names=None,
                          plot_retrieved_profiles=True,
                          plot_prior_profiles=False, 
                          plot_ref_profiles=True):
        
        """Plot gas profiles from the .prf file.

        Args:
            ax (matplotlib.Axes): The matplotlib.Axes object to plot to. If omitted then create a new Figure and Axes.
            pressure (bool): Whether to plot the gas profiles against pressure (if True) or height (if False).
            unit (str): One of '', 'ppm', 'ppb' or 'ppt'. 
            gas_names (list[str], optional): List of gas names to plot. If None, plot all gases.
            plot_retrieved_profiles (bool): Whether to plot the retrieved gas profiles
            plot_prior_profiles (bool): Whether to plot the prior gas profiles
            plot_ref_profiles (bool): Whether to plot the profile in the .ref file

        Returns:
            matplotlib.Figure: The Figure object to which the Axes belong
            matplotlib.Axes: The Axes object onto which the data was plotted"""
        
        # If no gas names specififed then get every gas profile
        if gas_names is None:
            gas_names = [x for x in self.retrieved_gases.columns if x not in ("height", "pressure", "temperature")]
        # Allow passing in of a single string instead of a list
        elif isinstance(gas_names, str):
            gas_names = [gas_names]

        # Get the appropriate y axes
        y = self.retrieved_gases["pressure"] if pressure else self.retrieved_gases["height"]
        y2 = self.core.ref.data["pressure"] if pressure else self.core.ref.data["height"]

        # Get the prior distributions if requested
        if plot_prior_profiles:
            priors = self.core.generate_prior_distributions()
            y3 = priors["pressure"] if pressure else priors["height"]

        # Determine the scaling factor
        if unit == "":
            scale = 1
        elif unit == "ppm":
            scale = 1e6
        elif unit == "ppb":
            scale = 1e9
        elif unit == "ppt":
            scale = 1e12
        else:
            raise ValueError("Invalid unit! - Must be one of '', 'ppm', 'ppb', 'ppt'")
        
        # Plot the profiles
        colors = itertools.cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
        for i, gas_name in enumerate(gas_names):
            c = next(colors)
            if plot_retrieved_profiles:
                ax.plot(self.retrieved_gases[gas_name]*scale, y, label=gas_name, c=c)
            if plot_ref_profiles:
                ax.plot(self.core.ref.data[gas_name]*scale, y2, label=f"{gas_name} Reference", c=c, ls="dashed")
            if plot_prior_profiles:
                ax.plot(priors[gas_name]*scale, y3, label=f"{gas_name} Prior", c=c, ls="dotted")

        # Set a limit on lowest VMR
        x1, x2 = ax.get_xlim()
        if x1*scale < 1e-12:
            ax.set_xlim(1e-12*scale, 1*scale)

        # Set unit label
        if unit == "":
            label = unit
        else:
            label = f"({unit})"
        ax.set_xscale("log")
        ax.set_xlabel(f"Volume Mixing Ratio {label}")
        ax.legend()

    def make_summary_plot(self, figsize=(11, 10)):
        """Make a summary plot with prior and retrieved spectra, error on the spectra, aerosol and chemical profiles, and chi-squared values.
        
        Args:
            figsize (int, int): matplotlib figure size
            
        Returns:
            matplotlib.Figure: The produced figure
            dict(str: matplotlib.Axes): The axes of the produced figure with labels:
                'A' for the spectrum
                'B' for the residuals
                'C' for the chi-sqaured plot
                'D' for the aerosol profiles
                'E' for the gas profiles"""
        
        names = []
        for label, profile in self.profiles.items():
            if isinstance(profile, profiles.GasProfile):
                names.append(label)

        fig, axs = plt.subplot_mosaic("AAA\nBBB\nCDE", 
                              gridspec_kw={"hspace": 0.25, "wspace": 0.35},
                              figsize=figsize)

        self.plot_spectrum(ax=axs["A"])
        self.plot_spectrum_residuals(ax=axs["B"])
        self.plot_chisq(ax=axs["C"])
        self.plot_aerosol_profiles(ax=axs["D"])
        self.plot_gas_profiles(ax=axs["E"], gas_names=names)

        if self.core.forward:
            fig.suptitle(f"Forward model in {self.core_directory.resolve()}", y=0.91)
        else:
            fig.suptitle(f"Retrieval in {self.core_directory.resolve()}", y=0.91)

        fig.savefig(self.core_directory / "plots/summary.png", bbox_inches="tight", dpi=400)

        return fig, axs

    def make_iterations_plot(self, figsize=(14, 10)):
        """Plot the state vector for each iteration of the retrieval. Each retreived parameter is plotted on a separate axis.
        
        Args:
            None
            
        Returns:
            matplotlib.Figure: The Figure object to which the Axes belong
            matplotlib.Axes: The Axes object onto which the data was plotted"""

        if self.core.forward:
            raise TypeError("Forward models cannot have interations plotted")

        data = self.itr.state_vectors

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
        x.savefig(self.core_directory / "plots" / name, bbox_inches="tight", **kwargs)

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


class SensitivityAnalysis:
    """Class for analysing the sensitivity cores generated by `eleos.cores.create_sensitivity_analysis`
    
    Attributes:
        parent_directory (str): The directory containing the sensitivity cores
        results (list[NemesisResult]): A list of NemesisResult objects for each core
        baseline (NemesisResult): The baseline result, an alias for results[0]
        params (pandas.DataFrame): A DataFrame  detailing which parameters were varied in each core and their value
    """

    def __init__(self, parent_directory):
        self.parent_directory = Path(parent_directory)
        self.results = load_multiple_cores(parent_directory)
        self.baseline = self.results[0]
        self.params = pd.read_csv(self.parent_directory / "sensitivity_analysis.txt")
    
    def _get_all_params(self):
        out = []
        prev = (None, None)
        for i, row in self.params.iterrows():
            x = (row["Profile Label"], row["Parameter"])
            if x != prev:
                out.append((row["Profile Label"], row["Parameter"]))
                prev = x
        return out

    def _get_params(self, profile_label, parameter):
        """Filter the params DataFrame to get the rows that match the given profile label and parameter."""
        return self.params[(self.params["Parameter"] == parameter) & (self.params["Profile Label"] == profile_label)]

    def get_results(self, profile_label, parameter):
        """Get the cores that varied the given parameter in the given profile.
        
        Args:
            profile_label: The label of the profile
            parameter: The parameter that was varied
            
        Returns:
            list[NemesisResult]: A list of NemesisResult objects that varied the given parameter in the given profile"""
        
        out = []
        for _, case in self._get_params(profile_label, parameter).iterrows():
            out.append(self.results[case["Core ID"] - 1])
        return out
        
    @plotting
    def plot_parameter(self, ax, profile_label, parameter):
        """Plot the sensitivity of the model to the given parameter in the given profile.
        
        Args:
            ax: The matplotlib.Axes object to plot to. If omitted then create a new Figure and Axes
            profile_label: The label of the profile to plot
            parameter: The variable to plot
            
        Returns:
            matplotlib.Figure: The Figure object to which the Axes belong
            matplotlib.Axes: The Axes object onto which the data was plotted"""
        
        def alpha_map(x, V, min_alpha):
            return 1 - (np.abs(x-1)/V) * (1 - min_alpha)

        df = self._get_params(profile_label, parameter)
        ress = self.get_results(profile_label, parameter)
        base = self.baseline.retrieved_spectrum.model

        for factor, r in zip(df["Factor"], ress):
            y = r.retrieved_spectrum.model / base
            ax.plot(r.retrieved_spectrum.wavelength, 
                    y, 
                    alpha=alpha_map(factor, 1-df["Factor"].min(), 0.25), 
                    color="#FF0000" if factor > 1 else "#0044FF",
                    label=factor)
            
        
        low, high = ax.get_ylim()
        bound = np.max((np.abs(low-1), np.abs(high-1)))
        ax.set_ylim(-bound+1, bound+1)
        ax.set_ylabel(f"Change from baseline")
        ax.set_xlabel("Wavelength (µm)")
        ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))
        ax.axhline(1, c="k", ls="dashed")

    def make_parameters_plot(self, ncol=3):
        p = self._get_all_params()
        nrow = int(np.ceil(len(p) / ncol))

        fig, axs = plt.subplots(nrow , ncol, figsize=(4*ncol, 1.5*nrow), sharex=True)
        axs = axs.flatten()

        for ax in axs[-(nrow*ncol - len(p)):]:
            ax.set_axis_off()

        for ax, name in zip(axs, p):
            self.plot_parameter(ax, *name)
            ax.set_ylabel("")
            ax.set_xlabel("")
            ax.text(0.99, 0.05, " ".join(name).replace("_", " "), transform=ax.transAxes, ha="right", va="bottom")

        fig.supxlabel("Wavelength (µm)")
        fig.supylabel("Radiance change from baseline", x=0.01)
        fig.tight_layout()
        fig.savefig("plots/sensitivity.png", dpi=300)
    
    def savefig(self, name, fig=None, **kwargs):
        """
        Save a matplotlib figure to a file in the parent_directory.

        Args:
            name (str): The name of the file to save the figure as.
            fig (matplotlib.Figure, optional): The figure to save. If None, the current figure will be saved. Default is None.
            **kwargs: Additional keyword arguments to pass to `savefig`.

        Returns:
            None
        """
        x = plt if fig is None else fig
        x.savefig(self.parent_directory / name, bbox_inches="tight", **kwargs)


def load_multiple_cores(parent_directory, raise_errors=True):
    """Read in all the cores in a given directory and return a list of NemesisResult objects.
    
    Args:
        parent_directory (str): The directory containing all the individual core directories
        raise_errors (bool): Whether to raise an error if a retieval failed (True) or to silently skip it (False) 
        
    Returns:
        list[NemesisResult]: A list containing the result object for each core"""
    
    parent_directory = Path(parent_directory)

    out = []
    for core in sorted(parent_directory.glob("core_*"), key=sort_key_paths):
        try:
            out.append(NemesisResult(core))
        except Exception as e:
            warnings.warn(f"Failed to load core {core}")
            if raise_errors:
                raise e
    return out


def load_best_cores(parent_directory, n, raise_errors=True):
    """Load n cores from the parent_directory with the lowest chi-squared values.

    Args:
        parent_directory (str): The directory containing all the individual core directories
        n (int): The number of cores to load
    
    Returns:
        list[NemesisResult]: A list containing the result object for each core, sorted by chi-squared value (so lowest chi-sqared is index 0 in the list)"""
    parent_directory = Path(parent_directory)

    dirs = [""]
    chisqs = [float("inf")]
    for core_directory in parent_directory.glob("core_*"):
        prc = parsers.NemesisPrc(core_directory / "nemesis.prc")
        try:
            chisq = prc.chisq[-1]
        except IndexError:
            continue
        if len(chisqs) < n:
            dirs.append(core_directory)
            chisqs.append(chisq)
        elif chisq < max(chisqs):
            i = np.argmax(chisqs)
            dirs[i] = core_directory
            chisqs[i] = chisq


    chis, dirs = zip(*sorted(zip(chisqs, dirs)))

    rs = []
    for d in dirs:
        try:
            rs.append(NemesisResult(d))
        except Exception as e:
            if raise_errors:
                raise e
            
    return rs


def sort_key_paths(path):
    x = str(path).split("_")[-1]
    return int(x)