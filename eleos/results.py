import shutil
import re
import itertools
from functools import wraps
from pathlib import Path
import pickle
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmasher
from cycler import cycler
from astropy import units as u
from tqdm import tqdm

from . import profiles
from . import utils
from . import cores
from . import constants
from . import parsers
from . import shapes
from . import mcmc

# set up warning formatter
warnings.formatwarning = lambda msg, *_: f"Warning: {msg}\n"

# set up custom mpl cycler
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
linestyles = ['-', '--', ':', '-.']
style_cycler = cycler(linestyle=linestyles) * cycler(color=colors)
plt.rcParams['axes.prop_cycle'] = style_cycler


class NemesisResult:
    """Class for storing and using the results of a NEMESIS retrieval.
    
    Attributes:
        core_directory (str):                   The directory of the core being analysed
        core (NemesisCore):                     The NemesisCore object that generated the core directory
        profiles (dict[Profile]):               A dictionary of all the retrieved Profile objects from the run. The keys are the labels given to the Profiles on creation (eg. GasProfiles have the form "<gas_name> <isotope_id>" such as "PH3 0")
        latitude (float):                       Latitude of the observed spectrum
        longitude (float):                      Longitude of the observed spectrum
        chi_sq (float):                         The chi-squared value of the retrieval
        elapsed_time (float):                   The time taken for the retrieval in decimal hours
        retrieved_spectrum (pandas.DataFrame):  A DataFrame containing the measured and modelled spectrum
        retrieved_aerosols (pandas.DataFrame):  A DataFrame containing the retrieved aerosol profiles
        retrieved_gases (pandas.DataFrame):     A DataFrame containing the retrieved chemical profiles
    """

    def __new__(cls, core_directory=None):
        """Constructor for NemesisResult. Dont worry about why this is a __new__ method and not a __init__ method.
        
        Args:
            core_directory (str): The directory of the core
        """   

        # Check for cached file (this speeds up loading by skipping all the expensive file parsing ops)
        if core_directory is None:
            return super().__new__(cls)
        if (Path(core_directory) / "results.pkl").exists():
            with open(Path(core_directory) / "results.pkl", mode="rb") as file:
                return pickle.load(file)
        
        # Otherwise, load core directory
        self = super().__new__(cls)
        self.core_directory = Path(core_directory)
        self.core = cores.load_core(self.core_directory)
        self.core.spectrum = parsers.NemesisSpx(self.core_directory / "nemesis.spx")
        self.profiles = self.core.profiles

        # Parse some files
        self.mre = parsers.NemesisMre(self.core_directory / "nemesis.mre")
        self.aerosol_prf = parsers.AerosolRef(self.core_directory / "aerosol.prf")
        self.nemesis_prf = parsers.NemesisRef(self.core_directory / "nemesis.prf")
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

        # Pickle the object
        self.save()

        return self

    def _add_results_to_profiles(self):
        for label, profile in self.profiles.items():
            if isinstance(profile.shape, shapes.Shape0):
                print("Warning, Shape0 is not fully supported (from NemesisResult._add_results_to_profiles)")
                continue
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

        prc = parsers.NemesisPrc(self.core_directory / "nemesis.prc")
        prc.write_chisqs(self.core_directory / "nemesis.chi")
        return prc.chisq

    def _convert_aerosol_units(self, data, p_ref=1):
        """
        Convert the aerosol densities in the aerosol.prf file to more useful units. This assumes that the cross-sections are normalised at a specific wavelength, such that the native units are cm2/g.
        This function converts them to optical depth per km, optical depth per bar, and pressure-integrated optical depth.

        Args:
            data (np.array): Aerosol densities in cm2/g
            p_ref (float):   The reference pressure in bar at which to calculate everything. Default is 1 bar and should not need to be changed.

        Returns:
            opacity_ODkm (np.array):  Aerosol densities in optical depth per km
            opacity_ODbar (np.array): Aerosol densities in optical depth per bar
            opacity_ODkm_int (np.array): Aerosol densities in optical depth integrated over height
            opacity_ODbar_int (np.array): Aerosol densities in optical depth integrated over pressure

        """

        # first, calculate the average molecular weight of the atmosphere
        masses = []
        for gas in self.core.reference.data.columns:
            if gas in ("height", "pressure", "temperature"):
                continue
            g = gas.split(" ")[0] # ignore isotope weights
            m = constants.GASES.loc[constants.GASES["name"] == g, "molar_mass"].iloc[0]
            masses.append(m)
        vmrs = self.core.reference.data.copy()
        x = vmrs.drop(columns=["height", "pressure", "temperature"])
        vmrs["mean_molar_wt"] = (x * masses).sum(axis=1)

        # then calculate the values at the target level
        target = vmrs.loc[utils.find_nearest(self.core.reference.data["pressure"], p_ref)[0]]
        pres = target["pressure"] * u.bar
        temp = target["temperature"] * u.K
        mol_wt = target["mean_molar_wt"] * u.g / u.mol
        gravity = constants.GRAVITY[self.core.planet] * u.m / u.s**2
        R = 8.314 * u.J / u.K / u.mol

        Hscale = ((R * temp) / (mol_wt * gravity)) 
        # print("Atmosphere scale height: ", Hscale.decompose())

        rho = ((mol_wt * pres) / (R * temp))
        # print("Atmospheric density: ", rho.decompose())

        opacity_ODkm = list(data) * u.cm**2/u.g * rho
        opacity_ODkm = opacity_ODkm.to(1 / u.km)
        # print("Maximum optical depth per km: ", max(opacity_ODkm))

        opacity_ODbar = opacity_ODkm * Hscale / pres
        opacity_ODbar = opacity_ODbar.to(1 / u.bar)
        # print("Maximum optical depth per bar: ", max(opacity_ODbar))

        opacity_ODbar_int = np.cumsum(opacity_ODbar)

        return opacity_ODkm.value, opacity_ODbar.value, opacity_ODbar_int.value

    def save(self):
        """Pickle the object for quicker loading later
        """
        with open(self.core_directory / "results.pkl", mode="wb") as file:
            pickle.dump(self, file)

    def print_summary(self, colors=False):
        print(f"Summary of retrieval in {self.core_directory}")
        print()
        print(f"Time taken: {utils.format_decimal_hours(self.elapsed_time)}")
        print(f"Number of iterations: {len(self.chi_sqs)}")
        print("Time per iteration: " + utils.format_decimal_hours(self.elapsed_time / len(self.chi_sqs)))
        print(f"Number of retrieved parameters: {len(self.mre.initial_state_vector)}")
        print(f"Chi squared value: {self.chi_sq}")
        print()
        
        for name, profile in self.profiles.items():
           profile.print_table(colors=colors, forward=self.core.forward)

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

    @utils.plotting
    def plot_chisq(self, ax):
        """Plot the chi-squared values as a function of iteration number
        
        Args:
            ax: The matplotlib.Axes object to plot to. If omitted then create a new Figure and Axes

        Returns:
            matplotlib.Figure: The Figure object to which the Axes belong
            matplotlib.Axes: The Axes object onto which the data was plotted
        """
        if self.core.forward:
            ax.scatter(0, self.chi_sq)
        else:
            ax.plot(self.chi_sqs)
        ax.axhline(y=1, ls="dashed")
        ax.set_xlabel("Iteration Number")
        ax.set_ylabel(R"$\chi^2$")

    @utils.plotting
    def plot_spectrum(self, ax, show_chisq=True, legend=True, log=False):
        """Plot the measured and model spectrum on a matplotlib Axes.
        
        Args:
            ax: The matplotlib.Axes object to plot to. If omitted then create a new Figure and Axes
            show_chisq (bool): Whether to display the chi-squared value of the fit
            legend (bool): Whether to draw the legend
            log (bool): Whether to use a log plot

        Returns:
            matplotlib.Figure: The Figure object to which the Axes belong
            matplotlib.Axes: The Axes object onto which the data was plotted"""
    
        if log:
            ax.set_yscale("log")

        ax.plot(self.retrieved_spectrum.wavelength, self.retrieved_spectrum.measured, lw=0.5, label="Measured" if legend else None)
        ax.fill_between(self.retrieved_spectrum.wavelength, self.retrieved_spectrum.measured-self.retrieved_spectrum.error, self.retrieved_spectrum.measured+self.retrieved_spectrum.error, alpha=0.5)

        ax.plot(self.retrieved_spectrum.wavelength, self.retrieved_spectrum.model, c="r", lw=0.5, label="Model" if legend else None)

        ax.set_xlabel("Wavelength (μm)")
        ax.set_ylabel("Radiance\n(μW cm$^{-2}$ sr$^{-1}$ μm$^{-1}$)")

        if legend:
            if self.chi_sq is not None and show_chisq:
                title = rf"$\chi^2 = {self.chi_sq:.3f}$"
            else:
                title = ""
            leg = ax.legend(title=title)
            for line in leg.get_lines():
                line.set_linewidth(1)

    @utils.plotting
    def plot_spectrum_residuals(self, ax, log=False):
        """Plot the spectrum residuals on a matplotlib Axes.
        
        Args:
            ax: The matplotlib.Axes object to plot to. If omitted then create a new Figure and Axes
            show_chisq (bool): Whether to display the chi-squared value of the fit
            log (bool): Whether to use a log plot

        Returns:
            matplotlib.Figure: The Figure object to which the Axes belong
            matplotlib.Axes: The Axes object onto which the data was plotted"""
        
        if log:
            residuals = np.log(self.retrieved_spectrum.model) - np.log(self.retrieved_spectrum.measured)

            ax.plot(self.retrieved_spectrum.wavelength, residuals, label="Residuals", lw=1)
            ax.fill_between(self.retrieved_spectrum.wavelength,
                            -self.retrieved_spectrum.error / self.retrieved_spectrum.measured, 
                            self.retrieved_spectrum.error / self.retrieved_spectrum.measured, 
                            alpha=0.5, label="Error")
            ax.set_ylabel("Log Residuals")

        else:
            residuals = self.retrieved_spectrum.model - self.retrieved_spectrum.measured

            ax.plot(self.retrieved_spectrum.wavelength, residuals, label="Residuals", lw=1)
            ax.fill_between(self.retrieved_spectrum.wavelength,
                            residuals-self.retrieved_spectrum.error, 
                            residuals+self.retrieved_spectrum.error, 
                            alpha=0.5, label="Error")
            ax.set_ylabel("Residuals\n(μW cm$^{-2}$ sr$^{-1}$ μm$^{-1}$)")

        ax.axhline(y=0, zorder=-1, c="k", ls="dashed")
        ax.set_xlabel("Wavelength (μm)")
        ax.legend()

    @utils.plotting_altitude
    def plot_temperature(self, ax, pressure):
        """Plot the prior and retrieved temperature profile on a matplotlib Axes.
        
        Args:
            ax: The matplotlib.Axes object to plot to. If omitted then create a new Figure and Axes
            pressure: Whether to plot the temperature profile against pressure (if True) or height (if False)

        Returns:
            matplotlib.Figure: The Figure object to which the Axes belong
            matplotlib.Axes: The Axes object onto which the data was plotted"""

        # Get the appropriate y axis data
        y = self.core.reference.pressure if pressure else self.core.reference.height
        
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
            ax.plot(self.core.reference.temp, y)

        ax.set_xlabel("Temperature (K)")

    @utils.plotting_altitude
    def plot_aerosol_profiles(self, 
                              ax,
                              pressure=True, 
                              unit="tau/bar", 
                              style="line", 
                              plot_prior_profiles=True,
                              hide_extra_labels=True):
        """Plot the retrieved aerosol profiles either in units of particles/gram of atmosphere or in units of optical thickness/bar
        against either height or pressure
        
        Args:
            ax: The matplotlib.Axes object to plot to. If omitted then create a new Figure and Axes
            unit: The unit to convert the aerosol profiles to. Valid values are:
                  'tau/bar' for optical depth per bar, 
                  'tau/km' for optical depth per bar, 
                  'tau' for pressure-integrated optical depth, 
                  'cm2/g' for the native prf units
            pressure: Whether to plot the aerosol profiles against pressure (if True) or height (if False)
            style: One of 'line' or 'step'
            plot_prior_profiles: Whether to plot the prior aerosol profiles
            hide_extra_labels (bool): Whether to show labels for eveything in the legend (False) or only the species colour (True)

        Returns:
            matplotlib.Figure: The Figure object to which the Axes belong
            matplotlib.Axes: The Axes object onto which the data was plotted"""
        
        if self.core.num_aerosol_modes == 0:
                return

        if unit not in {"tau/bar", "tau/km", "tau", "cm2/g"}:
            raise ValueError(
                f"Invalid unit {unit!r}. "
                "Must be one of 'tau/bar', 'tau/km', 'tau', or 'cm2/g'."
            )

        if style not in {"line", "step"}:
            raise ValueError(
                f"Invalid style {style!r}. Must be either 'line' or 'step'."
            )

        # Retrieve prior profiles only when requested.
        if plot_prior_profiles:
            priors, _ = self.core.generate_prior_profiles()
            priors = priors.data

        # Select the y-axis once.
        y_column = "pressure" if pressure else "height"

        # Labels and conversion attributes for each supported unit.
        unit_config = {
            "tau/bar": ("Optical depth / bar", 1),
            "tau/km": ("Optical depth / km", 0),
            "tau": ("Pressure-integrated optical depth", 2),
            "cm2/g": ("Aerosol density (cm$^2$ / gram)", None),
        }

        unit_label, conversion_index = unit_config[unit]
        unit_label += f" at {self.core.reference_wavelength:.1f}µm"

        # Plot each retrieved aerosol profile.
        for label in self.retrieved_aerosols.columns:
            if label in {"height", "pressure"}:
                continue

            retrieved = self.retrieved_aerosols[label]
            y_retrieved = self.retrieved_aerosols[y_column]

            if conversion_index is not None:
                retrieved = self._convert_aerosol_units(retrieved)[conversion_index]

            if style == "line":
                (line,) = ax.plot(retrieved, y_retrieved, label=label)
            else:
                (line,) = ax.step(retrieved, y_retrieved, label=label)

            # Plot prior using the same colour as the retrieved profile.
            if plot_prior_profiles:
                prior = priors[label]
                y_prior = priors[y_column]

                if conversion_index is not None:
                    prior = self._convert_aerosol_units(prior)[conversion_index]

                ax.plot(
                    prior,
                    y_prior,
                    color=line.get_color(),
                    ls="dotted",
                    label=f"{label} Prior" if not hide_extra_labels else "",
                )

        ax.set_xlabel(unit_label)
        ax.set_xscale("log")
        ax.set_xlim(1e-3, 1e2)
        ax.legend(loc="lower left")

    @utils.plotting_altitude
    def plot_gas_profiles(self, ax, 
                          pressure, 
                          unit="", 
                          gas_names=None,
                          plot_retrieved_profiles=True,
                          plot_prior_profiles=True, 
                          plot_ref_profiles=False,
                          hide_extra_labels=True):
        
        """Plot gas profiles from the .prf file.

        Args:
            ax (matplotlib.Axes): The matplotlib.Axes object to plot to. If omitted then create a new Figure and Axes.
            pressure (bool): Whether to plot the gas profiles against pressure (if True) or height (if False).
            unit (str): One of '', 'ppm', 'ppb' or 'ppt'. 
            gas_names (list[str], optional): List of gas names to plot. If None, plot all gases.
            plot_retrieved_profiles (bool): Whether to plot the retrieved gas profiles
            plot_prior_profiles (bool): Whether to plot the prior gas profiles
            plot_ref_profiles (bool): Whether to plot the profile in the .ref file
            hide_extra_labels (bool): Whether to show labels for eveything in the legend (False) or only the species colour (True)

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
        y2 = self.core.reference.data["pressure"] if pressure else self.core.reference.data["height"]

        # Get the prior distributions if requested
        if plot_prior_profiles:
            _, priors = self.core.generate_prior_profiles()
            priors = priors.data
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
                ax.plot(self.core.reference.data[gas_name]*scale, 
                        y2, label=f"{gas_name} .ref" if not hide_extra_labels else "",  
                        c=c, 
                        ls="dashed")
            if plot_prior_profiles:
                ax.plot(priors[gas_name]*scale, 
                        y3, 
                        label=f"{gas_name} Prior" if not hide_extra_labels else "", 
                        c=c, 
                        ls="dotted")

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

    def make_summary_plot(self, figsize=(11, 10), savefig=True, log=False, gas_kwargs=dict(), aerosol_kwargs=dict()):
        """Make a summary plot with prior and retrieved spectra, error on the spectra, aerosol and chemical profiles, and chi-squared values.
        
        Args:
            figsize (int, int): matplotlib figure size
            savefig (bool): Whether to save the figure
            log (bool): Whether to use a log plot for the spectrum and residuals
            gas_kwargs (dict): Keyword arguments for gas profile plotting
            aerosol_kwargs (dict): Keyword arguments for aerosol profile plotting
            
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
        axs["A"].sharex(axs["B"])

        self.plot_spectrum(ax=axs["A"], log=log)
        self.plot_spectrum_residuals(ax=axs["B"], log=log)
        self.plot_chisq(ax=axs["C"])
        self.plot_aerosol_profiles(ax=axs["D"], **aerosol_kwargs)
        self.plot_gas_profiles(ax=axs["E"], gas_names=names, **gas_kwargs)

        utils.label_axes(axs)

        # if self.core.forward:
        #     fig.suptitle(f"Forward model in {self.core_directory.resolve()}", y=0.91)
        # else:
        #     fig.suptitle(f"Retrieval in {self.core_directory.resolve()}", y=0.91)

        if savefig: 
            fig.savefig(self.core_directory / f"plots/summary{'_log' if log else ''}.png", bbox_inches="tight", dpi=400)


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


class NemesisResultGroup(list):
    """Subclass of `list` with a few extra routines to make extracting parameters easier"""
    def __init__(self, iterable=list()):
        super().__init__(iterable)

    def get_parameter(self, profile_name=None, parameter_name=None, return_results=True):
        """Get all the values for a given paramter in a profile across the cores.
        
        Args:
            profile_name (str or mcmc.Parameter): The label of the profile (or a mcmc.Paramter object)
            parameter_name (str): The name of the parameter within the given profile (if profile_name is a 
                                  mcmc.Paramter object this is not required)
                                  
        Returns:
            dict: A dictionary containing all the cores as keys and the associated parameter values"""
        
        if isinstance(profile_name, mcmc.Parameter):
            parameter_name = profile_name.param_name
            profile_name = profile_name.profile_name

        if profile_name.lower() == "chi_sq":
            out = dict()
            for res in self:
                out[res] = res.chi_sq
        elif profile_name is not None and parameter_name is not None:
            out = dict()
            for res in self:
                profile = res.profiles[profile_name]
                param = getattr(profile, parameter_name)
                out[res] = param
        else:
            raise ValueError("Must specify a profile name and parameter name, or just paramter_name='chi_sq'")

        if return_results:
            return out
        else:
            return list(out.values())

    def get_by_parameter(self, profile_name, parameter_name, value):
        """Get a list of all cores with a given parameter value

        Args:
            profile_name (str or mcmc.Parameter): The label of the profile (or a mcmc.Paramter object)
            parameter_name (str): The name of the parameter within the given profile (if profile_name is a 
                                  mcmc.Paramter object this is not required)
            value (int/float or tuple): The parameter value to search for. If it is a float, search for that exact value. 
                                    If it is a tuple, return all results between value[0] and value[1]

        Returns:
            NemesisResultGroup: A filtered list of the cores
        """

        if isinstance(profile_name, mcmc.Parameter):
            parameter_name = profile_name.param_name
            profile_name = profile_name.profile_name

        out = NemesisResultGroup()
        for res in self:
            profile = res.profiles[profile_name]
            param = getattr(profile, parameter_name)
            if isinstance(value, (float, int)):
                if param == value:
                    out.append(res)
            else:
                if value[0] < param < value[1]:
                    out.append(res)

                
        return out

    def get_central_wavelength(self):
        """Get the central wavelength in the spx file for each core

        Returns:
            dict: Dictionary containing cores as keys and central wavelengths as values
        """
        out = dict()
        for r in self:
            wvs = r.retrieved_spectrum.wavelength
            out[r.core] = np.median(wvs)
        return out

    def get_by_id(self, id_):
        """Gett the NemesisResult object for a given core ID"""
        for r in self:
            if r.core.id_ == id_:
                return r

    def get_best_fit(self, n=0):
        """Return the best fit; lowest chi squared."""
        chis = self.get_parameter("chi_sq")
        ranked = sorted(chis.items(), key=lambda x: x[1])
        if not 0 <= n < len(ranked):
            raise IndexError(f"Rank {n} out of range for {len(ranked)} fits")
        core,_ = ranked[n]
        return core

    def get_worst_fit(self, n=0):
        """Find the Fit with the highest chi squared value"""
        chis = self.get_parameter("chi_sq")
        ranked = sorted(chis.items(), key=lambda x: x[1], reverse=True)
        if not 0 <= n < len(ranked):
            raise IndexError(f"Rank {n} out of range for {len(ranked)} fits")
        core,_ = ranked[n]
        return core

    @utils.plotting
    def plot_grid_search(self, 
                         x_profile_name,
                         x_parameter_name,
                         x_values, 
                         y_profile_name,
                         y_parameter_name,
                         y_values,
                         ax,
                         cbar=True,
                         **kwargs):
        
        out = np.full((len(x_values), len(y_values)), np.nan)
        for i, x in enumerate(x_values):
            r = self.get_by_parameter(x_profile_name, x_parameter_name, x)
            for j, y in enumerate(y_values):
                z = r.get_by_parameter(y_profile_name, y_parameter_name, y)
                if len(z) == 1:
                    out[i,j] = z[0].chi_sq
                elif len(z) != 0:
                    raise ValueError(f"Multiple cores found with {x=} and {y=}")

        pcm = ax.pcolormesh(x_values, y_values, out.T, **kwargs)
        ax.set_xlabel(f"{x_profile_name} {x_parameter_name}")
        ax.set_ylabel(f"{y_profile_name} {y_parameter_name}")
        if cbar:
            plt.colorbar(pcm, ax=ax)
        return pcm


class SensitivityAnalysis:
    """Class for analysing the sensitivity cores generated by `eleos.cores.create_sensitivity_analysis`
    
    Attributes:
        parent_directory (str): The directory containing the sensitivity cores
        baseline (NemesisResult): The baseline result, an alias for results[0]
        params (pandas.DataFrame): A DataFrame detailing which parameters were varied in each core and their value
    """

    def __init__(self, parent_directory):
        """
        Args:
            parent_directory (str): The directory containing all the individual cores
            
        """
        self.parent_directory = Path(parent_directory)
        self.baseline = NemesisResult(self.parent_directory / "core_1")
        self.params = self._get_results()

    def _get_results(self):
        params = pd.read_csv(self.parent_directory / "sensitivity_analysis.txt")
        results = load_multiple_cores(self.parent_directory, failed='warn')
        result_by_core_id = {
            result.core.id_: result
            for result in results
        }
        params["Result"] = params["Core ID"].map(result_by_core_id)
        baseline_row = pd.DataFrame([{
            "Core ID": 1,
            "Profile Label": "Baseline",
            "Result": self.baseline,
        }])

        df = pd.concat([baseline_row, params], ignore_index=True)
        return df

    def _get_all_param_names(self):
        """Get a list of all the parameter names in the sensitivity analysis
        
        Returns:
            list[tuple]: A list of tuples containing the profile label and parameter name for each parameter in the sensitivity analysis"""

        lst = [
            (row["Profile Label"], row["Parameter"])
            for _, row in self.params.iterrows()
        ]
        return list(sorted(list(set(lst))))

    @utils.plotting
    def make_chisq_plot(self, param, ax, retrieved=True):
        if isinstance(param, (tuple, list)):
            param = mcmc.Parameter(*param)

        df = self.params[
            ((self.params["Profile Label"] == param.profile_name)
            & (self.params["Parameter"] == param.param_name))
            | (self.params["Core ID"] == 1)
        ]   

        values = []
        chis = []
        print(df["Result"])
        for res in df["Result"]:
            if retrieved:
                values.append(getattr(res.profiles[param.profile_name], "retrieved_" + param.param_name))
            else:
                values.append(getattr(res.profiles[param.profile_name], param.param_name))
            chis.append(getattr(res, "chi_sq"))
        ax.axvline(getattr(res.profiles[param.profile_name], param.param_name))
        
        ax.scatter(values, chis)
        if df["Log"].any():
            ax.set_xscale("log")
        ax.set_xlabel(f"{param.profile_name} {param.param_name}")
        ax.set_ylabel("Chi-squared value")

    @utils.plotting
    def plot_parameter(self, profile_label, parameter, ax, mode="difference"):
        """Plot the sensitivity of the model to the given parameter in the given profile.
        
        Args:
            ax: The matplotlib.Axes object to plot to. If omitted then create a new Figure and Axes
            profile_label: The label of the profile to plot
            parameter: The variable to plot
            mode (str): One of "absolute", "difference", "percentage": Absolute 
            
        Returns:
            matplotlib.Figure: The Figure object to which the Axes belong
            matplotlib.Axes: The Axes object onto which the data was plotted"""

        # get relevant rows
        df = self.params[
            (self.params["Profile Label"] == profile_label)
            & (self.params["Parameter"] == parameter)
        ]

        # create cmap for lines (red = increase in parameter value, blue is decrease)
        if df["Log"].all():
            cmapname = "PiYG_r"
        else:
            cmapname = "BrBG_r"
        cmap = mpl.colormaps[cmapname].resampled(256)
        blue = cmap(np.linspace(0.1, 0.4, 128))
        red = cmap(np.linspace(0.6, 0.9, 128))
        blue = blue[::-1]
        red = red[::-1]
        cmap = mpl.colors.LinearSegmentedColormap.from_list(
            "sensitivity",
            np.vstack([blue, red]),
        )
        norm = mpl.colors.TwoSlopeNorm(
            vmin=df["Factor"].min(),
            vcenter=1.0,
            vmax=df["Factor"].max(),
        )

        for i, row in df.iterrows():
            if row["Profile Label"] == profile_label and row["Parameter"] == parameter:
                result = row["Result"]
                if isinstance(result, float) and np.isnan(result):
                    continue
                if mode == "absolute":
                    y = result.retrieved_spectrum.model
                    b = self.baseline.retrieved_spectrum.model
                elif mode == "difference":
                    y = result.retrieved_spectrum.model - self.baseline.retrieved_spectrum.model
                    b = np.zeros(len(y))
                elif mode == "percentage":
                    y = (result.retrieved_spectrum.model - self.baseline.retrieved_spectrum.model) / self.baseline.retrieved_spectrum.model
                    b = np.zeros(len(y))
                    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=1))
                color = cmap(norm(row["Factor"]))
                ax.plot(result.retrieved_spectrum.wavelength, y, label=row['Factor'], color=color, ls="solid")
                ax.plot(result.retrieved_spectrum.wavelength, b, color="k", ls="solid")

    def make_parameters_plot(self, ncol=4, mode="difference"):
        p = [pp for pp in sorted(self._get_all_param_names()) if pp[0] != "Baseline"]

        nrow = int(np.ceil(len(p) / ncol))
        fig, axs = plt.subplots(nrow , ncol, figsize=(4*ncol, 1.5*nrow), sharex=True)
        axs = axs.flatten()

        for ax, name in zip(axs, p):
            self.plot_parameter(*name, ax=ax, mode=mode)
            ax.text(0.99, 0.05, " ".join(name).replace("_", " "), transform=ax.transAxes, ha="right", va="bottom")

        # show x ticks on plots at bottom of each column
        for ax in axs[len(p)-4:len(p)]:
            ax.tick_params(axis="x", which="both", labelbottom=True)

        # hide any unused axes
        for ax in axs[len(p):]:
            ax.set_visible(False)

        fig.supxlabel("Wavelength (µm)")
        if mode == "difference":
            fig.supylabel("Radiance change from baseline", x=0.01)
        elif mode == "absolute":
            fig.supylabel("Radiance", x=0.01)
        elif mode == "percentage":
            fig.supylabel("Percentage radiance change from baseline", x=0.01)

        fig.tight_layout()

    def savefig(self, name, fig=None, **kwargs):
        """
        Save a matplotlib figure to a file in the plots folder in the parent directory.

        Args:
            name (str): The name of the file to save the figure as.
            fig (matplotlib.Figure, optional): The figure to save. If None, the current figure will be saved. Default is None.
            **kwargs: Additional keyword arguments to pass to `savefig`.

        Returns:
            None
        """
        (self.parent_directory / "plots").mkdir(exist_ok=True)
        x = plt if fig is None else fig
        x.savefig(self.parent_directory / "plots" / name, bbox_inches="tight", **kwargs)


class GasAnalysis:
    def __init__(self, parent_directory):
        self.parent_directory = Path(parent_directory)
        self.gas_index = pd.read_csv(self.parent_directory / "gasindex.csv", keep_default_na=False, na_values=['NaN'])
        self.spectra = dict()
        itt = list(self.gas_index.iterrows())
        for _, (core_id, excluded_gas) in tqdm(itt, desc="Loading cores"):
            r = NemesisResult(self.parent_directory / f"core_{core_id}")
            self.spectra[excluded_gas] = r.retrieved_spectrum
            if excluded_gas == "None":
                self.baseline = r

    def _get_biggest_contributors(self):
        baseline = self.spectra["None"]
        biggas = []
        for i in range(len(baseline.wavelength)):
            max_diff = 0
            max_gas = None
            for gas, spectrum in self.spectra.items():
                if gas == "None":
                    continue
                diff = np.abs(spectrum.model.iloc[i] - baseline.model.iloc[i])
                if diff > max_diff:
                    max_diff = diff
                    max_gas = gas
            biggas.append(max_gas)
        return biggas

    @utils.plotting
    def plot_contributions_2D(self, ax):
        ax.set_ylabel("Gas contribution (arb. units)")
        reference_spectrum = self.spectra["None"]
        for gas, spectrum in self.spectra.items():
            if gas == "None":
                continue
            diff = spectrum.model - reference_spectrum.model
            ax.plot(spectrum.wavelength, diff, label=gas)
        ax.legend(ncols=2)

    @utils.plotting
    def plot_contributions_3D(self, ax, normalise=False, log_scale=False, cmap=cmasher.sunburst_r, wl_range=None, hide=list()):
        diffs = []
        reference_spectrum = self.spectra["None"]

        if wl_range is not None:
            mask = ((reference_spectrum.wavelength > wl_range[0]) & (reference_spectrum.wavelength < wl_range[1]))
        else:
            mask = np.ones(len(reference_spectrum.wavelength), dtype=bool)

        gases = []
        scores = []

        for gas, spectrum in self.spectra.items():
            if gas == "None" or gas in hide:
                continue
            gases.append(gas)
            diff = spectrum.model[mask] - reference_spectrum.model[mask]
            diff = np.clip(diff, 1e-20, np.inf)
            diffs.append(diff)

        gases = np.asarray(gases)
        diffs = np.stack(diffs)
        scores = np.nansum(diffs**2, axis=1)
        valid = scores > 0
        order = np.argsort(scores[valid])
        gases = gases[valid][order]
        diffs = diffs[valid][order]

        if normalise:
            factors = np.nanmax(np.abs(diffs), axis=1)
            factors /= np.nanmax(factors)
            factors = 1/factors

            diffs = np.multiply(
                diffs,
                factors[:, None],
                out=np.zeros_like(diffs),
            )

        if log_scale:
            norm = mpl.colors.LogNorm(
                vmin=max(np.nanmin(diffs), 1e-10),
                vmax=np.nanmax(diffs),
            )
        else:
            norm = None

        gas_labels = [constants.GAS_LABELS[l] for l in gases]

        ax.pcolormesh(
            reference_spectrum.wavelength[mask],
            gas_labels,
            diffs,
            cmap=cmap,
            norm=norm,
        )

        if normalise:
            def sci_label(f):
                if f == 0:
                    return r"$0$"
                if 0.01 < f < 100:
                    return f"{f:.2f}"
                exponent = int(np.floor(np.log10(abs(f))))
                mantissa = f / 10**exponent
                return rf"${mantissa:.1f}\times10^{{{exponent}}}$"

            ax_right = ax.twinx()
            ax_right.set_ylim(ax.get_ylim())
            ax_right.set_yticks(np.arange(len(gas_labels)))
            ax_right.set_yticklabels([sci_label(f) for f in factors])
            ax_right.set_ylabel("Normalization factor")

    @utils.plotting
    def plot_highlighted_spectrum(self, ax):
        baseline = self.spectra["None"]
        biggas = self._get_biggest_contributors()
        ax.plot(baseline.wavelength, baseline.model, c='k', lw=1, marker="o")
        ax.set_ylabel("Radiance\n(μW cm$^{-2}$ sr$^{-1}$ μm$^{-1}$)")
        ylim = ax.get_ylim()

        # Get unique gases and create a mapping to numeric indices
        gas_to_index = {gas: i for i, gas in enumerate(set(biggas))}
        
        # Create gradient data - map each gas to its index
        gradient = np.array([gas_to_index[gas] for gas in biggas]).reshape(1, -1)
        gradient = np.vstack((gradient, gradient))
        
        # Create a colormap with unique colors for each gas
        cmap = plt.cm.get_cmap('tab20', len(gas_to_index))
        
        # Add the gradient as background
        ax.pcolormesh(baseline.wavelength, [0, 1], gradient, cmap=cmap, zorder=-1)
    
        # Create a legend
        legend_elements = [mpl.patches.Patch(facecolor=cmap(i), label=gas) for gas, i in gas_to_index.items()]
        ax.legend(handles=legend_elements, loc='upper right')


class ScattershotMC(NemesisResultGroup):

    @classmethod 
    def load(cls, cd, **kwargs):
        ress = load_multiple_cores(cd, **kwargs)
        return cls(ress)

    def plot_retrieval_map(self, x_parameter, y_parameter, log_chi_sq=True, max_chi_sq=5, **kwargs):
        x_prior = self.get_parameter(x_parameter)
        y_prior = self.get_parameter(y_parameter)
        x_ret   = self.get_parameter(x_parameter.ret())
        y_ret   = self.get_parameter(y_parameter.ret())
        chi = self.get_parameter("chi_sq")


        keep = {result for result, value in chi.items() if value < max_chi_sq}
        x_prior = [v for r, v in x_prior.items() if r in keep]
        y_prior = [v for r, v in y_prior.items() if r in keep]
        x_ret   = [v for r, v in x_ret.items() if r in keep]
        y_ret   = [v for r, v in y_ret.items() if r in keep]
        chi     = [v for r, v in chi.items() if r in keep]

        fig, ax = plt.subplots(1, 3, figsize=(8,3), width_ratios=(15,15,1))

        ax[0].scatter(x_prior, y_prior, **kwargs)
        norm = mpl.colors.Normalize() if not log_chi_sq else mpl.colors.LogNorm()
        s = ax[1].scatter(x_ret, y_ret, c=chi, norm=norm, cmap="magma",  **kwargs)

        ax[1].yaxis.set_label_position("right")
        ax[1].yaxis.tick_right()

        ax[0].set_xlabel(f"Prior {x_parameter}")
        ax[0].set_ylabel(f"Prior {y_parameter}")
        ax[1].set_xlabel(f"Retrieved {x_parameter}")
        ax[1].set_ylabel(f"Retrieved {y_parameter}")

        fig.colorbar(s, cax=ax[2], label=r"$\chi^2$ value")
        fig.tight_layout()

        pos0 = ax[0].get_position()
        pos1 = ax[1].get_position()
        gap = 0.05  # desired gap in figure coordinates
        ax[1].set_position([
            pos0.x1 + gap,
            pos1.y0,
            pos1.width,
            pos1.height
        ])

        return fig, ax[:-1]
    

def load_multiple_cores(parent_directory, failed='warn', show_progress=True):
    """Read in all the cores in a given directory and return a list of NemesisResult objects.
    
    Args:
        parent_directory (str): The directory containing all the individual core directories
        failed (str): Whether to raise an error if a retieval failed ('raise') or to warn and skip it ('warn'), or silently skip ('skip') 
        
    Returns:
        NemesisResultGroup: A subclass of list containing the result object for each core"""
    
    parent_directory = Path(parent_directory)

    assert failed in ("raise", "warn", "skip"), "failed parameter must be one of 'raise', 'warn' or 'skip'"

    out = []
    iterator = sorted(parent_directory.glob("core_*"), key=sort_key_paths)
    for core in tqdm(iterator, desc="Loading cores"):
        try:
            out.append(NemesisResult(core))
        except Exception as e:
            if failed == "warn":
                warnings.warn(f"Failed to load core {core}")
            
            if failed == "raise":
                raise e
            
    return NemesisResultGroup(out)


def load_best_cores(parent_directory, n, **kwargs):
    """Load n cores from the parent_directory with the lowest chi-squared values.

    Args:
        parent_directory (str): The directory containing all the individual core directories
        n (int): The number of cores to load
        **kwargs: Other keyword args passed to load_multiple_cores()
    
    Returns:
        NemesisResultGroup: A subclass of list containing the result object for each core, sorted by chi-squared value (so lowest chi-sqared is index 0 in the list)
    """
    
    ress = load_multiple_cores(parent_directory, **kwargs)
    ress = sorted(ress, key=lambda x: x.chi_sq)

    return ress[:n]


def load_parallelised_cores(parent_directory):
    """Given a directory of cores generated with cores.parallelise_forward, load them and recombine into a
    single NemesisResult object. The chi-sq of the core will be set to None, and the chi_sqs attribute will 
    be set to the chi squared values of each chunk (ie, chi squared as a coarse function of wavelength)
    
    Args:
        parent_directory (str): Path to the directory of cores
        
    Returns:
        NemesisResult: The recombined results object
    """

    ress = load_multiple_cores(parent_directory, failed="raise")
    retspec = None
    chis = []
    for r in ress:
        if retspec is None:
            retspec = r.retrieved_spectrum
        else:
            retspec = pd.concat([retspec, r.retrieved_spectrum], ignore_index=True)
        chis.append(r.chi_sq)
    ress[0].retrieved_spectrum = retspec
    ress[0].chi_sq = None
    ress[0].chi_sqs = chis
    return ress[0]


def sort_key_paths(path):
    x = str(path).split("_")[-1]
    return int(x)