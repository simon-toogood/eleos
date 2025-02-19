"""This module contains the classes for creating Profile objects."""

from collections import defaultdict
from pathlib import Path
import re

from . import constants
from . import shapes
from . import utils
from . import results
from . import spx


class Profile:
    def __init__(self, label=None):
        """This is the base class for all Profile objects. It should never be instantiated directly, use a subclass such as ``GasProfile`` or ``TemperatureProfile``"""
        self.retrieved = False
        self.core = None
        self.label = label

    def __str__(self):
        if not self.retrieved:
            headers = ["Prior", "Error"]
        else:
            headers = ["Prior", "Error", "Retrieved", "Error"]

        names = self._get_displayable_attributes()
        data = []
        for name_group in names:
            data.append([])
            for i in range(len(headers)):
                try:
                    data[-1].append(getattr(self, name_group[i]))
                except:
                    data[-1].append("NA")
        
        return utils.generate_ascii_table(self.label, headers, [n[0] for n in names], data)

    @classmethod
    def from_previous_retrieval(cls, core_directory, label=None):
        """Create a Profile object using the retrieved parameters from a previous retrieval as priors. Use either id or label to specify the profile to use in the previous retrieval.
        
        Args:
            core_directory: The core directory of the previous retrieval
            label: The label of the profile to use in the previous retrieval"""
        
        res = results.NemesisResult(core_directory)

        prev_profile = None
        for profile in res.profiles:
            if profile.label == label and cls == type(profile):
                prev_profile = profile
                break
            
        if prev_profile is None:
            raise ValueError(f"Profile with label {label} not found in previous retrieval")
        
        new_profile = cls._create_profile_from_previous_retrieval(prev_profile)
        new_profile.shape._set_prior_to_retrieved()

        return new_profile

   
        # Toggle retrieved flag
        self.retrieved = True

    def _clear_result(self):
        """Clear the retrieved values from the profile
        
        Args:
            None
            
        Returns:
            None"""
        
        for name in self.shape.NAMES:
            for attr in [f"retrieved_{name}", f"retrieved_{name}_error"]:
                delattr(self.shape, attr)
        self.retrieved = False


class TemperatureProfile(Profile):
    def __init__(self, filepath, label=None):
        """Create a temperature profile from a prior file
        
        Args:
            filepath: The filepath of the prior temperature profile
            label (str): A label to associate with this Profile. By default it is "Temperature"
        """
        super().__init__(label)
        if self.label is None:
            self.label = "Temperature"
        self.shape = shapes.Shape0(filepath=filepath)

    def __repr__(self):
        return f"<TemperatureProfile [{self.create_nemesis_string()}]>"

    def _get_displayable_attributes(self):
        return []

    def _add_result(self, df):
        self.shape.data = df
        self.retrieved = True

    def create_nemesis_string(self):
        """Create the NEMESIS code that represents the temperature profile. Temperature profiles currently only support mode 0 0 0 so this function's return value is always constant
        
        Args:
            None
            
        Returns:
            str: '0 0 0'"""
        return f"0 0 0"
    
    def generate_apr_data(self):
        """Generate the section of the .apr file for this temperature profile
        
        Args:
            None
            
        Returns:
            str: The string to write to the .apr file"""
        return "0 0 0 - Temp\n" + self.shape.generate_apr_data()

    def get_name(self):
        if self.label is None:
            return "Temperature"
        else:
            return self.label


class GasProfile(Profile):
    def __init__(self, 
                 gas_name=None, 
                 gas_id=None, 
                 isotope_id=0, 
                 shape=None, 
                 label=None):
        """Create a profile for a given gas (optionally an isotopologue) with a given shape
        
        Args:
            gas_name (str): The name of the gas (eg 'CH4'). Specify either this OR gas_id
            gas_id (int): The radtrans ID of the gas (eg. 6). Specify either this OR gas_name
            isotope_id (int): The ID of the isotopologue to use. Use 0 for a mix of all at terrestrial abundance
            shape (Shape): A Shape object to use for the profile shape
            label (str): A label to associate with this Profile. By default it is "<gas_name> <isotope_id>" (eg. "PH3 0")
        """
        
        super().__init__(label=label)

        self.isotope_id = isotope_id

        if shape is None:
            raise ValueError("shape attribute must be specified")
        self.shape = shape
        self.shape.share_parameters(self)

        if not ((gas_id is None) ^ (gas_name is None)):
            raise ValueError("Specifiy exactly one of gas_name or gas_id (not both!)")
        if gas_name is None:
            self.gas_id = gas_id
            self.gas_name = constants.GASES.loc[constants.GASES.radtrans_id == gas_id].name.iloc[0]
        elif gas_id is None:
            self.gas_name = gas_name
            self.gas_id = constants.GASES.loc[constants.GASES.name == gas_name].radtrans_id.iloc[0]
        
        if label is None:
            self.label =  f"{self.gas_name} {self.isotope_id}"

    def __repr__(self):
        return f"<GasProfile {self.gas_name} [{self.create_nemesis_string()}]>"

    def _get_displayable_attributes(self):
        def sort_key(name):
            prefix = name.startswith("retrieved_")
            error = name.endswith("_error")
            return (prefix * 2 + error, name)
        
        groups = defaultdict(list)

        # Populate dictionary
        for name in self.__dict__:
            core_name = re.sub(r"^retrieved_", "", name).split('_')[0]
            groups[core_name].append(name)

        # Convert to list of grouped names
        grouped = [sorted(x, key=sort_key) for x in list(groups.values())]

        return [g for g in grouped if g[0] not in ("label", "gas_id", "isotope_id", "retrieved", "core")]

    def _create_profile_from_previous_retrieval(prev_profile):
        return GasProfile(gas_name=prev_profile.gas_name, isotope_id=prev_profile.isotope_id, shape=prev_profile.shape, label=prev_profile.label)

    def _add_result(self, df):
        """Take in a DataFrame created by reading in the .mre file (results.NemesisResult.read_mre) and assign the correct attributes
        
        Args:
            df: pandas.DataFrame with columns "prior", "prior_error", "retrieved", "retrieved_error" and a row for each parameter
            
        Returns:
            None"""
        
        # Get order of parameters
        names = self.shape.NAMES
        assert len(names) == len(df)

        # Set attributes of the child Shape object
        for name, (_, row) in zip(names, df.iterrows()):
            for title, value in zip(row.index, row.values):
                if "prior" in title:
                    attrname = title.replace("prior", name)
                elif "retrieved" in title:
                    attrname = title.replace("retrieved", f"retrieved_{name}")
                setattr(self.shape, attrname, value)

        self.retrieved = True
        self.shape.share_parameters(self)

    def create_nemesis_string(self):
        """Create the NEMESIS code that represents the gas profile (eg. 23 0 1)
        
        Args:
            None
            
        Returns:
            str: The NEMESIS code"""
        return f"{self.gas_id} {self.isotope_id} {self.shape.ID}"
    
    def generate_apr_data(self):
        """Generate the section of the .apr file for this gas profile
        
        Args:
            None
            
        Returns:
            str: The string to write to the .apr file"""
        return self.create_nemesis_string() + " - " + self.gas_name + "\n" + self.shape.generate_apr_data()


class AerosolProfile(Profile):
    def __init__(self, 
                 shape, 
                 radius, 
                 variance,
                 real_n,
                 imag_n, 
                 retrieve_optical=False, 
                 radius_error=None,
                 variance_error=None,
                 imag_n_error=None,
                 label=None):
        """Create a profile for a given aerosol with a given shape
        
        Args:
            shape (Shape): A Shape object to use for the profile shape
            label (str): A label to associate with this Profile. By default it is "Aerosol <aerosol_id>" (eg. "Aerosol 1") """            
        
        super().__init__(label=label)
        if label is None:
            self.label = f"Aerosol {self.aerosol_id}"

        self.aerosol_id = "UNASSIGNED"
        self.shape = shape
        self.shape.share_parameters(self)

        self.radius = radius
        self.variance = variance
        self.real_n = real_n
        self.imag_n = imag_n
        self.retrieve_optical = retrieve_optical

        if retrieve_optical:
            self.radius_error = radius_error
            self.variance_error = variance_error
            self.imag_n_error = imag_n_error
            self.NAMES = self.NAMES + ["radius", "variance", "imag_n"]

    def __repr__(self):
        return f"<AerosolProfile {self.aerosol_id} [{self.create_nemesis_string()}]>"

    def _get_displayable_attributes(self):
        def sort_key(name):
                prefix = name.startswith("retrieved_")
                error = name.endswith("_error")
                return (prefix * 2 + error, name)
        
        groups = defaultdict(list)

        # Populate dictionary
        for name in self.__dict__:
            core_name = re.sub(r"^retrieved_", "", name).split('_')[0]
            groups[core_name].append(name)

        # Convert to list of grouped names
        grouped = [sorted(x, key=sort_key) for x in list(groups.values())]

        return [g for g in grouped if g[0] not in ("label", "retrieved", "core", "aerosol_id")]

    def _add_result(self, df, df_444=None):
        """Take in a DataFrame created by reading in the .mre file (results.NemesisResult.read_mre) and assign the correct attributes
        
        Args:
            df: pandas.DataFrame with columns "prior", "prior_error", "retrieved", "retrieved_error" and a row for each parameter
            df_444: Same as df, but this is for the 444 profile if retrieving optical properties
        
        Returns:
            None"""

        # Get order of parameters
        names = self.shape.NAMES
        assert len(names) == len(df)

        # Set attributes of the child Shape object
        for name, (_, row) in zip(names, df.iterrows()):
            for title, value in zip(row.index, row.values):
                if "prior" in title:
                    attrname = title.replace("prior", name)
                elif "retrieved" in title:
                    attrname = title.replace("retrieved", f"retrieved_{name}")
                setattr(self.shape, attrname, value)
    
        if df_444 is not None:
            for name, (_, row) in zip(["radius", "variance", "imag_n"], df_444.iterrows()):
                for title, value in zip(row.index, row.values):
                    if "prior" in title:
                        attrname = title.replace("prior", name)
                    elif "retrieved" in title:
                        attrname = title.replace("retrieved", f"retrieved_{name}")
                    setattr(self, attrname, value)

        self.retrieved = True
        self.shape.share_parameters(self)

    def create_nemesis_string(self):
        """Create the NEMESIS code that represents the aerosol profile (eg. -1 0 32)
        
        Args:
            None
            
        Returns:
            str: The NEMESIS code"""
        return f"-{self.aerosol_id} 0 {self.shape.ID}"
    
    def generate_apr_data(self):
        """Generate the section of the .apr file for this aerosol profile
        
        Args:
            None
            
        Returns:
            str: The string to write to the .apr file"""
        
        aerosol_part = self.create_nemesis_string() + f" - {self.label}\n" + self.shape.generate_apr_data()

        if self.retrieve_optical:
            imagn_part = f"\n444 {self.aerosol_id} 444 - {self.label}\ncloudf{self.aerosol_id}.dat"
            return aerosol_part + imagn_part
        else:
            return aerosol_part

    def _generate_cloudfn_dat(self, directory):
        # Get first wavelength
        directory = Path(directory)
        wl = spx.read(directory / "nemesis.spx").geometries[0].wavelengths[0]

        # Get number of wavelengths in xsc file
        with open(directory / "nemesis.xsc") as file:
            lines = file.read().split("\n")
            nwave = (len(lines)-2) // 2
            refwave = float(lines[1].split()[0])
        
        # Generate cloudfN.dat
        with open(directory / f"cloudf{self.aerosol_id}.dat", mode="w+") as file:
            utils.write_nums(file, self.radius, self.radius_error)
            utils.write_nums(file, self.variance, self.variance_error)
            file.write(f"{nwave}    -1\n")
            utils.write_nums(file, refwave, self.real_n)
            utils.write_nums(file, refwave)
            for line in lines:
                vals = line.split()
                if len(vals) < 2:
                    continue
                utils.write_nums(file, float(vals[0]), self.imag_n, self.imag_n_error)
