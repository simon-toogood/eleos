"""This module contains the classes for creating Profile objects."""

from itertools import zip_longest
from collections import defaultdict
from pathlib import Path

from . import constants
from . import shapes
from . import utils
from . import results
from . import spx

# idea for the future: combine AerosolProfile and ImagRefractiveIndexProfile into a single class with a toggle to retrieve
# the radius, variance, imag refractive index

# not happy with the implementation for setting priors from previous retrieval...

class Profile:
    def __init__(self, label=None):
        """This is the base class for all Profile objects. It should never be instantiated directly, use a subclass such as ``GasProfile`` or ``TemperatureProfile``"""
        self.retrieved = False
        self.core = None
        self.label = label

    def __str__(self):
        # generate table title
        try:
            self.label
        except:
            self.label = None

        attrs = []
        for name in self.shape.NAMES:
            if self.retrieved:
                attrs.append([f"{name}", f"{name}_error", f"retrieved_{name}", f"retrieved_{name}_error"])
                headers = ["Prior", "Error", "Retrieved", "Error"]
            else:
                attrs.append([f"{name}", f"{name}_error"])
                headers = ["Prior", "Error"]

        # extract all the attribute values
        values = []
        for attr_list in attrs:
            values.append([])
            for attr in attr_list:
                values[-1].append(getattr(self.shape, attr))
    
        return utils.generate_ascii_table(self.get_name(), headers, self.shape.NAMES, values)

    def compact_str(self):
        out = f"ID: {self.shape.ID}\n"
        for k, v in self.__dict__.items():
            if k == "shape":
                for name in self.shape.NAMES:
                    out += f"{name}: {getattr(self.shape, name)}±{getattr(self.shape, f'{name}_error')}\n"
                    if self.retrieved:
                        f"{getattr(self.shape, f'retrieved_{name}')}±{getattr(self.shape, f'retrieved_{name}_error')}\n"
            elif k in ("core", "retrieved", "shape"):
                pass
            else:
                out += f"{k}: {v}\n"
        return out.rstrip("\n")

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
    def __init__(self, filepath, **kwargs):
        """Create a temperature profile from a prior file
        
        Args:
            filepath: The filepath of the prior temperature profile
            label: (optional) An arbitrary label to associate with this profile
        """
        super().__init__(**kwargs)
        self.shape = shapes.Shape0(filepath=filepath)

    def __repr__(self):
        return f"<TemperatureProfile [{self.create_nemesis_string()}]>"

    def _add_result(self, df):
        self.shape.data = df
        self.retrieved = True

    def _create_profile_from_previous_retrieval(prev_profile):
        raise NotImplementedError()

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
    def __init__(self, gas_name=None, gas_id=None, isotope_id=0, shape=None, **kwargs):
        """Create a profile for a given gas (optionally an isotopologue) with a given shape
        
        Args:
            gas_name: The name of the gas (eg 'CH4'). Specify either this OR gas_id
            gas_id: The radtrans ID of the gas (eg. 6). Specify either this OR gas_name
            isotope_id: The ID of the isotopologue to use. Use 0 for a mix of all at terrestrial abundance
            shape: A Shape object to use for the profile shape
        """
        
        super().__init__(**kwargs)
        self.isotope_id = isotope_id
        if shape is None:
            raise ValueError("shape attribute must be specified")
        self.shape = shape
        if not ((gas_id is None) ^ (gas_name is None)):
            raise ValueError("Specifiy exactly one of gas_name or gas_id (not both!)")
        if gas_name is None:
            self.gas_id = gas_id
            self.gas_name = constants.GASES.loc[constants.GASES.radtrans_id == gas_id].name.iloc[0]
        elif gas_id is None:
            self.gas_name = gas_name
            self.gas_id = constants.GASES.loc[constants.GASES.name == gas_name].radtrans_id.iloc[0]
        self.label = self.gas_name

    def __repr__(self):
        return f"<GasProfile {self.gas_name} [{self.create_nemesis_string()}]>"

    def _create_profile_from_previous_retrieval(prev_profile):
        return GasProfile(gas_name=prev_profile.gas_name, isotope_id=prev_profile.isotope_id, shape=prev_profile.shape, label=prev_profile.label)

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

    def get_name(self):
        if self.label is None:
            return f"{self.gas_name} {self.isotope_id}"
        else:
            return self.label


class AerosolProfile(Profile):
    def __init__(self, 
                 shape, 
                 radius, 
                 variance,
                 refractive_index, 
                 retrieve_optical=False, 
                 radius_error=None,
                 variance_error=None,
                 imag_refractive_index_error=None,
                 **kwargs):
        """Create a profile for a given aerosol with a given shape
        
        Args:
            shape: A Shape object to use for the profile shape
            label: (optional) An arbitrary label to associate with this profile"""
            
        super().__init__(**kwargs)
        self.shape = shape

        self.radius = radius
        self.variance = variance
        self.refractive_index = refractive_index
        self.retrieve_optical = retrieve_optical

        if retrieve_optical:
            self.radius_error = radius_error
            self.variance_error = variance_error
            self.imag_refractive_index_error = imag_refractive_index_error

    def __repr__(self):
        return f"<AerosolProfile {self.aerosol_id} [{self.create_nemesis_string()}]>"

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
            utils.write_nums(file, refwave, self.refractive_index.real)
            utils.write_nums(file, refwave)
            for line in lines:
                vals = line.split()
                if len(vals) < 2:
                    continue
                utils.write_nums(file, float(vals[0]), self.refractive_index.imag, self.imag_refractive_index_error)
    


    def get_name(self):
        if self.label is None:
            return f"Aerosol {self.aerosol_id}"
        else:
            return self.label
