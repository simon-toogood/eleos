"""This module contains the classes for creating Profile objects."""

from itertools import zip_longest
from collections import defaultdict

from . import constants
from . import shapes
from . import utils


class Profile:
    def __init__(self, label=None):
        """This is the base class for all Profile objects. It should never be instantiated directly, use a subclass such as ``GasProfile`` or ``TemperatureProfile``"""
        self.retrieved = False
        self.core = None
        self.label = label

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

    def __str__(self):
        # generate table title
        try:
            self.label
        except:
            self.label = None
        if self.label is not None:
            title = self.label
        elif isinstance(self, GasProfile):
            title = self.gas_name
        elif isinstance(self, (AerosolProfile, ImagRefractiveIndexProfile)):
            title = self.__class__.__name__ + f" ID={self.aerosol_id}"
        else:
            title = self.__class__.__name__

        # get shape attributes
        attrs = self.shape.__dict__.keys()

        # group attributes together by parameter
        grouped = defaultdict(list)
        for attr in attrs:
            base = attr.removeprefix("retrieved_").removesuffix("_error")
            grouped[base].append(attr)

        # extract all the attribute values
        values = []
        for name, attrs in grouped.items():
            values.append([])
            for i, attr in zip_longest(range(4), attrs):
                if attr is not None:
                    values[-1].append(getattr(self.shape, attr))
                else:
                    values[-1].append("NA")
    
        return utils.generate_ascii_table(title, ["Prior", "Error", "Retrieved", "Error"], grouped.keys(), values)


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

    def _add_result(self, df):
        self.shape.data = df
        self.retrieved = True


class GasProfile(Profile):
    def __init__(self, gas_name=None, gas_id=None, isotope_id=0, shape=None, **kwargs):
        """Create a profile for a given gas (optionally an isotopologue) with a given shape
        
        Args:
            gas_name: The name of the gas (eg 'CH4'). Specify either this OR gas_id
            gas_id: The radtrans ID of the gas (eg. 6). Specify either this OR gas_name
            isotope_id: The ID of the isotopologue to use. Use 0 for a mix of all at terrestrial abundance
            shape: A Shape object to use for the profile shape
            label: (optional) An arbitrary label to associate with this profile"""
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

    def __repr__(self):
        return f"<GasProfile {self.gas_name} [{self.create_nemesis_code()}]>"

    def create_nemesis_code(self):
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
        return self.create_nemesis_code() + " - " + self.gas_name + "\n" + self.shape.generate_apr_data()


class AerosolProfile(Profile):
    def __init__(self, shape, **kwargs):
        """Create a profile for a given aerosol with a given shape
        
        Args:
            aerosol_id: The ID of the aerosol
            shape: A Shape object to use for the profile shape
            label: (optional) An arbitrary label to associate with this profile"""
        super().__init__(**kwargs)
        self.shape = shape

    def __repr__(self):
        return f"<AerosolProfile {self.aerosol_id} [{self.create_nemesis_code()}]>"

    def create_nemesis_code(self):
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
        return self.create_nemesis_code() + f" - CP{self.aerosol_id}\n" + self.shape.generate_apr_data()


class ImagRefractiveIndexProfile(Profile):
    def __init__(self, shape, **kwargs):
        """Create a profile that retrieves the imaginary part of a clouds refractive index. Will require a corresponding AerosolProfile
        to be defined. This is currently limited to constant refractive index as a function of wavelength.

        Args:
            aerosol_id (int): The ID of the aerosol
            shape (Shape444): The associated Shape444 object containing the aerosol properties
            label: (optional) An arbitrary label to associate with this profile
            """

        super().__init__()
        assert isinstance(shape, shapes.Shape444)
        self.shape = shape

    def __repr__(self):
        return f"<ImagRefractiveIndexProfile [{self.create_nemesis_string()}]>"

    def create_nemesis_string(self):
        """Create the NEMESIS code that represents the profile.
        
        Args:
            None
            
        Returns:
            str: The NEMESIS string"""
        return f"444 {self.aerosol_id} 444"
    
    def generate_apr_data(self):
        return f"{self.create_nemesis_string()} - n{self.aerosol_id}\n{self.shape.generate_apr_data()}"
