"""This module contains the classes for creating Profile objects."""

from . import constants
from . import shapes
from . import utils


class Profile:
    def __init__(self):
        """This is the base class for all Profile objects. It should never be instantiated directly, use a subclass such as ``GasProfile`` or ``TemperatureProfile``"""
        self.retrieved = False

    def add_result(self, df):
        """Take in a DataFrame created by reading in the .mre file (results.NemesisResult.read_mre) and assign the correct attributes
        
        Args:
            df: pandas.DataFrame with columns "prior", "prior_error", "retrieved", "retrieved_error" and a row for each parameter
            
        Returns:
            None"""
        
        # Get order of parameters
        names = self.shape.NAMES
        if names == ["None"]:
            # In this case, it is a temperature retrieval
            # In the future if I add other profiles that also require a text file then this is going to break spectacularly
            self.shape.data = df
        else:
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


class TemperatureProfile(Profile):
    def __init__(self, filepath):
        """Create a temperature profile from a prior file
        
        Args:
            filepath: The filepath of the prior temperature profile
        """
        super().__init__()
        self.shape = shapes.Shape0(filepath=filepath)

    def __repr__(self):
        return f"<TemperatureProfile [{self.create_nemesis_string()}]>"

    def __str__(self):
        return f"TemperatureProfile:\n    File: {self.shape.filepath}\n    Retrieved: {self.retrieved}\n" + utils.indent(str(self.shape))

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


class GasProfile(Profile):
    def __init__(self, gas_name=None, gas_id=None, isotope_id=0, shape=None):
        """Create a profile for a given gas (optionally an isotopologue) with a given shape
        
        Args:
            gas_name: The name of the gas (eg 'CH4'). Specify either this OR gas_id
            gas_id: The radtrans ID of the gas (eg. 6). Specify either this OR gas_name
            isotope_id: The ID of the isotopologue to use. Use 0 for a mix of all at terrestrial abundance
            shape: A Shape object to use for the profile shape"""
        super().__init__()
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

    def __str__(self):
        return f"GasProfile:\n    Species: {self.gas_name}\n    Isotope: {self.isotope_id}\n" + utils.indent(str(self.shape))

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
    def __init__(self, aerosol_id=None, shape=None):
        """Create a profile for a given aerosol with a given shape
        
        Args:
            aerosol_id: The ID of the aerosol
            shape: A Shape object to use for the profile shape"""
        super().__init__()
        self.aerosol_id =abs(aerosol_id)
        self.shape = shape

    def __repr__(self):
        return f"<AerosolProfile {self.aerosol_id} [{self.create_nemesis_code()}]>"

    def __str__(self):
        return f"AerosolProfile:\n    ID: {self.aerosol_id}\n    Retrieved: {self.retrieved}\n" + utils.indent(str(self.shape))

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


def create_profile_from_code(code):
    """Create a Profile object from a NEMESIS code (eg. '23 0 1'). The Profile.shape attribute will be a refernce to a Shape subclass, 
    NOT an instantiated Shape[N] object as there is no a priori data. To create a fully instantiated Profile and Shape, see ``profiles.create_profile_from_apr``

    Args:
        code (str): The NEMESIS code in string form seperated by whitespace (eg. '23 0 1')

    Returns:
        Profile: The instantiated Profile object
    """
    tokens = [int(x) for x in code.split()]

    # Special case for temp profile
    if tokens == [0, 0, 0]:
        return TemperatureProfile(None)
    
    # Create the Shape
    shape = shapes.get_shape_from_id(tokens[2])

    # Create the Profile
    if tokens[0] > 0:
        return GasProfile(gas_id=tokens[0], isotope_id=tokens[1], shape=shape)
    elif tokens[0] < 0:
        return AerosolProfile(aerosol_id=abs(tokens[0]), shape=shape)
    else:
        raise NotImplementedError
    

def create_profile_from_apr(string):
    """Creates a Profile object based on the .apr string representation
    
    Args:
        string (str): The section of the .apr file that contains the parameters
        
    Returns:
        Profile: A fully instantiated Profile object"""
    lines = string.split("\n")
    code = lines[0].split(" - ")[0]
    profile = create_profile_from_code(code)
    params = [y for x in lines[1:] for y in x.split()]
    if isinstance(profile, TemperatureProfile):
        profile.filepath = params[0]
        return profile
    else:
        profile.shape = profile.shape(*params)
        return profile