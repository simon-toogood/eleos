import constants
import shapes


class Profile:
    def __init__(self):
        """Do not instantiate directly, use a subclass"""
        pass

    def add_result(self, result):
        self.result = result


class TemperatureProfile(Profile):
    def __init__(self, filepath):
        """Create a temperature profile from a prior file
        
        Args:
            filepath: The filepath of the prior temperature profile
        """
        super().__init__()
        self.shape = shapes.Shape0(filepath=filepath)

    def __str__(self):
        return f"TemperatureProfile[{self.create_nemesis_string()}]"

    def create_nemesis_string(self):
        return f"0 0 0"
    
    def generate_apr_data(self):
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

    def __str__(self):
        return f"GasProfile({self.name}) [{self.create_nemesis_string()}]"

    def create_nemesis_string(self):
        return f"{self.gas_id} {self.isotope_id} {self.shape.ID}"
    
    def generate_apr_data(self, *args):
        return self.create_nemesis_string() + " - " + self.gas_name + "\n" + self.shape.generate_apr_data()


class AerosolProfile(Profile):
    def __init__(self, aerosol_id=None, shape=None):
        """Create a profile for a given aerosol with a given shape
        
        Args:
            aerosol_id: The ID of the aerosol
            shape: A Shape object to use for the profile shape"""
        super().__init__()
        self.aerosol_id =abs(aerosol_id)
        self.shape = shape

    def __str__(self):
        return f"AerosolProfile([{self.create_nemesis_string()}])"
    
    def create_nemesis_string(self):
        return f"-{self.aerosol_id} 0 {self.shape.ID}"
    
    def generate_apr_data(self):
        return self.create_nemesis_string() + f" - CP{self.aerosol_id}\n" + self.shape.generate_apr_data()


def create_profile_from_code(code):
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
        string: The section of the .apr file that contains the parameters
        
    Returns:
        Profile: A Profile-subclassed object"""
    lines = string.split("\n")
    code = lines[0].split(" - ")[0]
    profile = create_profile_from_code(code)
    params = [y for x in lines[1:] for y in x.split()]
    if isinstance(profile, TemperatureProfile):
        profile.filepath = params[0]
        return profile
    else:
        profile.shape.__init__(profile.shape, *params)
        return profile