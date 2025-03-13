"""This module contains the classes for creating Profile objects. There are currently 3 profiles available through Eleos:
temperature (TemperatureProfile)
gases (GasProfile)
aerosols (AerosolProfile)"""

from collections import defaultdict
from pathlib import Path
import re
import numpy as np

from . import constants
from . import shapes
from . import utils
from . import results
from . import parsers
from . import spx


class Profile:
    def __init__(self, label=None):
        """This is the base class for all Profile objects. It should never be instantiated directly, use a subclass such as ``GasProfile`` or ``TemperatureProfile``"""
        self.retrieved = False
        self.core = None
        self.label = label

    def print_table(self, colors=True, **kwargs):
        if not self.retrieved:
            headers = ["Prior", "Error"]
        else:
            headers = ["Prior", "Error", "Retrieved", "Error", "Change"]

        names = self._get_displayable_attributes()
        data = []
        for name_group in names:
            data.append([])

            # Add prior (and retrieved) columns
            for i in range(4 if self.retrieved else 2):
                try:
                    data[-1].append(getattr(self, name_group[i]))
                except:
                    data[-1].append(" ")
        
            # Add the difference column
            if self.retrieved:
                try:
                    pct_diff = 100*(getattr(self, name_group[2]) - getattr(self, name_group[0])) / getattr(self, name_group[0])
                    end = "\x1b[0m"
                    if colors and (abs(pct_diff) > 100):
                        color = "\x1b[41m"
                    elif colors and ( 20 < abs(pct_diff) <= 100):
                        color = "\x1b[43m"
                    elif colors and (abs(pct_diff) <= 20):
                        color = "\x1b[42m"
                    else:
                        color = ""
                        end = ""
                    data[-1].append(f"{color}{pct_diff:+.1f}%{end}")
                except:
                    data[-1].append(" ")
        
        table = utils.generate_ascii_table(f"{self.label} {self.shape}", headers, [n[0] for n in names], data)
        print(table, **kwargs)
        return table

    def __setattr__(self, name, value):
        # Override the setattr method to allow sharing of attributes between the Profile and Shape objects for convenience
        if name == "NAMES":
            return super().__setattr__(name, value)

        try:
            if hasattr(self.shape, name):
                setattr(self.shape, name, value)
        except AttributeError:
            pass
        return super().__setattr__(name, value)

    @classmethod
    def from_previous_retrieval(cls, result, label):
        """Create a Profile object using the retrieved parameters from a previous retrieval as priors. Use either id or label to specify the profile to use in the previous retrieval.
        
        Args:
            result (NemesisResult): The result object from the previous retrieval
            label (str): The label of the profile to use in the previous retrieval"""

        prev = None
        for profile_label, profile in result.profiles.items():
            if profile_label == label and cls == type(profile):
                prev = profile
                break
            
        if prev is None:
            raise ValueError(f"Profile with label {label} not found in previous retrieval")
        
        new_profile = cls._create_profile_from_previous_retrieval(prev)
        new_profile.shape._set_prior_to_retrieved()
        new_profile.retrieved = False

        return new_profile

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
        """Create a temperature profile from a prior file. NOT IMPLEMENTED YET
        
        Args:
            filepath: The filepath of the prior temperature profile
            label (str): A label to associate with this Profile. By default it is "Temperature"
        """
        raise NotImplementedError("Temperature profiles are not fully implemented yet")
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
        """Create a profile for a given gas (optionally an isotopologue) with a given shape. Call signatures:

        Create profile with the name of the gas (ie 'NH3')
        GasProfile(gas_name, isotope_id, shape)

        Create profile with the ID of the gas (eg. 11)
        GasProfile(gas_id, isotope_id, shape)
        
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
        self.NAMES = self.shape.NAMES

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

        return [g for g in grouped if g[0] not in ("label", "gas_id", "isotope_id", "retrieved", "core", "shape", "NAMES")]

    def _create_profile_from_previous_retrieval(prev):
        return GasProfile(gas_name=prev.gas_name, 
                          isotope_id=prev.isotope_id, 
                          shape=prev.shape, 
                          label=prev.label)

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
                 n_lookup=None,
                 real_n=None,
                 imag_n=None, 
                 retrieve_optical=False, 
                 radius_error=None,
                 variance_error=None,
                 imag_n_error=None,
                 label=None):
        """Create a profile for an aerosol layer with a given shape and particle/optical properties. Call signatures:

        Constant refractive index over range (not retrieved):
        AerosolProfile(shape, radius, variance, real_n, imag_n)

        Constant refractive index over range (retrieved):
        AerosolProfile(shape, radius, radius_error, variance, variance_error, real_n, imag_n, imag_n_error)

        Use refractive index from lookup table (not retrieved):
        AerosolProfile(shape, radius, variance, n_lookup)

        Use refractive index from lookup table (retrieved):
        AerosolProfile(shape, radius, radius_error, variance, variance_error, n_lookup)
        
        Args:
            shape (Shape): A Shape object to use for the profile shape
            label (str): A label to associate with this Profile. By default it is "Aerosol <aerosol_id>" (eg. "Aerosol 1")       
            shape (Shape): A Shape object that describes the profile shape
            radius (float): Particle radius in microns
            variance (float): Varaince of the particle size distribution in microns (gamma distribution)
            n_lookup (str): If given, use the refractive indicies of this gas (eg. 'CH4')
            real_n (float): If n_lookup is not given, then use this as the real part of the refractive index
            imag_n (float): If n_lookup is not given, then use this as the imaginary part of the refractive index,
            retrieve_optical (bool): Whether to retrieve the optical porperties of the aerosol (this adds a 444 profile in NEMESIS)
            radius_error (float): If retrieve_optical is set, the error used for the particle radius prior
            variance_error (float): If retrieve_optical is set, the error used for the particle size variance prior
            imag_n_error (float): If retrieve_optical is set, the error used for the imaginary part of the refractive index
            label (str): A label to assosiate with this profile
        """

        super().__init__(label=label)
        if label is None:
            self.label = f"Aerosol {self.aerosol_id}"

        # Assign basic parameters
        self.aerosol_id = "UNASSIGNED"
        self.shape = shape
        self.shape.share_parameters(self)
        self.retrieve_optical = retrieve_optical

        # Set particle properties
        self.radius = radius
        self.variance = variance

        # Set either n_lookup or ral_n and imag_n
        if n_lookup is None:
            self.real_n = real_n
            self.imag_n = imag_n
            self.lookup = False
        else:
            self.n_lookup = n_lookup
            self.lookup = True

        # Set the prior errors if retrieving
        if retrieve_optical:
            self.radius_error = radius_error
            self.variance_error = variance_error
            self.imag_n_error = imag_n_error
            if radius_error is None or variance_error is None or imag_n_error is None:
                raise ValueError("Cannot retireve optical properties without specified errors. Did you remember to set radius_error, variance_error, or imag_n_error?")
        else:
            if radius_error is not None or variance_error is not None or imag_n_error is not None:
                raise ValueError("Cannot specify errors for radius/variance.imag_n without retrieving optical peroperties. Did you forget to set retrieve_optical=True?")

        self.NAMES = self.shape.NAMES + ["radius", "variance", "imag_n"]

    def __repr__(self):
        return f"<AerosolProfile {self.label} [{self.create_nemesis_string()}]>"

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

        return [g for g in grouped if g[0] not in ("label", "retrieved", "core", "aerosol_id", "NAMES", "retrieve_optical", "shape", "lookup")]

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

    def _create_profile_from_previous_retrieval(prev):
        if prev.retrieve_optical:
            if prev.lookup: 
                return AerosolProfile(shape=           prev.shape,
                                      radius=          prev.retrieved_radius,
                                      radius_error=    prev.retrieved_radius_error,
                                      variance=        prev.retrieved_variance,
                                      variance_error=  prev.retrieved_variance_error,
                                      n_lookup=        prev.n_lookup,
                                      label=           prev.label,
                                      retrieve_optical=True)
            else:
                return AerosolProfile(shape=           prev.shape,
                                      radius=          prev.retrieved_radius,
                                      radius_error=    prev.retrieved_radius_error,
                                      variance=        prev.retrieved_variance,
                                      variance_error=  prev.retrieved_variance_error,
                                      real_n=          prev.real_n,
                                      imag_n=          prev.retrieved_imag_n,
                                      imag_n_error=    prev.retrieved_imag_n_error,
                                      label=           prev.label,
                                      retrieve_optical=True)
        else:
            if prev.lookup:
                return AerosolProfile(shape=           prev.shape,
                                      radius=          prev.radius,
                                      variance=        prev.variance,
                                      n_lookup=        prev.n_lookup,
                                      label=           prev.label,
                                      retrieve_optical=False)

            else:
                return AerosolProfile(shape=           prev.shape,
                                      radius=          prev.radius,
                                      variance=        prev.variance,
                                      real_n=          prev.real_n,
                                      imag_n=          prev.imag_n,
                                      label=           prev.label,
                                      retrieve_optical=False)

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
        directory = Path(directory)

        # Get the wavelengths in the xsc file
        xsc_parser = parsers.NemesisXsc(directory / "nemesis.xsc")
        nwave = len(xsc_parser.xsc.wavelength)
        refwave = self.core.reference_wavelength

        # Get the refractive indicies in the Makephase output or from the given attributes
        if self.lookup:
            refindexes = parsers.MakephaseOut(directory / "makephase.out").data
            imag_ns = refindexes[self.label + " imag"]
            wi,_ = utils.find_nearest(refindexes.wavelength, refwave)
            ref_real_n = refindexes[self.label + " imag"].iloc[wi]
        else:
            imag_ns = [self.imag_n for x in range(nwave)]
            ref_real_n = self.real_n

        # Generate cloudfN.dat
        with open(directory / f"cloudf{self.aerosol_id}.dat", mode="w+") as file:
            utils.write_nums(file, self.radius, self.radius_error)
            utils.write_nums(file, self.variance, self.variance_error)
            file.write(f"{nwave}    {0 if self.lookup else -1}\n")
            utils.write_nums(file, refwave, ref_real_n)
            utils.write_nums(file, refwave)
            for imag_n, wavelength in zip(imag_ns, xsc_parser.xsc.wavelength):
                utils.write_nums(file, wavelength, imag_n, self.imag_n_error)
