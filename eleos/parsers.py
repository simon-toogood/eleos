"""This module provides parsing objects for reading some NEMESIS files, such as nemesis.ref"""

import pandas as pd
import itertools as it
import io
import numpy as np
from pathlib import Path

from . import utils
from . import constants
from . import profiles as profiles_


## TODO: Move all parsing routines here, .itr, .prc etc...


class NemesisRef:
    """Parser for nemesis.ref and nemesis.prf
    
    Attributes:
        amform:    Number of latitudes (always 1)
        planet_id: ID of the planet being analysed.
        latitude:  Latitude at which the .ref file applies
        num_layers:Number of layers used 
        num_gases: Number of gases in the file
        gas_names: List of the names of all the gases
        data:      Contains the pressure, temperature, and VMR profiles of each gas at each height
    """
    def __init__(self, filepath):
        self.filepath = Path(filepath)
        self._extra_header = True
        self.read()

    def read(self):
        with open(self.filepath) as file:
            lines = file.read().split("\n")
            if self._extra_header:
                del lines[1]
                
        self.amform, = utils.get_ints_from_string(lines[0])
        planet_id, latitude, num_layers, num_gases = utils.get_floats_from_string(lines[1])
        self.planet_id = int(planet_id)
        self.latitude = latitude
        self.num_layers = int(num_layers)
        self.num_gases = int(num_gases)

        self.gas_names = []
        for l in lines[2:2+int(self.num_gases)]:
            gas_id, isotope_id = utils.get_ints_from_string(l)
            gas_name = constants.GASES[constants.GASES.radtrans_id == gas_id].name.iloc[0]
            self.gas_names.append(f"{gas_name} {isotope_id}")

        self.data = pd.read_table(self.filepath, 
                                  skiprows=3+int(self._extra_header)+self.num_gases, 
                                  sep="\s+", 
                                  header=None)
        self.data.columns = ["height", "pressure", "temperature"] + self.gas_names


class NemesisPrf(NemesisRef):
    def __init__(self, filepath):
        self.filepath = Path(filepath)
        self._extra_header = False
        self.read()


class NemesisMre:
    """Parser for the nemesis.mre file
    
    Attributes:
        ispec (int): Don't know
        ngeom (int): Number of geometries (should be 1)
        latitude (float): Latitude of the observation
        longitude (float): Longitude of the observation
        retrieved_spectrum pd.DataFrame: DataFrame containing the measured spectrum + all error sources and the fitted model spectra and its errors
        retrieved_parameters List[pd.DataFrame]: List of DataFrames containing the retrieved parameters from each Profile"""
    def __init__(self, filepath):
        self.filepath = Path(filepath)
        self.read()

    def _parse_header_line(self, line, num_fields, cast_to):
        fields = [cast_to(x) for x in line.split()[:num_fields]]
        if num_fields == 1:
            return fields[0]
        else:
            return fields

    def read(self):
        with open(self.filepath) as file:
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
        self.ispec, self.ngeom, _,_,_ = self._parse_header_line(header[1], num_fields=5, cast_to=int)
        self.latitude, self.longitude = self._parse_header_line(header[2], num_fields=2, cast_to=float)

        # Read in the fitted spectrum as a DataFrame
        self.retrieved_spectrum = pd.read_table(self.filepath, 
                                                names=["wavelength", "measured", "error", "pct_error", "model", "pct_diff"],
                                                index_col=0, sep="\s+", skiprows=5, nrows=blocks[0]-7)

        # Read in each retrieved parameter 
        self.retrieved_parameters = []
        self.initial_state_vector = []
        with open(self.filepath) as file:
            for start, end in it.pairwise(blocks):
                data = utils.read_between_lines(file, start, end)
                df = pd.read_table(io.StringIO(data), skiprows=4, sep="\s+", names=["i", "ix", "prior", "prior_error", "retrieved", "retrieved_error"])
                df.drop(["i", "ix"], axis=1, inplace=True)
                self.initial_state_vector += list(df.prior)
                self.retrieved_parameters.append(df)


class NemesisXsc:
    """Parser for the nemesis.xsc file
    
    Attributes:
        xsc (pd.DataFrame): The aerosol cross-sections as a function of wavelength for each aerosol mode
        ssa (pd.DataFrame): The single scattering albedos as a function of wavelength for each aerosol modes"""
    
    def __init__(self, filepath):
        self.filepath = Path(filepath)
        self.read()

    def read(self):
        waves = []
        ssas = []
        xscs = []
        with open(self.filepath) as file:
            for i, line in enumerate(file):
                if i == 0:
                    continue
                if i % 2 == 1:
                    wavelength, *x = utils.get_floats_from_string(line)
                    waves.append(wavelength)
                    xscs.append(x)
                else:
                    s = utils.get_floats_from_string(line)
                    ssas.append(s)

        self.ssa = pd.DataFrame(ssas)
        self.ssa.insert(0, column="wavelength", value=waves)
        self.xsc = pd.DataFrame(xscs)
        self.xsc.insert(0, column="wavelength", value=waves)
            

class NemesisItr:
    """Parser for nemesis.itr. Also requires a NemesisMre parser
    
    Attributes:
        state_vectors: pd.DataFrame containing the linear state vector for each iteration"""
    
    def __init__(self, filepath, mre=None):
        self.filepath = Path(filepath)

        if mre is None:
            self.mre = NemesisMre(self.filepath.parent / "nemesis.mre")
        else:
            self.mre = mre

        self.read()

    def read(self):
        data = []
        count = -1
        with open(self.filepath) as file:
            for i, line in enumerate(file.read().split("\n")):
                if line == " ":
                    count = 3
                else:
                    count -= 1
                if count == 0:
                    d = []
                    for i, v in enumerate(line.split()):
                        d.append(float(v))
                    data.append(d)

        exps = []
        for i, value in enumerate(data[0]):
            flag = np.isclose(self.mre.initial_state_vector[i], np.exp(value))
            exps.append(flag)

        self.state_vectors = pd.DataFrame(data)
        for i, column in enumerate(self.state_vectors.columns):
            if exps[i]:
                self.state_vectors[column] = np.exp(self.state_vectors[column])

    def add_column_names(self, profiles):
        names = []
        for label, profile in profiles.items():
            names += [f"{label} {n}" for n in profile.shape.NAMES]
            if isinstance(profile, profiles_.AerosolProfile):
                names += [f"{label} radius", f"{label} variance", f"{label} imag_n"]
        self.state_vectors.columns = names


class NemesisPrc:
    def __init__(self, filepath):
        self.filepath = Path(filepath)
        self.read()
    
    def read(self):
        self.chisq = []
        with open(self.filepath) as file:
            for line in file:
                if "chi" in line and "should" not in line:
                    self.chisq.append(utils.get_floats_from_string(line)[0])


class MakephaseOut:
    def __init__(self, filepath):
        self.filepath = Path(filepath)
        self.read()
        self.add_aerosol_names()

    def read(self):
        xsc = NemesisXsc(Path(self.filepath).parent / "nemesis.xsc")
        wavelengths = xsc.xsc.wavelength

        with open(self.filepath) as file:
            isprev = False
            start = False
            out = []
            for line in file:
                if not isprev and not start:
                    out.append([])
                    start = True
                if line.startswith(" Refractive index:"):
                    isprev = True
                    start = False
                    out[-1].append(utils.get_floats_from_string(line))
                else:
                    isprev = False
        
        out = pd.DataFrame(np.hstack(out))
        out.insert(0, "wavelength", wavelengths)
        self.data = out

    def add_aerosol_names(self):
        # Get the names of the aerosol layers
        names = []
        with open(Path(self.filepath).parent / "aerosol_names.txt") as file:
            n = file.read().split("\n")
            for m in n:
                names.append(m + " real")
                names.append(m + " imag")
        
   
class AerosolPrf:
    """Parser for the aerosol.prf file
    
    Attributes:
        data: pd.DataFrame containing the aerosol density as a function of height. The units of aerosol density are particles per gram of atmosphere"""
    
    def __init__(self, filepath):
        self.filepath = Path(filepath)
        self.read()
        self.add_aerosol_names()

    def read(self):
        self.data = pd.read_table(self.filepath, sep="\s+", skiprows=2, header=None)
        num = len(self.data.columns) - 1
        header = ["height"] + [f"aerosol_{x}" for x in range(1, num+1)]
        self.data.columns = header

    def add_aerosol_names(self):
        with open(Path(self.filepath).parent / "aerosol_names.txt") as file:
            names = file.read().split("\n")
            self.data.columns = ["height"] + names

