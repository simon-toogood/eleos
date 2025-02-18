"""This module provides parsing objects for reading some NEMESIS files, such as nemesis.ref"""

import pandas as pd
import itertools as it
import io
from pathlib import Path

from . import utils
from . import constants

## TODO: Move all parsing routines here, .itr, .prc etc...


class NemesisRef:
    """Parser for nemesis.ref
    
    Attributes:
        amform:
        planet_id: """
    def __init__(self, filepath):
        self.filepath = Path(filepath)
        self._extra_header = True
        self.read()

    def read(self):
        with open(self.filepath) as file:
            lines = file.read().split("\n")
            if self._extra_header:
                del lines[1]
        
        print(lines[:10])
        
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
        retireved_spectrum pd.DataFrame: DataFrame containing the measured spectrum + all error sources and the fitted model spectra and its errors
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
        retrievals = []
        with open(self.filepath) as file:
            for start, end in it.pairwise(blocks):
                data = utils.read_between_lines(file, start, end)
                df = pd.read_table(io.StringIO(data), skiprows=4, sep="\s+", names=["i", "ix", "prior", "prior_error", "retrieved", "retrieved_error"])
                df.drop(["i", "ix"], axis=1, inplace=True)
                retrievals.append(df)
        self.retrieved_parameters = retrievals


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


class AerosolPrf:
    def __init__(self, filepath):
        self.filepath = Path(filepath)
        self.read()

    def read(self):
        self.data = pd.read_table(self.filepath, sep="\s+", skiprows=2)
        num = len(self.data.columns) - 1
        header = ["height"] + [f"aerosol_{x}" for x in range(1, num+1)]
        self.data.columns = header
