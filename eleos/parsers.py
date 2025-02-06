"""This module provides parsing objects for reading some NEMESIS files, such as nemesis.ref"""

import pandas as pd
import itertools as it
import io

from . import utils
from . import constants

## TODO: Move all parsing routines here, .itr, .prc etc...


class NemesisRef:
    def __init__(self, filepath):
        self.filepath = filepath
        self.read()

    def read(self):
        with open(self.filepath) as file:
            lines = file.read().split("\n")
        
        self.amform, = utils.get_ints_from_string(lines[0])
        planet_id, latitude, num_layers, num_gases = utils.get_floats_from_string(lines[2])
        self.planet_id = int(planet_id)
        self.latitude = latitude
        self.num_layers = int(num_layers)
        self.num_gases = int(num_gases)

        self.gas_names = []
        for l in lines[3:3+int(self.num_gases)]:
            gas_id, isotope_id = utils.get_ints_from_string(l)
            gas_name = constants.GASES[constants.GASES.radtrans_id == gas_id].name.iloc[0]
            self.gas_names.append(f"{gas_name} {isotope_id}")

        self.data = pd.read_table(self.filepath, skiprows=4+self.num_gases, sep="\s+", header=None)
        self.data.columns = ["height", "pressure", "temperature"] + self.gas_names

    def write(self, filepath=None):
        if filepath is None:
            filepath = self.filepath

        with open(filepath, mode="w+") as file:
            file.write(str(self.amform) + "\n")
            file.write("1\n")
            file.write(f"{self.planet_id} {self.latitude:.2f} {self.num_layers} {self.num_gases}\n")
            
            out = self.data.copy()
            for gas_name in self.gas_names:
                name, isotope_id = gas_name.split(" ")
                gas_id = constants.GASES[constants.GASES.name == name].radtrans_id.iloc[0]
                file.write(f"{gas_id} {isotope_id}\n")
                out[gas_name] = self.data[gas_name].apply(lambda x: f"{x:.5e}")

            file.write(out.to_string(index=False, col_space=13))

    def set_pressure_limits(self, pmin=-float("inf"), pmax=float("inf")):
        self.data = self.data[(self.data.pressure > pmin) & (self.data.pressure < pmax)]
        self.num_layers = len(self.data)

    def remove_gas(self, gas_name):
        raise NotImplementedError()
    

class NemesisMre:
    def __init__(self, filepath):
        self.filepath = filepath
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

