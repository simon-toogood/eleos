"""This module contains any constants used across the module. These attributes are not intended to be used from outside the module but may be useful in some circumstances."""

import pandas as pd
import os
from pathlib import Path


def _read_radtrans_gas_id():
    data = []
    with open(PATH / "data/radtrans_ids.csv") as file:
        for line in file:
            x = line.rstrip("\n").split(",")
            data.append(x[:7] + [x[7:],])
    df = pd.DataFrame(data)
    df.columns = ["radtrans_id", "name", "H04", "G03", "H12", "H16", "N_iso", "isotopes"]
    df["G03"].astype(int)
    df = df.astype({"radtrans_id":int, "name": str, "H04": int, "G03": int, "H12": int, "H16": int, "N_iso": int})
    return df

PATH = Path(os.path.dirname(__file__)) #: The absolute path of the prepackaged data directory

GASES = _read_radtrans_gas_id() #: The database of gas names and IDs used by NEMESIS

DISTANCE = {"jupiter": 5.2,
             "saturn": 9.546,
             "uranus": 19.2,
             "neptune": 30.0,
             "titan": 9.546,} #: The distances from the Sun for all the bodies with significant atmospheres

GRAVITY = {"jupiter": 24.79, 
           "saturn": 10.44, 
           "uranus": 8.69, 
           "neptune": 11.15, 
           "titan":1.352} #: The surface gravity for all the bodies with significant atmospheres

MAKEPHASE_GASES = {"H2O": 3,
                   "NH3": 4,
                   "Tholins": 5,
                   "CH4": 6,
                   "NH4SH": 7,
                   "N2H4": 8,
                   "H2SO4": 9} #: The IDs used by Makephase for looking up refractive indicies of the gases