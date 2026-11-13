"""This module contains any constants used across the module. These attributes are not intended to be used from outside the module but may be useful in some circumstances."""

import pandas as pd
import os
from pathlib import Path


def _read_ch4_lines():
    df = pd.read_csv(PATH / "data/misc/ch4_emission.csv")
    df["wavelength"] = 10000 / df["nu"] # convert cm-1 to um because screw doing anything in wavenumbers
    df["intensity"] = df["sw"]
    df = df[["wavelength", "intensity"]]
    df.sort_values(by="intensity", ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df


def _get_nemesis_path():
    """Get the path of the NEMESIS installation. If your installation doesnt follow the structure
    NEMESIS/
        bin/
            gfort/
                Nemesis"
        lib/
        obj/
        radtrancode/
    then replace this function with the approriate version (directory name changes are fine, level changes are not)
    """

    r = os.popen('which Nemesis').read().strip()
    p = Path(r)
    return (p / "../../../").resolve()


PATH = Path(os.path.dirname(__file__)) #: The absolute path of the prepackaged data directory

NEMESIS_PATH = _get_nemesis_path()

GASES = pd.read_csv(PATH / "data/radtran_ids.csv") #: The database of gas names and IDs used by NEMESIS

CH4_LINES = _read_ch4_lines() #: The database of methane emission lines

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
                   "N2H4": 8} #: The IDs used by Makephase for looking up refractive indicies of the gases

GAS_LABELS = {
    'c2h4': r'$\mathrm{C_2H_4}$',
    'c4h2': r'$\mathrm{C_4H_2}$',
    'c3h4': r'$\mathrm{C_3H_4}$',
    'c6h6': r'$\mathrm{C_6H_6}$',
    'c3h8': r'$\mathrm{C_3H_8}$',
    'co2': r'$\mathrm{CO_2}$',
    'ch3': r'$\mathrm{CH_3}$',
    'h2': r'$\mathrm{H_2}$',
    'hcn': r'$\mathrm{HCN}$',
    'c2h6_all': r'$\mathrm{C_2H_6}$',
    'ch4_13C': r'$^{13}\mathrm{CH_4}$',
    'c2h2': r'$\mathrm{C_2H_2}$',
    'ch4': r'$\mathrm{CH_4}$',
    'co': r'$\mathrm{CO}$',
    'ash3': r'$\mathrm{AsH_3}$',
    'geh4': r'$\mathrm{GeH_4}$',
    'ch3d': r'$\mathrm{CH_3D}$',
    'nh3_all': r'$\mathrm{NH_3}$',
    'h2o': r'$\mathrm{H_2O}$',
    'ph3': r'$\mathrm{PH_3}$',
    'ch4 1': r'$\mathrm{CH_4}$',
    'ch4 2': r'$\mathrm{^{13}CH_4}$',
    'ch4 3': r'$\mathrm{^{12}CH_3D}$',
    'c2h4 0': r'$\mathrm{C_2H_4}$',
    'c4h2 0': r'$\mathrm{C_4H_2}$',
    'c3h4 0': r'$\mathrm{C_3H_4}$',
    'c6h6 0': r'$\mathrm{C_6H_6}$',
    'c3h8 0': r'$\mathrm{C_3H_8}$',
    'co2 0': r'$\mathrm{CO_2}$',
    'ch3 0': r'$\mathrm{CH_3}$',
    'h2 0': r'$\mathrm{H_2}$',
    'hcn 0': r'$\mathrm{HCN}$',
    'c2h6 0': r'$\mathrm{C_2H_6}$',
    'c2h2 0': r'$\mathrm{C_2H_2}$',
    'co 0': r'$\mathrm{CO}$',
    'ash3 0': r'$\mathrm{AsH_3}$',
    'geh4 0': r'$\mathrm{GeH_4}$',
    'ch3d 0': r'$\mathrm{CH_3D}$',
    'nh3 0': r'$\mathrm{NH_3}$',
    'h2o 0': r'$\mathrm{H_2O}$',
    'ph3 0': r'$\mathrm{PH_3}$',
}