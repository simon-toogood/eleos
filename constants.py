import pandas as pd


def _read_radtrans_gas_id():
    data = []
    with open("data/radtrans_ids.csv") as file:
        for line in file:
            x = line.rstrip("\n").split(",")
            data.append(x[:7] + [x[7:],])
    df = pd.DataFrame(data)
    df.columns = ["radtrans_id", "name", "H04", "G03", "H12", "H16", "N_iso", "isotopes"]
    return df


GASES = _read_radtrans_gas_id()

DISTANCES = {"jupiter": 5.2,
             "saturn": 9.546,
             "uranus": 19.2,
             "neptune": 30.0,
             "titan": 9.546,}