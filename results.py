import io
import itertools as it
import pandas as pd
import matplotlib.pyplot as plt

import profiles


class NemesisResult:
    def __init__(self, name, core_directory):
        self.name = name
        self.core_directory = core_directory
        self.profiles = []
        self.read_mre()

    def _parse_header_line(self, line, num_fields, cast_to):
        fields = [cast_to(x) for x in line.split()[:num_fields]]
        if num_fields == 1:
            return fields[0]
        else:
            return fields

    def read_mre(self):
        mre_file = self.core_directory + self.name + ".mre"
        with open(mre_file) as file:
            mre_data = file.read().split("\n")

        header = []
        block_ends = []
        for i, line in enumerate(mre_data):
            # Read the first 3 lines to the header array
            if i < 3:
                header.append(line)
            # Find the boundaries between result blocks
            elif "Variable" in line:
                block_ends.append(i)
        block_ends.append(i)

        # Set some attributes from the header info    
        self.num_retrievals = self._parse_header_line(header[0], num_fields=1, cast_to=int)
        self.ispec, self.ngeom, self.ny1, self.nx, self.ny2 = self._parse_header_line(header[1], num_fields=5, cast_to=int)
        self.latitude, self.longitude = self._parse_header_line(header[2], num_fields=2, cast_to=float)

        # Read in the data as a DataFrame
        self.fitted_spectrum = pd.read_table(mre_file, 
                                             names=["wavelength", "measured", "error", "pct_error", "model", "pct_diff"],
                                             index_col=0, sep="\s+", skiprows=5, nrows=block_ends[0]-6)

        # Read in each retrieved parameter
        for start_idx, end_idx in it.pairwise(block_ends):
            profile = profiles.read_profile_string(mre_data[start_idx+1])
            result = pd.read_table(io.StringIO("\n".join(mre_data[start_idx+3: end_idx])), sep="\s+")
            profile.add_result(result)
            self.profiles.append(profile)


    def plot_spectrum(self):
        fig, ax = plt.subplots(1, 1)
        ax.plot(self.fitted_spectrum.wavelength, self.fitted_spectrum.measured, c="k", lw=0.5, label="Measured")
        ax.plot(self.fitted_spectrum.wavelength, self.fitted_spectrum.model, c="r", lw=0.5, label="Model")
        ax.set_xlabel("Wavelength (μm)")
        ax.set_ylabel("Radiance (μW cm$^{-2}$ sr$^{-1}$ μm$^{-1}$)")
        ax.legend()
        fig.tight_layout()
        return fig, ax



