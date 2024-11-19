import io
import itertools as it
import pandas as pd
import matplotlib.pyplot as plt

import profiles
import utils
import shapes


class NemesisResult:
    def __init__(self, core_directory):
        """Class for storing the results of a NEMESIS retrieval"""
        self.core_directory = core_directory
        self.read_apr()
        self.read_mre()

    def _parse_header_line(self, line, num_fields, cast_to):
        """Take in a header line from the .mre file and split each token and cast to a dtype.
        Eg. "10     2     4     7.2" -> [10.0, 2.0, 4.0, 7.2]"""
        fields = [cast_to(x) for x in line.split()[:num_fields]]
        if num_fields == 1:
            return fields[0]
        else:
            return fields
        
    def read_apr(self):
        """Read in the nemesis.apr and create a list of Profile objects as an object attribute
        
        Args:
            None
            
        Returns:
            None"""
        with open(self.core_directory+"nemesis.apr", mode="r") as file:
            blocks = []
            # First read of file: get line numbers of block starts
            for i, line in enumerate(file):
                if i == 0:
                    pass
                elif i == 1:
                    num_profiles = int(line)
                elif " - " in line:
                    blocks.append(i)
            blocks.append(i+1)

            # Second read: get each profile from the block
            self.profiles = []
            for start, end in it.pairwise(blocks):
                data = utils.read_between_lines(file, start, end)
                profile = profiles.create_profile_from_apr(data)
                self.profiles.append(profile)
        
    def read_mre(self):
        """Read in the nemesis.mre file and """
        mre_file = self.core_directory + "nemesis.mre"
        with open(mre_file) as file:
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
        self.num_retrievals = self._parse_header_line(header[0], num_fields=1, cast_to=int)
        self.ispec, self.ngeom, self.ny1, self.nx, self.ny2 = self._parse_header_line(header[1], num_fields=5, cast_to=int)
        self.latitude, self.longitude = self._parse_header_line(header[2], num_fields=2, cast_to=float)

        # Read in the data as a DataFrame
        self.fitted_spectrum = pd.read_table(mre_file, 
                                             names=["wavelength", "measured", "error", "pct_error", "model", "pct_diff"],
                                             index_col=0, sep="\s+", skiprows=5, nrows=blocks[0]-7)

        # Read in each retrieved parameter
        with open(mre_file) as file:
            for profile, (start, end) in zip(self.profiles, it.pairwise(blocks)):
                data = utils.read_between_lines(file, start, end)
                df = pd.read_table(io.StringIO(data), skiprows=3, sep="\s+", names=["i", "value", "err", ])
                print(df)
                #profile.shape.result_matrix = data

# index  nothing  prior error  retrieved error
    def plot_spectrum(self):
        fig, ax = plt.subplots(1, 1)
        ax.plot(self.fitted_spectrum.wavelength, self.fitted_spectrum.measured, c="k", lw=0.5, label="Measured")
        ax.plot(self.fitted_spectrum.wavelength, self.fitted_spectrum.model, c="r", lw=0.5, label="Model")
        ax.set_yscale("log")
        ax.set_xlabel("Wavelength (μm)")
        ax.set_ylabel("Radiance (μW cm$^{-2}$ sr$^{-1}$ μm$^{-1}$)")
        ax.legend()
        fig.tight_layout()
        return fig, ax



