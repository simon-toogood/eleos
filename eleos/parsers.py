"""This module provides parsing objects for reading some NEMESIS files, such as nemesis.ref"""

import pandas as pd
import itertools as it
import io
import numpy as np
from pathlib import Path
import astropy.units as u
import struct

import warnings
warnings.formatwarning = lambda msg, *_: f"Warning: {msg}\n"

from . import utils
from . import constants
from . import profiles as profiles_


## TODO: Move all parsing routines here, .itr, .prc etc...
## I also want to be ablet o completely reconstruct the core object from these files so no pickling of NemesisCore needed


class Parser:
    def __init__(self, filepath):
        self.filepath = Path(filepath)
        self.read()


class NemesisRef(Parser):
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
        self._extra_header = True
        super().__init__(filepath)

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


class NemesisMre(Parser):
    """Parser for the nemesis.mre file
    
    Attributes:
        ispec (int): Don't know
        ngeom (int): Number of geometries (should be 1)
        latitude (float): Latitude of the observation
        longitude (float): Longitude of the observation
        retrieved_spectrum pd.DataFrame: DataFrame containing the measured spectrum + all error sources and the fitted model spectra and its errors
        retrieved_parameters List[pd.DataFrame]: List of DataFrames containing the retrieved parameters from each Profile"""
    
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


class NemesisXsc(Parser):
    """Parser for the nemesis.xsc file
    
    Attributes:
        xsc (pd.DataFrame): The aerosol cross-sections as a function of wavelength for each aerosol mode
        ssa (pd.DataFrame): The single scattering albedos as a function of wavelength for each aerosol modes"""
    
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
            

class NemesisItr(Parser):
    """Parser for nemesis.itr. Also requires a NemesisMre parser
    
    Attributes:
        state_vectors (pd.DataFrame): Linear state vectors for each iteration"""
    
    def __init__(self, filepath, mre=None):
        if mre is None:
            self.mre = NemesisMre(Path(filepath).parent / "nemesis.mre")
        else:
            self.mre = mre
        super().__init__(filepath)

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
            names += [f"{label} {n}" for n in profile.VARIABLES]
        self.state_vectors.columns = names


class NemesisPrc(Parser):
    """Parser for nemesis.prc
    
    Attributes:
        chisq (List[float]): List of chi-squared values"""
    def read(self):
        self.chisq = []
        with open(self.filepath) as file:
            for line in file:
                if "chi" in line and "should" not in line:
                    self.chisq.append(utils.get_floats_from_string(line)[0])


class NemesisSpx(Parser):
    """Parser for the nemesis.spx file

    Attributes:
        lat (float):              Latitude of the observation
        lon (float):              longitude of the observation
        phase (float):            Phase angle of the observation
        emission (float):         Emission angle of the observation
        azimuth (float):          Azimuth angle of the observation
        wgeom (float):            wgeom
        nconv (int):              nconv
        nav (int):                nav
        fwhm (float):             FWHM of the instrument
        wavelengths (np.ndarray): Wavelengths of the observed spectra in um
        spectrum (np.ndarray):    The observed spectra in uW/cm2/sr/um
        error (np.ndarray):       The error on the observed spectra in uW/cm2/sr/um
    """

    def read(self):
        """Read the spx file in
        
        Args:
            None
            
        Returns:
            None
        """
        with open(self.filepath) as f:
            self.fwhm, *_ = (float(x) for x in f.readline().split())
            self.nconv = int(float(f.readline()))
            self.nav = int(float(f.readline()))
            self.lat, self.lon, self.phase, self.emission, self.azimuth, self.wgeom = (
                float(x) for x in f.readline().split()
            )
            self.wavelengths = np.full(self.nconv, np.nan)
            self.spectrum = np.full(self.nconv, np.nan)
            self.error = np.full(self.nconv, np.nan)
            for i in range(self.nconv):
                w, s, e = (float(x) for x in f.readline().split())
                self.wavelengths[i] = w
                self.spectrum[i] = s
                self.error[i] = e

    def write(self, filepath, exclude_nans=True):
        """
        Write an SPX file to disk

        Args:
            filepath (str): The filepath of the new .spx file

        Returns:
            None
        """

        lines = []
        lines.append(f'{self.fwhm}\t{self.lat}\t{self.lon}\t1')
        if exclude_nans:
            mask = ~(np.isnan(self.spectrum) | np.isnan(self.error))
        else:
            mask = np.ones_like(self.spectrum).astype(bool)
        nconv = len(self.wavelengths[mask])
        lines.append(f'{nconv}')
        lines.append(f'{self.nav}')
        lines.append(f'{self.lat}\t{self.lon}\t{self.phase}\t{self.emission}\t{self.azimuth}\t{self.wgeom}')

        if np.max(self.spectrum) > 1e-4:
            warnings.warn(f"This spectrum looks extremely bright (max value {np.max(self.spectrum)}) - have you converted from MJy/sr to W/cm2/sr/um?")

        for w, s, e in zip(self.wavelengths[mask], self.spectrum[mask], self.error[mask]):
            lines.append(f'{w: 12f}  {s: 12e}  {e: 12e}')
        lines.append('')  # end file with a newline

        with open(filepath, 'w') as f:
            f.write('\n'.join(lines))

    @property
    def wavelength(self):
        return self.wavelengths
    
    @wavelength.setter
    def wavelength(self, value):
        self.wavelengths = np.array(value)
    
    @property
    def errors(self):
        return self.error
    
    @errors.setter
    def errors(self, value):
        self.error = np.array(value)

    @classmethod
    def from_lists(cls,
                   wavelengths, 
                   spectrum, 
                   error, 
                   lat, lon, phase, emission, azimuth, fwhm=0,
                   convert=True):
        spx = cls.__new__(cls)
        if convert:
            spx.spectrum = NemesisSpx.convert_MJysr_to_Wcm2srum(wavelengths, spectrum)
            spx.error = NemesisSpx.convert_MJysr_to_Wcm2srum(wavelengths, error)
        else:
            spx.spectrum = np.array(spectrum)
            spx.error = np.array(error)
        spx.wavelengths = np.array(wavelengths)
        spx.lat = lat
        spx.lon = lon
        spx.phase = phase
        spx.emission = emission
        spx.azimuth = azimuth
        spx.fwhm = fwhm
        spx.nconv = len(spx.wavelengths)
        spx.nav = 1
        spx.wgeom = 1
        return spx

    @staticmethod
    def convert_MJysr_to_Wcm2srum(
        wavelengths: np.ndarray,
        spectrum: np.ndarray,
    ) -> np.ndarray:
        """
        Convert a spectrum or cube from MJy/sr to W/cm2/sr/micron.
        """
        while len(wavelengths.shape) < len(spectrum.shape):
            # work for cubes passed to sp
            wavelengths = np.expand_dims(wavelengths, -1)

        spx_MJysr = spectrum * u.MJy / u.sr
        spx_Wcm2srum = spx_MJysr.to(
            u.W / (u.cm * u.cm) / u.sr / u.micron,
            equivalencies=u.spectral_density(wavelengths * u.micron),
        )

        return spx_Wcm2srum.value

    @staticmethod
    def convert_Wcm2srum_to_MJysr(
        wavelengths: np.ndarray,
        spectrum: np.ndarray,
    ) -> np.ndarray:
        """
        Convert a spectrum or cube from W/cm2/sr/micron to MJy/sr.
        """
        while len(wavelengths.shape) < len(spectrum.shape):
            # work for cubes passed to sp
            wavelengths = np.expand_dims(wavelengths, -1)

        spx_Wcm2srum = spectrum * u.W / (u.cm * u.cm) / u.sr / u.micron
        spx_MJysr = spx_Wcm2srum.to(
            u.MJy / u.sr,
            equivalencies=u.spectral_density(wavelengths * u.micron),
        )
        return spx_MJysr.value


class NemesisInp(Parser):   
    """Parser for nemesis.inp
    
    Attributes:
        wavelength (bool):      Whether to use wavelength or wavenumber
        scattering (bool):      Whether to use a scattering run
        linebyline (bool):      Whether to use line-by-line or corrrelated-k method
        woff (float):           Wavelength offset
        fmerror_path (str):     Path to the file contianing forward modelling errors for each wavelength
        num_iterations (int):   Maximum number of iterations 
        min_phi_change (float): Minimim percentage change in phi before terminating
        num_spectra (int):      Number of spectra to retrieve
        start_spectra (int):    ID of spectra to start with
        sub_previous (bool):    Whether to use a previous retrieval
        output_format (int):    Format of the output files
    """
    
    def read(self):
        with open(self.filepath) as file:
            lines = file.read().split("\n")
            self.wavelength, self.scattering, self.linebyline = map(bool, utils.get_ints_from_string(lines[0]))
            self.woff, = utils.get_floats_from_string(lines[1])
            self.fmerror_path = lines[2]
            self.num_iterations, = utils.get_ints_from_string(lines[3])
            self.min_phi_change, = utils.get_floats_from_string(lines[4])
            self.num_spectra, self.start_spectra = utils.get_ints_from_string(lines[5])
            self.sub_previous, = map(bool, utils.get_ints_from_string(lines[6]))
            self.output_format, = utils.get_ints_from_string(lines[7])


class MakephaseOut(Parser):
    def __init__(self, filepath):
        super().__init__(filepath)
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
        
   
class AerosolPrf(Parser):
    """Parser for the aerosol.prf file
    
    Attributes:
        data: pd.DataFrame containing the aerosol density as a function of height. The units of aerosol density are particles per gram of atmosphere"""
    
    def __init__(self, filepath):
        super().__init__(filepath)
        self.add_aerosol_names()

    def read(self):
        self.data = pd.read_table(self.filepath, sep="\s+", skiprows=2, header=None)
        num = len(self.data.columns) - 1
        header = ["height"] + [f"aerosol_{x}" for x in range(1, num+1)]
        self.data.columns = header

    def add_aerosol_names(self):
        with open(Path(self.filepath).parent / "aerosol_names.txt") as file:
            names = file.read().split("\n")
            if names == [""]:
                return
            self.data.columns = ["height"] + names


class kTable(Parser):
    """Parser for ktables
    
    Atributes:
        ktable (np.ndarray):      4D array of absorption coefficients with axes (wavelength, pressure, temperature, g-ordinate).
        wavelength (np.ndarray):  Array of spectral points
        pressure (np.ndarray):    Array of pressure levels
        temperature (np.ndarray): Array of temperature levels
        g_ordinates (np.ndarray): Array of g-ordinates used in correlated-k integration
        g_weights (np.ndarray):   Corresponding quadrature weights for g-ordinates
        irec0 (int):              Number of header records
        npoint (int):             Number of spectral points
        vmin (float):             Starting value of the spectral axis
        delv (float):             Spectral spacing. If negative, grid is non-uniform and read from file.
        fwhm (float):             Full Width at Half Maximum for spectral bins
        np_ (int):                Number of pressure levels
        nt (int):                 Number of temperature levels
        ng (int):                 Number of g-ordinates
        idgas (int):              Identifier for the gas species 
        isogas (int):             Identifier for the isotope of the gas
    """
    
    @staticmethod
    def read_float(file):
        return struct.unpack('f', file.read(4))[0]

    @staticmethod
    def read_long(file):
        return struct.unpack('i', file.read(4))[0]

    def read(self):
        with open(self.filepath, 'rb') as f:

            # Read header
            self.irec0  = self.read_long(f)
            self.npoint = self.read_long(f)
            self.vmin   = self.read_float(f)
            self.delv   = self.read_float(f)
            self.fwhm   = self.read_float(f)
            self.np_    = self.read_long(f)
            self.nt     = self.read_long(f)
            self.ng     = self.read_long(f)
            self.idgas  = self.read_long(f)
            self.isogas = self.read_long(f)

            # Read g ordinates and weights
            self.g_ordinates = np.array([self.read_float(f) for _ in range(self.ng)])
            self.g_weights = np.array([self.read_float(f) for _ in range(self.ng)])

            # Skip two records
            _ = self.read_float(f)
            _ = self.read_float(f)

            # Construct pressure/temperature grid
            self.pressure = np.array([self.read_float(f) for _ in range(self.np_)])
            self.temperature = np.array([self.read_float(f) for _ in range(self.nt)])

            # Read in wavelengths
            nrec = 10 + 2 * self.ng + 2 + self.np_ + self.nt
            if self.delv < 0.0:
                self.wavelength = np.array([self.read_float(f) for _ in range(self.npoint)])
                nrec += self.npoint
            else:
                self.wavelength = self.vmin + self.delv * np.arange(self.npoint)

            # Skip records to reach k-table data
            jrec = self.irec0 - nrec - 1
            for _ in range(jrec):
                _ = self.read_float(f)

            # Read k-table
            total_values = self.npoint * self.np_ * self.nt * self.ng
            ktable_flat = np.frombuffer(f.read(total_values * 4), dtype=np.float32)
            self.ktable = ktable_flat.reshape((self.npoint, self.np_, self.nt, self.ng))

