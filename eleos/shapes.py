""" 
Each Shape object represents a differnet model in NEMESIS. Eg. ``Shape1`` represents profile shape 1 which is a constant VMR 
up to a knee pressure then a fractional scale height. For a ful list of possible models see the NEMESIS manual. Please note: 
due to the large number of profile shapes NEMESIS supports, 

Adding a Shape subclass
=======================

All subclasses of Shape must implement:
    - ID: Integer ID of the shape in NEMESIS
    - NAMES: The names of the parameters in order they appear in the mre file (do not include error parameters here)
    - __init__: Instantiates the parent class and defines:
        - nemesis_code: The number of the profile in the NEMESIS manual
        - Any parameters associated with the model (eg. knee pressure, fsh etc). These MUST be defined in the order that they appear in the .apr file
    - generate_apr_data: Returns a string of the profile shape parameters in the format
                         that the .apr file expects.

Additionally, some shapes require additional files in the core (eg. tempapr.dat). If this is the case,
override the create_required_files method.

Finally, add the new ShapeN class to the ALL_SHAPES list at the end of the file.

Class template::

    @shapeclass
    class ShapeN(Shape):
        "copy the desciption from the NEMESIS manual here (use triple quotes!)"
        ID: ClassVar[int] = N
        NAMES:  ClassVar[list[str]] = ["arg_name_1", "arg_name_2", ...]
        arg_name_1: type_1
        arg_name_2: type_2
        arg_name_2_error: type_2
        ...

        def generate_apr_data(self):
            return a string for apr file

        # optionally
        def create_required_files(self, directory):
            shutil.copy(..., directory)

"""

from dataclasses import dataclass
from typing import ClassVar
import shutil

from . import spx, utils


def shapeclass(*args, **kwargs):
    """Syntatic sugar for ``@dataclass(repr=False)``"""
    return dataclass(*args, **kwargs, repr=False)
    


class Shape:
    def __init__(self):
        """Do not instantiate directly, use a subclass"""
        pass

    def __repr__(self):
        out = f"Shape:\n    ID: {self.ID}"
        for name, value in self.__dict__.items():
            out += f"\n    {name}: {value}"    
        return out
    
    def generate_apr_data(self):
        """Generate the section of the .apr file that stores the shape parameters
        
        Args:
            None
            
        Returns:
            str: The string to write to the .apr file"""
        pass

    def create_required_files(self, directory):
        """Some Shapes require additional files to be created/copied into the core directory"""
        pass

@shapeclass
class Shape0(Shape):
    """Profile is to be treated as continuous over the pressure range of runname.ref, the
    next line of the .apr file should then contain a filename, which specifies the a
    priori profile as a function of height and should have the same number of levels
    as the .ref file."""
    ID: ClassVar[int] = 0
    NAMES: ClassVar[list[str]] = ["None"]
    filepath: str

    def create_required_files(self, directory):
        shutil.copy(self.filepath, directory)

    def generate_apr_data(self):
        return self.filepath.split("/")[-1]

@shapeclass
class Shape1(Shape):
    """Profile is to be represented as a deep VMR up to a certain ‘knee’ pressure, and
then a defined fractional scale height. The next line of the .apr file then contains
the ‘knee’ pressure, followed by the a priori deep and fractional scale height
values together with their estimated errors. """
    ID: ClassVar[int] = 1
    NAMES: ClassVar[list[str]] = ["deep_vmr", "fsh"]
    knee_pressure: float
    deep_vmr: float
    deep_vmr_error: float
    fsh: float
    fsh_error: float

    def generate_apr_data(self):
        return f"{self.knee_pressure}\n{self.deep_vmr} {self.deep_vmr_error}\n{self.fsh} {self.fsh_error}"

@shapeclass
class Shape2(Shape):
    """Profile is to be represented by a simple scaling of the corresponding profile
runname.ref (for T, v.m.r.), aerosol.ref (for aerosol density), parah2.ref (for
para-H2 fraction) or fcloud.ref (for fractional cloud cover). The next line of the
.apr file then contains the a priori factor and error."""
    ID: ClassVar[int] = 2
    NAMES: ClassVar[list[str]] = ["scale_factor"]
    scale_factor: float
    scale_factor_error: float

    def generate_apr_data(self):
        return f"{self.scale_factor} {self.scale_factor_error}"

@shapeclass
class Shape4(Shape):
    """Very similar to Shape1 in that the profile is to be
represented as a deep value up to a certain 'knee' pressure, and then a defined
fractional scale height. However, in this case the knee pressure is also a variable
parameter and thus must be supplied with an error estimate."""
    ID: ClassVar[int] = 4
    NAMES: ClassVar[list[str]] = ["knee_pressure", "deep_vmr", "fsh"]
    knee_pressure: float
    knee_pressure_error: float
    deep_vmr: float
    deep_vmr_error: float
    fsh: float
    fsh_error: float

    def generate_apr_data(self):
        return f"{self.knee_pressure} {self.knee_pressure_error}\n{self.deep_vmr} {self.deep_vmr_error}\n{self.fsh} {self.fsh_error}"

@shapeclass
class Shape32(Shape):
    """Similar to model 8 in that profile is a cloud profile represented by a variable
base pressure, specific density at the level and fractional scale height. The next
line of the .apr file then contains the a priori base pressure, followed by the a
priori opacity and fractional scale height values together with their estimated
errors. All quantities are taken as logs so negative fractional scale heights are
not allowed. Difference from Model 8 is that cloud density at pressures greater
than the base pressure is set to drop exponentially with increasing pressure with
a scale height of 1km, rather than just being set to zero. This makes it easier for
NEMESIS to actually find an optimal value of the knee pressure."""
    ID:  ClassVar[int] = 32
    NAMES: ClassVar[list[str]] = ["opacity", "fsh", "base_pressure"]
    base_pressure: float
    base_pressure_error: float
    opacity: float
    opacity_error: float
    fsh: float
    fsh_error: float

    def generate_apr_data(self):
        return f"{self.base_pressure} {self.base_pressure_error}\n{self.opacity} {self.opacity_error}\n{self.fsh} {self.fsh_error}"
        
@shapeclass
class Shape48(Shape):
    """As model 32, in that profile is a cloud profile represented by a variable base
pressure, specific density at the level, fractional scale height, but also a variable
top pressure. The next line of the .apr file then contains the a priori base
pressure, followed by the a priori top pressure, opacity and fractional scale
height values together with their estimated errors."""
    ID:  ClassVar[int] = 48
    NAMES: ClassVar[list[str]] = ["base_pressure", "top_pressure", "fsh", "opacity"]
    base_pressure: float
    base_pressure_error: float
    top_pressure: float
    top_pressure_error: float
    opacity: float
    opacity_error: float
    fsh: float
    fsh_error: float

    def generate_apr_data(self):
        return f"{self.base_pressure} {self.base_pressure_error}\n{self.top_pressure} {self.top_pressure_error}\n{self.opacity} {self.opacity_error}\n{self.fsh} {self.fsh_error}"
        
@shapeclass
class Shape444(Shape):
    """If VARIDENT(IVAR,1) is equal to 444 then the parameter described is a retrieval
of the imaginary part of a cloud’s complex refractive index spectrum. The cloud
particle identifier is given by VARIDENT(IVAR,2). The next line contains the
name of a separate input file, which contains the following information. Line 1
contains the mean radius of the particle size distribution and error (assumes standard
size distribution), while line 2 gives the variance of the size distribution and error.
Line 3 gives the number of wavelengths/wavenumbers for which imaginary
refractive index spectrum is tabulated, together with the correlation length of this a
priori spectrum. Line 4 gives a reference wavelength/wavenumber and the value of
the real part of the refractive index at that reference. Line 5 gives the
wavelength/wavenumber to which the extinction cross-section spectrum should be
normalised. Following lines contain the wavelengths/wavenumbers and the a priori
values of the imaginary refractive index spectrum and errors. In this model, the code
the real part of the refractive index spectrum is calculated with a Kramers-Kronig
analysis and then the Mie scattering properties of the particles calculated."""
    ID: ClassVar[int] = 444
    NAMES: ClassVar[list[str]] = ["radius", "variance", "refractive_index"]
    radius: float
    radius_error: float
    variance: float
    variance_error: float
    refractive_index: complex
    refractive_index_error: float

    def generate_apr_data(self):
        return f"cloudf{self.aerosol_id}.dat"

    def create_required_files(self, directory):
         # Get first wavelength
        wl = spx.read(directory+"nemesis.spx").geometries[0].wavelengths[0]

        # Get number of wavelengths in xsc file
        with open(directory+"nemesis.xsc") as file:
            lines = file.read().split("\n")
            nwave = (len(lines)-2) // 2
            refwave = float(lines[1].split()[0])
        
        # Generate cloudfN.dat
        with open(directory + f"cloudf{self.aerosol_id}.dat", mode="w+") as file:
            utils.write_nums(file, self.radius, self.radius_error)
            utils.write_nums(file, self.variance, self.variance_error)
            file.write(f"{nwave}    -1\n")
            utils.write_nums(file, refwave, self.refractive_index.real)
            utils.write_nums(file, refwave)
            for line in lines:
                vals = line.split()
                if len(vals) < 2:
                    continue
                print(float(vals[0]), repr(self.refractive_index.imag), repr(self.refractive_index_error))
                utils.write_nums(file, float(vals[0]), self.refractive_index.imag, self.refractive_index_error)
    

def get_shape_from_id(id_):
    """Given a shape ID integer, return a reference to the class corresponding to that ID. Note that this returns a class, not an instantiated object."""
    id_ = int(id_)
    for shape_class in ALL_SHAPES:
        if shape_class.ID == id_:
            return shape_class


ALL_SHAPES = [Shape0, Shape1, Shape2, Shape4, Shape32, Shape48, Shape444] #: A list of all the Shape classes