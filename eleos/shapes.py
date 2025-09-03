""" 
Each Shape object represents a differnet model in NEMESIS. Eg. ``Shape1`` represents profile shape 1 which is a constant VMR 
up to a knee pressure then a fractional scale height. For a ful list of possible models see the NEMESIS manual. Please note: 
due to the large number of profile shapes NEMESIS supports, 

Adding a Shape subclass
=======================

All subclasses of Shape must implement:
    - ID: Integer ID of the shape in NEMESIS
    - CONSTANTS: The names of parameters that are not fitted by NEMESIS but required for the model
    - VARIABLES: The names of the parameters in order they appear in the mre file (do not include error parameters here). **These MUST be defined in the order that they appear in the .apr file!!**
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
        CONSTANTS: ClassVar[list[str]] = ["arg_name_1", ...]
        VARIABLES: ClassVar[list[str]] = ["arg_name_2", "arg_name_3", ...]
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

from . import parsers
from . import utils


def shapeclass(*args, **kwargs):
    """Syntatic sugar for ``@dataclass(repr=False)``"""
    return dataclass(*args, **kwargs, repr=False)
    


class Shape:
    def __init__(self):
        """Do not instantiate directly, use a subclass"""
        pass

    def __repr__(self):
        return f"<Shape{self.ID}>"
    
    def _set_prior_to_retrieved(self):
        """Set the prior values to the retrieved values and delete the retrived_* attributes."""
        for base in self.VARIABLES:
            for name in [base, f"{base}_error"]:
                try:
                    setattr(self, name, getattr(self, "retrieved_"+name))
                    delattr(self, "retrieved_"+name)
                except AttributeError:
                    print(f"Clear attribute failed for {name} on {self}")

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
    as the .ref file.
    
    In Eleos, one can specify either a filepath or give the values and error to be used directly.

    Parameters
    ----------
    filepath : str, optional
        Path to the input profile file.
    values : list[float], optional
        Profile values to be written if no file is given.
    errors : list[float], optional
        Error values corresponding to `values`.
    """
    ID: ClassVar[int] = 0
    CONSTANTS: ClassVar[list[str]] = []
    VARIABLES: ClassVar[list[str]] = []
    filepath: str = None
    values: list[float] = None
    errors: list[float] = None

    def create_required_files(self, directory):
        if self.filepath is not None:
            shutil.copy(self.filepath, directory)
            self.values = []
            self.errors = []
            with open(self.filepath) as file:
                lines = file.read().split("\n")[1:]
                for line in lines:
                    p, v, e = utils.get_floats_from_string(line)
                    self.values.append(v)
                    self.errors.append(e)
        else:
            self.filepath = directory / "shape0.dat"
            ref = parsers.NemesisRef(directory / "nemesis.ref")
            with open(self.filepath, mode="w+") as file:
                file.write(f"{len(ref.data.pressure)} 1\n#")
                for p, v, e in zip(ref.data.pressure, self.values, self.errors):
                    file.write(f"{p:08e}  {v:08e}  {e:08e}\n")

    def generate_apr_data(self):
        return self.filepath.name


@shapeclass
class Shape1(Shape):
    """Profile is to be represented as a deep VMR up to a certain ‘knee’ pressure, and
    then a defined fractional scale height.

    Parameters
    ----------
    knee_pressure : float
        Pressure at which the transition occurs.
    deep_vmr : float
        Deep volume mixing ratio.
    deep_vmr_error : float
        Error estimate for deep_vmr.
    fsh : float
        Fractional scale height.
    fsh_error : float
        Error estimate for fsh.
    """
    ID: ClassVar[int] = 1
    CONSTANTS: ClassVar[list[str]] = ["knee_pressure"]
    VARIABLES: ClassVar[list[str]] = ["deep_vmr", "fsh"]
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
    runname.ref, aerosol.ref, parah2.ref or fcloud.ref.

    Parameters
    ----------
    scale_factor : float
        Scaling factor applied to the reference profile.
    scale_factor_error : float
        Error estimate for scale_factor.
    """
    ID: ClassVar[int] = 2
    CONSTANTS: ClassVar[list[str]] = []
    VARIABLES: ClassVar[list[str]] = ["scale_factor"]
    scale_factor: float
    scale_factor_error: float

    def generate_apr_data(self):
        return f"{self.scale_factor} {self.scale_factor_error}"

@shapeclass
class Shape4(Shape):
    """Very similar to Shape1, but knee pressure is also a variable parameter.

    Parameters
    ----------
    knee_pressure : float
        Pressure at which the transition occurs.
    knee_pressure_error : float
        Error estimate for knee_pressure.
    deep_vmr : float
        Deep volume mixing ratio.
    deep_vmr_error : float
        Error estimate for deep_vmr.
    fsh : float
        Fractional scale height.
    fsh_error : float
        Error estimate for fsh.
    """
    ID: ClassVar[int] = 4
    CONSTANTS: ClassVar[list[str]] = []
    VARIABLES: ClassVar[list[str]] = ["deep_vmr", "fsh", "knee_pressure"]
    knee_pressure: float
    knee_pressure_error: float
    deep_vmr: float
    deep_vmr_error: float
    fsh: float
    fsh_error: float

    def generate_apr_data(self):
        return f"{self.knee_pressure} {self.knee_pressure_error}\n{self.deep_vmr} {self.deep_vmr_error}\n{self.fsh} {self.fsh_error}"

@shapeclass
class Shape20(Shape):
    """Similar to Shape1, but profile forced to a small number above tropopause.

    Parameters
    ----------
    knee_pressure : float
        Pressure at which the transition occurs.
    tropopause_pressure : float
        Pressure representing the tropopause.
    deep_vmr : float
        Deep volume mixing ratio.
    deep_vmr_error : float
        Error estimate for deep_vmr.
    fsh : float
        Fractional scale height.
    fsh_error : float
        Error estimate for fsh.
    """
    ID: ClassVar[int] = 20
    CONSTANTS: ClassVar[list[str]] = ["knee_pressure", "tropopause_pressure"]
    VARIABLES: ClassVar[list[str]] = ["deep_vmr", "fsh"]
    knee_pressure: float
    tropopause_pressure: float
    deep_vmr: float
    deep_vmr_error: float
    fsh: float
    fsh_error: float

    def generate_apr_data(self):
        return f"{self.knee_pressure} {self.tropopause_pressure}\n{self.deep_vmr} {self.deep_vmr_error}\n{self.fsh} {self.fsh_error}"

@shapeclass
class Shape32(Shape):
    """Cloud profile with variable base pressure, opacity, and fractional scale height.

    Parameters
    ----------
    base_pressure : float
        Cloud base pressure (bar).
    base_pressure_error : float
        Error estimate for base_pressure.
    opacity : float
        Cloud opacity (log-scaled).
    opacity_error : float
        Error estimate for opacity.
    fsh : float
        Fractional scale height (log-scaled).
    fsh_error : float
        Error estimate for fsh.
    """
    ID:  ClassVar[int] = 32
    CONSTANTS: ClassVar[list[str]] = []
    VARIABLES: ClassVar[list[str]] = ["opacity", "fsh", "base_pressure"]
    base_pressure: float
    base_pressure_error: float
    opacity: float
    opacity_error: float
    fsh: float
    fsh_error: float

    def generate_apr_data(self):
        return f"{self.base_pressure} {self.base_pressure_error}\n{self.opacity} {self.opacity_error}\n{self.fsh} {self.fsh_error}"
        
@shapeclass
class Shape37(Shape):
    """Cloud with constant opacity/bar between two pressure levels.

    Parameters
    ----------
    bottom_pressure : float
        Bottom pressure level (bar).
    top_pressure : float
        Top pressure level (bar).
    opacity : float
        Cloud opacity per bar.
    opacity_error : float
        Error estimate for opacity.
    """
    ID: ClassVar[int] = 37
    CONSTANTS: ClassVar[list[str]] = ["bottom_pressure", "top_pressure"]
    VARIABLES: ClassVar[list[str]] = ["opacity"]
    bottom_pressure: float
    top_pressure: float
    opacity: float
    opacity_error: float

    def generate_apr_data(self):
        return f"{self.bottom_pressure} {self.top_pressure}\n{self.opacity} {self.opacity_error}"

@shapeclass
class Shape47(Shape):
    """Cloud centred at specified pressure with variable FWHM and opacity.

    Parameters
    ----------
    central_pressure : float
        Pressure where cloud distribution peaks (bar).
    central_pressure_error : float
        Error estimate for central_pressure (bar).
    pressure_width : float
        Full width at half maximum in log pressure units.
    pressure_width_error : float
        Error estimate for pressure_width.
    opacity : float
        Total cloud opacity.
    opacity_error : float
        Error estimate for opacity.
    """
    ID: ClassVar[int] = 47
    CONSTANTS: ClassVar[list[str]] = []
    VARIABLES: ClassVar[list[str]] = ["opacity", "central_pressure", "pressure_width"]
    central_pressure: float
    central_pressure_error: float
    pressure_width: float
    pressure_width_error: float
    opacity: float
    opacity_error: float

    def generate_apr_data(self):
        return f"{self.opacity} {self.opacity_error}\n{self.central_pressure} {self.central_pressure_error}\n{self.pressure_width} {self.pressure_width_error}"

@shapeclass
class Shape48(Shape):
    """Cloud profile with variable base, top pressure, opacity, and scale height.

    Parameters
    ----------
    base_pressure : float
        Cloud base pressure (bar).
    base_pressure_error : float
        Error estimate for base_pressure (bar).
    top_pressure : float
        Cloud top pressure (bar).
    top_pressure_error : float
        Error estimate for top_pressure (bar).
    opacity : float
        Cloud opacity.
    opacity_error : float
        Error estimate for opacity.
    fsh : float
        Fractional scale height.
    fsh_error : float
        Error estimate for fsh.
    """
    ID:  ClassVar[int] = 48
    CONSTANTS: ClassVar[list[str]] = []
    VARIABLES: ClassVar[list[str]] = ["opacity", "fsh", "base_pressure", "top_pressure"]
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


def get_shape_from_id(id_):
    """Given a shape ID integer, return a reference to the class corresponding to that ID. Note that this returns a class, not an instantiated object.
    
    Args:
        id_: Shape ID to get
        
    Returns:
        Shape: Class object of the corresponding shape"""
    id_ = int(id_)
    for shape_class in ALL_SHAPES:
        if shape_class.ID == id_:
            return shape_class


ALL_SHAPES = [Shape0, Shape1, Shape2, Shape4, Shape20, Shape32, Shape37, Shape47, Shape48] #: A list of all the Shape classes