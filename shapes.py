from dataclasses import dataclass
from typing import ClassVar
import shutil


""" HOW TO ADD A SHAPE: 
All subclasses of Shape must implement:
    - ID: Integer ID of the shape in NEMESIS
    - NAMES: The names of the parameters in order they appear in the mre file (do not include error parameters here)
    - __init__: Instantiates the parent class and defines:
        - nemesis_code: The number of the profile in the NEMESIS manual
        - Any parameters associated with the model (eg. knee pressure, fsh etc). These MUST be defined in the order that they appear in the .apr file
    - generate_apr_data: Returns a string of the profile shape parameters in the format
                         that the .apr file expects.

Additionally, some shapes require additional files in the core (eg. tempapr.dat). If this is the case,
override the copy_required_files method.

Finally, add the new ShapeN class to the ALL_SHAPES list at the end of the file.

Class template:

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

"""


def shapeclass(*args, **kwargs):
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

    def copy_required_files(self, directory):
        """Some Shapes require additional files to be copied into the core directory"""
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

    def copy_required_files(self, directory):
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
    NAMES: ClassVar[list[str]] = ["base_pressure", "fsh", "opacity"]
    base_pressure: float
    base_pressure_error: float
    opacity: float
    opacity_error: float
    fsh: float
    fsh_error: float

    def generate_apr_data(self):
        return f"{self.base_pressure} {self.base_pressure_error}\n{self.opacity} {self.opacity_error}\n{self.fsh} {self.fsh_error}"
        

def get_shape_from_id(id_):
    id_ = int(id_)
    for shape_class in ALL_SHAPES:
        if shape_class.ID == id_:
            return shape_class




ALL_SHAPES = [Shape0, Shape1, Shape32]


if __name__ == "__main__":
    x = Shape1(1,2,3,4,5)
    print(x)