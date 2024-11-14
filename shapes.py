from dataclasses import dataclass
from typing import ClassVar


""" HOW TO ADD A SHAPE: 
All subclasses of Shape must implement:
    - ID: Integer ID of the shape in NEMESIS
    - __init__: Instantiates the parent class and defines:
        - nemesis_code: The number of the profile in the NEMESIS manual
        - Any parameters associated with the model (eg. knee pressure, fsh etc)
    - generate_apr_data: Returns a string of the profile shape parameters in the format
                         that the .apr file expects.

eg:

@dataclass
class ShapeN(Shape):
    ID:  ClassVar[int] = N
    arg_name_1: type_1
    arg_name_2: type_2
    ...

    def generate_apr_data(self):
        return a string for apr file
        
Finally, add the new ShapeN class to the ALL_SHAPES list at the end of the file

"""


class Shape:
    def __init__(self):
        """Do not instantiate directly, use a subclass"""
        pass


@dataclass
class Shape0(Shape):
    ID: ClassVar[int] = 0
    filepath: str

    def generate_apr_data(self):
        return self.filepath


@dataclass
class Shape1(Shape):
    ID: ClassVar[int] = 1
    knee_pressure: float
    deep_vmr: float
    deep_vmr_error: float
    fsh: float
    fsh_error: float

    def generate_apr_data(self):
        return f"{self.knee_pressure}\n{self.deep_vmr} {self.deep_vmr_error}\n{self.fsh} {self.fsh_error}"


@dataclass
class Shape32(Shape):
    ID:  ClassVar[int] = 32
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
