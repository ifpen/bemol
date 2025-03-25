"""Definition of simple rotor classes.

Instances of the rotors available in the rotors folder are generated
considering the default properties.

To use it: 'bemol.rotor.nameOfTheRotor'

"""

from types import SimpleNamespace
import yaml
import os
from pathlib import Path

import numpy as np

from . import airfoil
from . import section


class Rotor(object):
    """Class to define a rotor.

    Blades are considered identical and constant.
    
    Parameters
    ----------
    folder :
        path to rotor folder. Must contain the `blade.dat` with the blade
        defintion, a `rotor.yml` file with general rotor properties and the
        airfoils polars (`*.foil`) 
    airfoil : 
        class for defining the airfoil. airfoil.BaseAirfoil by default
    
    """

    def __init__(self,folder,airfoil=airfoil.BaseAirfoil):

        blade_file = Path(folder)/'blade.dat'
        if not blade_file.is_file():
            raise ValueError(f'Blade definition file not found in model folder ({folder}).')
        
        blade_data = np.genfromtxt(blade_file,skip_header=1,dtype=str)

        self.radius = np.asarray(blade_data[:,0],dtype=float)
        self.twist = np.asarray(blade_data[:,1],dtype=float)
        self.chords = np.asarray(blade_data[:,2],dtype=float)

        list_sections = blade_data[:,3]
        self.airfoils = SimpleNamespace()
        self.sections = []
        # get the airfoils
        for i in range(len(list_sections)):
            airfoil_name = list_sections[i]
            airfoil_file = Path(folder)/'airfoils'/f'{airfoil_name}.foil'
            if not hasattr(self.airfoils,airfoil_name):
                # no duplicates
                setattr(self.airfoils,airfoil_name,airfoil(airfoil_file))
            # defining section objects to store the information of each section
            self.sections.append(
                section.Section(
                    radius=self.radius[i],twist=self.twist[i],chord=self.chords[i],
                    airfoil=getattr(self.airfoils,airfoil_name)
                    )
                )
        rotor_file = Path(folder)/'rotor.yml'
        if not blade_file.is_file():
            raise ValueError(f'Rotor definition file not found in model folder ({folder}).')
        with open(rotor_file) as stream:
            try:
                data_rotor = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        for name, value in data_rotor.items():
            setattr(self,name,value)
    

    def __getitem__(self,item):
        """Return section."""
        return self.sections[item]
    

    def __iter__(self):
        """Return iterator to list of sections."""
        return self.sections.__iter__()


# creating instances of predefined rotors with default properties
script_path = os.path.abspath(__file__)
for available_rotor in Path(script_path).parent.glob('rotors/*/'):
    rotor_name = available_rotor.name
    # TODO: remove use of globals, not good practice!
    globals()[rotor_name] = Rotor(available_rotor)
del script_path, rotor_name, available_rotor