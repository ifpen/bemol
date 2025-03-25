"""

Section class.

"""

class Section(object):
    """Section class.
    
    Object to define a section. A section is a node of the blade defined
    by its radius, chord, twist and airfoil.
    
    Parameters
    ----------
    radius : float
        section radius, m
    twist : float
        section twist angle, rad
    chord : float
        section chord, m
    airfoil :
        section airfoil. Object that contains functions for returning the polar
        coefficients


    TODO: inplement getter and setters
    """

    def __init__(self,radius:float,twist:float,chord:float,airfoil):
        self.radius = radius
        self.twist = twist
        self.chord = chord
        self.airfoil = airfoil