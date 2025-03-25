"""

Classes defining airfoils and methods for working with polars.

"""

from pathlib import Path

import numpy as np
from scipy import interpolate
from scipy import optimize


class BaseAirfoil(object):
    """Base class to define an airfoil.
    
    Parameters
    ----------
    file
        path to .foil file with the polars definition.
    interpolator : callable, optional
        function for defining the interpolation method to look for the
        aerodynamic coefficients given an angle of attack. Default is 
        scipy interpolate.interp1d with default parameters.
    
    """
    
    def __init__(self,file,interpolator=interpolate.interp1d):
        airfoil_file = Path(file)
        if not airfoil_file.is_file():
            raise ValueError(f'Airfoil (polar) file not found ({file}).')
        self.file = file

        polar_data = np.genfromtxt(airfoil_file)
        self._alpha = np.radians(polar_data[:,0])
        self._cl = polar_data[:,1]
        self._cd = polar_data[:,2]

        self._interp_cl = interpolator(self._alpha,self._cl)
        self._interp_cd = interpolator(self._alpha,self._cd)
    

    @property
    def polar(self):
        """Return alpha, lift and drag coefficients."""
        return self._alpha, self._cl, self._cd
    

    @property
    def alpha(self):
        """Angles of attack of polar in radians."""
        return self._alpha
    

    def cl(self,aoa):
        """Lift coefficients of polar at angle of attack."""
        return self._interp_cl(aoa)
    

    def cd(self,aoa):
        """Drag coefficients of polar at angle of attack."""
        return self._interp_cd(aoa)
    


class DynStallAirfoil(BaseAirfoil):
    """Airfoil class with dynamic stall correction.
    
    Under construction.
    
    """

    def __init__(self,file):
        super().__init__(file)

        self.splineLift = interpolate.interp1d(self._alpha,self._cl)
        self.splineDrag = interpolate.interp1d(self._alpha,self._cd)

        self.splineStaticAttachment = 0.
        self.splineFullySeparatedLift = 0.

        self.cl0Slope = 0.
        self.aoaZero = 0.

        self.estimateAoaZero()
        self.estimateLiftSlope()
        self.estimateAttachmentDegree(self.alpha)
        self.estimateFullySeparatedLift(self.alpha)


    def estimateFullySeparatedLift(self, attackAngles):

        clFS = []
        for attackAngle in attackAngles:

            coeff = self.cl(attackAngle)
            clFA = self.cl0Slope * (attackAngle - self.aoaZero)
            staticAttachment = self.splineStaticAttachment(attackAngle)

            limit = 0.975
            if (staticAttachment < limit):
                cl_s = (coeff - staticAttachment * clFA)/(1.-min(1.,staticAttachment))
            else:
                # Numerical issues...
                cl_s = coeff / 2.
            clFS.append(cl_s)

        self.splineFullySeparatedLift = interpolate.interp1d(attackAngles,clFS)
        return


    def estimateAoaZero(self):
        self.aoaZero = optimize.brentq(self.splineLift,np.radians(-10.),0.,args=())
        return


    def estimateLiftSlope(self):
        dA = np.radians(3.)
        self.cl0Slope = (self.splineLift(self.aoaZero+dA) - self.splineLift(self.aoaZero)) / dA
        return


    def estimateAttachmentDegree(self, attackAngles):
        staticAttachment = []
        for attackAngle in attackAngles:
            coeff = self.cl(attackAngle)
            clFA = self.cl0Slope*(attackAngle-self.aoaZero)
            # Since clFA is an estimate, sometimes for aoa<0 lift < clFA
            static_attachment = (np.sqrt(np.abs(coeff/clFA))*2.0 - 1.0)**2.0
            staticAttachment.append( max(0.0,min(static_attachment,1.0)) )

        self.splineStaticAttachment = interpolate.interp1d(attackAngles,staticAttachment)

        return

