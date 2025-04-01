"""

Secondary correction models

"""

import numpy as np
from . import bem


INV_TWO_PI = float(2.0/np.pi)


class HubTipLoss:

    class Dummy:
        """Empty hub/tip loss correction model."""

        def __init__(self) -> None:
            return

        def __call__(self,*args,**kwargs):
            return 1.0
    
    class Prandtl:
        """Prandt hub/tip loss correction model.
        
        Parameters
        ----------
        epsilon: float, optional
            minimal value of F to avoid errors close to/at the bounds of the
            blade. Default is 1e-12.
        
        """
        def __init__(self,epsilon=1e-12) -> None:
            self._epsilon = epsilon
            return

        def __call__(self,radius,nBlades,hubRadius,tipRadius,inflowAngle):

            fTip = nBlades/2.0*((tipRadius - radius) / (radius*np.abs(np.sin(inflowAngle))))
            fTip = INV_TWO_PI*np.arccos(np.exp(-fTip))

            fRoot = nBlades/2.0*((radius - hubRadius) / (hubRadius*np.abs(np.sin(inflowAngle))))
            fRoot = INV_TWO_PI*np.arccos(np.exp(-fRoot))

            return max(fTip*fRoot,self._epsilon)


class SkewAngle:

    class Dummy:
        """Empty skewed wake model."""

        def __init__(self) -> None:
            return

        def __call__(self,axialInduction,*args,**kwargs):
            return axialInduction
    
    class Burton:
        """Burton skewed wake model."""

        def __init__(self,) -> None:
            return

        def __call__(self,axialInduction,yawAngle):
            return (0.6 * axialInduction + 1.0) * yawAngle


class DynamicInflow:

    class Dummy:
        """Empty dynamic inflow model."""

        def __init__(self) -> None:
            return

        def __call__(self,axialInduction,*args,**kwargs):
            return axialInduction
    
    class Knudsen:
        """Knudsen dynamic inflow model."""

        def __init__(self,alphaDynamic=0.3,tauScale=3.0) -> None:
            self._alpha0 = alphaDynamic
            self.alphaDynamic = alphaDynamic
            self.tauScale = tauScale

        def __call__(self,axialInduction,Ux,radius,tStep):
            if (tStep == 0.):
                raise ValueError(
                    f'Using a dynamic inflow model with tStep = {tStep} does not make sense!'
                    )
            # lower bound for wind velocity
            kappa = 1. / (self.tauScale * radius / max(1.0,Ux) )
            alphaDynamic = tStep*kappa*(axialInduction - self.alphaDynamic) + self.alphaDynamic
            self.alphaDynamic = alphaDynamic

            return alphaDynamic
        
        def restart(self):
            """Restart alpha to the starting value."""
            self.alphaDynamic = self._alpha0


class YawModel:

    class Dummy:
        """Empty skewed wake correction model."""
        
        def __init__(self) -> None:
            return

        def __call__(self,axialInduction,*args,**kwargs):
            return axialInduction
        
    class PittAndPeters:
        """Pitt & Peters skewed wake model."""

        def __init__(self,factor=15.*np.pi/64.) -> None:
            self.factor = factor
        
        def __call__(self,axialInduction,wakeSkewAngle,azimuthAngle,radius,tipRadius):
            return axialInduction*(
                1. + self.factor*np.tan(wakeSkewAngle/2.)*radius/tipRadius*np.sin(azimuthAngle)
                )
        

class TurbulentWakeState:

    class Dummy:
        """Empty high-induction correction model.
        
        Calls the base thrust coefficient calculation method.
        """
        
        def __init__(self) -> None:
            return

        def __call__(self,axialInduction,lossFactor,skewAngle,*args,**kwargs):
            return bem.BaseBEM.CT(axialInduction,lossFactor,skewAngle)


    class Buhl:
        """Buhl's empirical high-induction correction model.
        
        TODO: need to find the ref to this method.
        
        """

        def __init__(self,beta:float=0.4,da:float=0.02) -> None:
            self.beta = beta
            self.da = da
        
        def __call__(self,axialInduction,lossFactor,skewAngle):
            f0 = bem.BaseBEM.CT(self.beta,lossFactor,skewAngle)
            fp0 = (
                bem.BaseBEM.CT(self.beta + self.da,lossFactor,skewAngle)
                - bem.BaseBEM.CT(self.beta - self.da,lossFactor,skewAngle)
                ) / (2.*self.da)
            f1 = 2.*np.cos(skewAngle)
            k2 = (f1-f0-fp0*(1.0 - self.beta)) / (1.0 - self.beta)**2.
            k1 = fp0 - 2.*k2*self.beta
            k0 = f1 - k1 - k2

            return k0 + k1*axialInduction + k2*axialInduction**2.