


import numpy as np
from scipy import optimize
from scipy.optimize import minimize

from . import bem
from . import section
from . import tools

PI = float(np.pi)
PI_HALF = float(np.pi/2.0)
PI_QUARTER = float(np.pi/4.0)
TWO_PI = float(2.0*np.pi)


class NingUncoupled(bem.BaseBEM):
    """Ning uncoupled.
    
    Solves the BEM equations for a given section.
    
    Parameters
    ----------
    args
        parent class (baseBEM) arguments.
    kwargs
        parent class (baseBEM) key arguments.
    epsilon : float, optional
        scalar for setting up intervals of BEM algorithm solution.
    
    Attributes
    ----------
    epsilon: float
        scalar for setting up intervals of BEM algorithm solution.
    _axial_induction: float
        axial induction factor.
    _tangential_induction: float
        tangential induction factor.
    _Ux: float
        axial (along airfoil) velocity, m/s
    _Uy: float
        tangential (across airfoil) velocity, m/s
    _angle: float
        sectional angle (twsit + pitch), radians
    _funLift: callable
        function that return the lift coefficient from the angle of attack.
    _funDrag: callable
        function that return the drag coefficient from the angle of attack.
    _chord: float
        sectional chord, m.
    _radius: float
        sectional radius, m.

    """

    def __init__(self,*args,epsilon=1e-5,**kwargs):
        super().__init__(*args,**kwargs)

        self.epsilon = epsilon

        self._Ux = 0.0
        self._Uy = 0.0
        self._axial_induction = 0.0
        self._tangential_induction = 0.0

        self._angle = 0.0
        self._funLift = None
        self._funDrag = None
        self._chord = 0.0
        self._radius = 0.0


    def residuals(self,inflowAngle:float,) -> float:
        """Residuals of uncoupled Ning algorithm.
        
        Solution of the base BEM equations. Only the hub/tip loss correction is
        considered, if provided.

        Parameters
        ----------
        inflowAngle: float
            inflow angle in radians
        
        Return
        ------
        residual value
        
        """
        
        attackAngle = inflowAngle - self._angle

        lift = self._funLift(attackAngle)
        drag = self._funDrag(attackAngle)

        Cx = lift*np.cos(inflowAngle) + drag*np.sin(inflowAngle)
        Cy = lift*np.sin(inflowAngle) - drag*np.cos(inflowAngle)

        F = self.corrections.hubTipLoss(
            self._radius,self.rotor.nBlades,
            self.rotor.hubRadius,self.rotor.tipRadius,inflowAngle
            )

        sigmaPrime = self.rotor.nBlades * self._chord / ( TWO_PI*self._radius )
        kappa = sigmaPrime * Cx / (4.*F*np.sin(inflowAngle)**2.)
        kappaPrime = sigmaPrime * Cy / (4.*F*np.sin(inflowAngle)*np.cos(inflowAngle))

        axialInd = 0.0
        tangentialInd = 0.0
        maxAxialInd = 2.0
        maxTangentialInd = 2.0

        isValidSolution = True
        isMomentumRegion = False
        if ( (inflowAngle >= 0.0 and self._Ux >= 0.0) or (inflowAngle < 0.0 and self._Ux < 0.0)):
            isMomentumRegion = True

        if isMomentumRegion == True:
            if kappa <= 2.0/3.0:
                axialInd = kappa / (1.+kappa)
                if kappa < -1.0:
                    isValidSolution = False
            else:
                twoF_Kappa = 2.0*F*kappa
                gamma1 = twoF_Kappa - (10.0/9.0 - F)
                gamma2 = twoF_Kappa - F*(4.0/3.0 - F)
                gamma3 = twoF_Kappa - (25.0/9.0 - 2.0*F)

                if np.abs(gamma3) < self.epsilon:
                    axialInd = 1.0 - 1.0/(2.0*np.sqrt(gamma2))
                else:
                    axialInd = (gamma1 - np.sqrt(np.abs(gamma2))) / gamma3
        else:
            # propeller brake region
            axialInd = kappa / (kappa - 1.0)
            if kappa <= 1.0:
                isValidSolution = False
            elif (axialInd > maxAxialInd):
                axialInd = maxAxialInd

        if self._Ux <= 0.0:
            kappaPrime = -kappaPrime

        tangentialInd = kappaPrime / (1. - kappaPrime)
        if np.abs(tangentialInd) > maxTangentialInd:
            tangentialInd = maxTangentialInd*np.sign(tangentialInd)

        if isValidSolution == False:
            axialInd = 0.0
            tangentialInd = 0.0

        if isMomentumRegion:
            residuals = np.sin(inflowAngle)/(1.0 - axialInd) \
                        - self._Ux/self._Uy*np.cos(inflowAngle)/(1.0 + tangentialInd)
        else:
            residuals = np.sin(inflowAngle)/(1.0 - kappa) \
                        - self._Ux/self._Uy*np.cos(inflowAngle)/(1.0 + tangentialInd)

        self._axial_induction = axialInd
        self._tangential_induction = tangentialInd

        return residuals


    def solve(self,section:section.Section,azimuth:float,pitch:float,
              velocity:list=[0.0,0.0,0.0],angles=[0.0,0.0],tStep=0.0):
        """Solve the uncoupled Ning algorithm.
        
        Parameters
        ----------
        section : bemol.section.Section
            section object storing geometrical information (chord, radius, twist
            and airfoil polar).
        azimuth : float
            azimuthal angle of blade, radians.
        pitch : float
            blade pitch angle, in radians.
        velocity : list, optional
            axial, tangential (normally omega * r) and radial velocity, m/s.
        angles : list, optional
            inflow yaw and tilt angles, radians.
        tStep: float, optional
            timestep duration, seconds.
        
        Returns
        -------
        sectional forces and induction factors.
        
        """

        # tilt not used!
        yaw = angles[0]

        self._axial_induction = 0.0
        self._tangential_induction = 0.0

        angle = section.twist + pitch
        chord = section.chord
        radius = section.radius
        funDrag = section.airfoil.cd
        funLift = section.airfoil.cl

        # radial velocity not used
        Ux = velocity[0]
        Uy = velocity[1]

        # update the flow state before calculating the residuals
        self.update(
            Ux=Ux,Uy=Uy,
            angle=angle,funLift=funLift,funDrag=funDrag,
            chord=chord,radius=radius,
            )
        
        residualEpsilon = self.residuals(self.epsilon)
        residualPiOvTwo = self.residuals(PI_HALF)

        if residualEpsilon * residualPiOvTwo < 0.0:
            inflowAngle = optimize.brentq(self.residuals,self.epsilon,PI_HALF,)
        else:
            residualMinusEpsilon = self.residuals(-self.epsilon)
            residualMinPiOvFour = self.residuals(-PI_QUARTER)

            if residualMinusEpsilon*residualMinPiOvFour < 0.0:
                # propeller break region
                inflowAngle = optimize.brentq(
                    self.residuals,-PI_QUARTER,-self.epsilon,
                    )
            else:
                inflowAngle = optimize.brentq(
                    self.residuals,PI_HALF,PI - self.epsilon,
                    )

        # yawModel: apply to axial induction only
        wakeSkewAngle = self.corrections.skewAngle(
            self._axial_induction,yaw)
        self._axial_induction = self.corrections.yawModel(
            self._axial_induction,wakeSkewAngle,azimuth,radius,
            self.rotor.hubRadius,self.rotor.tipRadius
            )
        self._axial_induction = self.corrections.dynamicInflow(
            self._axial_induction,Ux,radius,tStep)

        uxRelative = Ux * (1.0 - self._axial_induction)
        uthetaRelative = Uy * (1.0 + self._tangential_induction)
        inflowAngle = np.arctan2(uxRelative, uthetaRelative)

        attackAngle = inflowAngle - angle
        liftCoeff = funLift(attackAngle)
        dragCoeff = funDrag(attackAngle)

        normalCoeff = liftCoeff * np.cos(attackAngle) + dragCoeff * np.sin(attackAngle)
        tangentialCoeff = -liftCoeff * np.sin(attackAngle) + dragCoeff * np.cos(attackAngle)

        uRelative = np.sqrt(uxRelative**2. + uthetaRelative**2.)

        normalForce = 0.5*self.rho*uRelative**2.*chord*normalCoeff
        tangentialForce = 0.5*self.rho*uRelative**2.*chord*tangentialCoeff

        return normalForce, tangentialForce, self._axial_induction, self._tangential_induction




class NingCoupled(bem.BaseBEM):
    """Coupled BEM algorithm from Ning.
    
    Corrections of the uncoupled algorithm are defined as the same as in the
    coupled algorithm. To check if this is the desired behavior in case of
    dynamic inflow correction.

    Attributes
    ----------
    _uncoupled  : bemol.ning.NingUncoupled
        uncoupled solver instance.
    _uInfty : float
        inflow (wind) velocity, m/s.
    _UxPrime : float
        complete (shear, turbulence, structural vibration) axial velocity, m/s
    _UyPrime : float
        complete (shear, turbulence, structural vibration) tangential velocity, m/s
    _skew : float
        wake skew angle, radians.
    _yaw : float
        inflow yaw angle, radians.
    _azimuth: float
        azimuthal angle, radians.

    See `NingUncoupled` for the remaining attributes.
    """

    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)

        # replace the last argument (corrections) by the instance of the
        # main solver
        if len(args) > 2:
            args = list(args)
            args[2] = list(self.corrections)
        if 'corrections' in kwargs:
            kwargs['corrections'] = self.corrections
        
        self._uncoupled = NingUncoupled(*args,**kwargs)

        self._uInfty = 0.0
        self._UxPrime = 0.0
        self._UyPrime = 0.0
        self._skew = 0.0 
        self._yaw = 0.0
        self._azimuth = 0.0
        
        self._axial_induction = 0.0
        self._tangential_induction = 0.0

        self._angle = 0.0
        self._funLift = None
        self._funDrag = None
        self._chord = 0.0
        self._radius = 0.0


    def pre(self,section,omega,wind,precone) -> dict:
        """Calculate prime velocities at a given radius."""
        UxPrime = wind*np.cos(precone)
        UyPrime = omega*section.radius*np.cos(precone)
        return dict(UxPrime=UxPrime,UyPrime=UyPrime,)


    def residuals(self,inductions,) -> float:
        """Residuals of coupled Ning algorithm.
        
        Parameters
        ----------
        inductions : list
            list of axial and tangential inductions.

        Returns
        -------
        square sum of the residuals.

        """

        residuals = np.zeros(2)

        axialInduction = inductions[0]
        tangentialInduction = inductions[1]

        wakeSkewAngle = self.corrections.skewAngle(axialInduction,self._skew)

        # To get the inflow angle for tip loss evaluation
        uxHub = self._UxPrime*(np.cos(self._skew)-axialInduction) \
                + self._UyPrime*tangentialInduction*np.sin(wakeSkewAngle)*np.cos(self._azimuth)*(1.+np.sin(wakeSkewAngle)*np.sin(self._azimuth))
        uyHub = self._UyPrime*(1.+tangentialInduction*np.cos(wakeSkewAngle)*(1.+np.sin(wakeSkewAngle)*np.sin(self._azimuth))) \
                + self._UxPrime*np.cos(self._azimuth)*(axialInduction*np.tan(wakeSkewAngle/2.)-np.sin(self._skew))

        W2 = uxHub**2. + uyHub**2.
        sigmaPrime = self.rotor.nBlades * self._chord / (TWO_PI*self._radius)

        inflowAngle = np.arctan2(uxHub, uyHub)
        attackAngle = inflowAngle - self._angle
        lift = self._funLift(attackAngle)
        drag = self._funDrag(attackAngle)

        Cx = lift * np.cos(inflowAngle) + drag * np.sin(inflowAngle)
        Cy = lift * np.sin(inflowAngle) - drag * np.cos(inflowAngle)

        # tip loss evaluation
        lossFactor = self.corrections.hubTipLoss(
            self._radius,self.rotor.nBlades,
            self.rotor.hubRadius,self.rotor.tipRadius,inflowAngle
            )

        CtElement = W2 / self._uInfty ** 2. * sigmaPrime * Cx
        CqElement = W2 / self._uInfty ** 2. * sigmaPrime * (
            Cy*np.cos(wakeSkewAngle) - Cx*np.sin(wakeSkewAngle)*np.cos(self._azimuth)
            )

        residuals[0] = CtElement - self.CtMomentum(
            axialInduction,tangentialInduction,lossFactor,self._skew
            )
        residuals[1] = CqElement - self.CqMomentum(
            axialInduction,tangentialInduction,self._UxPrime,self._UyPrime,
            lossFactor,wakeSkewAngle,self._azimuth,self._yaw
            )

        return residuals[0]*residuals[0] + residuals[1]*residuals[1]
    

    def solve(self,section:section.Section,azimuth:float,pitch:float,
              velocity:list,angles:list=[0.0,0.0],
              uInfty=None,UxPrime=None,UyPrime=None,
              skew=0.0,tStep=0.0):
        """Solve the coupled Ning algorithm.

        If velocities are not explicitely given, consider the inflow velocity
        as the first component of the `velocity` argument.
        
        Parameters
        ----------
        uInfty: float, optional
            inflow (wind) velocity, m/s.
        UxPrime: float, optional
            complete (shear, turbulence, structural vibration) axial velocity, m/s,
        UyPrime: float, optional
            complete (shear, turbulence, structural vibration) tangential velocity, m/s,
        skew: float, optional
            wake skew angle, radians.

        See `NingUncoupled.solve` for other arguments.

        Returns
        -------
        sectional forces and induction factors.
        
        """

        uInfty = velocity[0] if uInfty is None else uInfty
        UxUncoupled = velocity[0]
        UyUncoupled = velocity[1]
        UxPrime = velocity[0] if UxPrime is None else UxPrime
        UyPrime = velocity[1] if UyPrime is None else UyPrime
        
        yaw = angles[0]
        
        ## Initialize minimization algorithm and solve
        # initialize inductions using Uncoupled algorithm
        self._uncoupled.update(
            Ux=UxUncoupled,Uy=UyUncoupled,
            angle=section.twist + pitch,
            funLift=section.airfoil.cl,funDrag=section.airfoil.cd,
            chord=section.chord,radius=section.radius,
            )
        
        _, _, axInd, tanInd = self._uncoupled.solve(
            section,azimuth,pitch,
            velocity=[UxUncoupled,UyUncoupled],angles=angles,
            tStep=tStep)

        inductionsInit = np.array((axInd,tanInd))

        angle = section.twist + pitch
        funDrag = section.airfoil.cd
        funLift = section.airfoil.cl

        self.update(
            uInfty=uInfty,UxPrime=UxPrime,UyPrime=UyPrime,
            skew=skew,yaw=yaw,azimuth=azimuth,
            chord=section.chord,radius=section.radius,angle=angle,
            funLift=funLift,funDrag=funDrag
            )

        # solve for axial and tangential inductions
        res = minimize(
            self.residuals,inductionsInit,
            method='Powell',bounds=((axInd-0.2, axInd+0.2),(tanInd-0.2, tanInd+0.2)),
            tol=1e-6,options={'disp': False}
            )

        ## Postprocessing
        inductions = res.x
        axialInduction = inductions[0]
        tangentialInduction = inductions[1]

        wakeSkewAngle = self.corrections.skewAngle(axialInduction,skew)
        axialInduction = self.corrections.dynamicInflow(
            axialInduction,UxUncoupled,section.radius,tStep
            )

        # To get the inflow angle for tip loss evaluation
        uxHub = UxPrime*(np.cos(skew)-axialInduction) \
                + UyPrime*tangentialInduction*np.sin(wakeSkewAngle)*np.cos(azimuth)*(1.+np.sin(wakeSkewAngle)*np.sin(azimuth))
        uyHub = UyPrime*(1.+tangentialInduction*np.cos(wakeSkewAngle)*(1.+np.sin(wakeSkewAngle)*np.sin(azimuth))) \
                + UxPrime*np.cos(azimuth)*(axialInduction*np.tan(wakeSkewAngle/2.)-np.sin(skew))

        inflowAngle = np.arctan2(uxHub, uyHub)

        attackAngle = inflowAngle - angle
        liftCoeff = funLift(attackAngle)
        dragCoeff = funDrag(attackAngle)

        normalCoeff = liftCoeff * np.cos(attackAngle) + dragCoeff * np.sin(attackAngle)
        tangentialCoeff = -liftCoeff * np.sin(attackAngle) + dragCoeff * np.cos(attackAngle)

        uRelative = np.sqrt(uxHub**2.+uyHub**2.)

        normalForce = 0.5*self.rho*uRelative**2.*section.chord*normalCoeff
        tangentialForce = 0.5*self.rho*uRelative**2.*section.chord*tangentialCoeff

        return normalForce, tangentialForce, axialInduction, tangentialInduction








