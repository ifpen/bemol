

import inspect

import numpy as np

from . import rotor
from . import secondary
from . import tools




class Corrections(object):
    """Class for storing corrections.

    Consider all the corrections available in the secondary
    models module. If not given by the user considers the base (empty)
    correction with the default parameters.

    Parameters
    ----------
    corrections : optional
        dictionary or list with the secondary corrections, either classes of
        instances - in case of custom parameters. If dictionary is
        given, the key must be the name of the correction as defined in
        secondary.py with a lower first letter.
    
    """
    def __init__(self,corrections:dict={}):
        for name, obj in inspect.getmembers(secondary):
            if inspect.isclass(obj):
                _name = name[0].lower() + name[1:]
                if type(corrections) is dict:
                    if _name in corrections:
                        # instantiation of correction with default values
                        corr = corrections[_name]
                        corr = corr() if isinstance(corr,type) else corr
                        setattr(self,_name,corr)
                    else:
                        setattr(self,_name,obj.Dummy())
                else:
                    effect_name = obj.__qualname__
                    setattr(self,_name,obj.Dummy())
                    # loop for all effects to check if any of the input
                    # corrections are inner of the available corrections
                    for corr in corrections:
                        # instantiation of correction with default values
                        corr = corr() if isinstance(corr,type) else corr
                        correction_name = type(corr).__qualname__
                        if effect_name in correction_name:
                            setattr(self,_name,corr)
                            break


    def __iter__(self):
        """Iterate corrections."""
        for value in self.__dict__.values():
            yield value


class BaseBEM:
    """Base BEM class.

    Parameters
    ----------
    rotor : bemol.rotor.Rotor
        instance of a rotor.
    rho : float
        air density in kg/m3, optional
    corrections : optional
        dictionary or list with the secondary corrections, either classes of
        instances - in case of custom parameters. If dictionary is
        given, the key must be the name of the correction as defined in
        secondary.py with a lower first letter.
        
        Cannot given a custom correction model for the moment, they must be
        coded as in secondary.py

    TODO: filter the corrections so 2 of the same type cannot be used, raise
          error or pick first and raise warning.

    """
    
    def __init__(self,rotor:rotor.Rotor,rho:float=1.225,corrections:dict=None):
        
        self.rotor = rotor
        self.rho = rho
        self.n = len(rotor.sections)

        if corrections is None: corrections = {}
        self.corrections = Corrections(corrections)
        


    @staticmethod
    def CT(a,F:float=1.0,chi:float=0.0) -> float:
        """Base thrust coefficient formulation.
        
        When there is no loss and skew is null: Ct = 4a(a-1)

        Parameters
        ----------
        a : float
            axial induction factor
        F : float, optional
            loss factor
        chi : float, optional
            wake skew angle in radians

        Returns
        -------
        thrust coefficient.

        """
        return 4.0*a*F*np.sqrt(1.0 - a*(2.*np.cos(chi) - a))


    def CtMomentum(self,a,aprime,F,chi) -> float:
        """Calculate thrust coefficient.

        The aprime (tangential induction) is not used.

        Parameters
        ----------
        a : float
            axial induction factor
        aprime : float
            tangential induction factor
        F : float
            loss factor
        chi : float
            wake skew angle in radians
        
        Returns
        -------
        thrust coefficient.

        """
        beta = 0.4

        if a <= beta:
            # Momentum region
            return self.CT(a,F,chi)
        elif a > 1.0:
            # Propeller brake region
            print('Unexpected propeller brake state...')
            return 4.*a*F*(a - 1.0)*np.cos(chi)
        else:
            # Empirical region
            return self.corrections.turbulentWakeState(a,F,chi)


    def CqMomentum(self,a:float,aprime:float,Ux:float,Uy:float,
                   F:float,chi:float,psi:float,gamma:float) -> float:
        """Calculate torque coefficient.
        
        Parameters
        ----------
        a : float
            axial induction
        aprime : float
            tangential induction
        Ux : float
            axial velocity, m/s
        Uy : float
            tangential velocity, m/s
        F : float, optional
            loss factor
        y : float, optional
            skew angle in radians
        chi : float
            wake skew angle in radians
        psi : float
            azimuth angle in radians
        gamma : float
            yaw angle in radians

        Returns
        -------
        torque coefficient.

        """
        # yaw instead of skew
        return 4.0*(Uy/Ux)*aprime*F*(np.cos(gamma) - a)*(
            np.cos(psi)**2. + np.cos(chi)**2.*np.sin(psi)**2.
            )
    

    def update(self,**kwargs):
        """Update the variables indicating the state of the section.
        
        Considers that the attributes are defined with a leading underscore
        (protected). If the attribute is not found raises an error.

        When defining new model always define the attributes of interest
        in the constructor.

        Parameters
        ----------
        kwargs
            key arguments with the parameters to update.
         
        """
        for key, val in kwargs.items():
            if hasattr(self,f'_{key}',):
                setattr(self,f'_{key}',val)
            else:   
                raise ValueError(f'This solver has no attribute {key}')
    

    def pre(self,*args,**kwargs) -> dict:
        """Function to update inputs of solver for given azimuth and section.

        Function that calculates/updates parameters needed for the call of
        the steady solution. It is used to update values that evolve with
        the section (radius, chord, airfoil), azimuthal angle, angular
        velocity or any other parameter.
        
        To be replaced for each model that needs it.

        """
        return {}


    def steady(self,azimuth:float,pitch:float,wind:float,omega:float,*args,
               angles:list=[0.0,0.0],precone:float=0.0,elements=slice(None),
               **kwargs):
        """Solving the BEM equations for a given section for steady condition.
        
        Considers that the flow is steady, return the forces and induction
        factors.

        Parameters
        ----------
        azimuth : float
            azimuthal angle, radians.
        pitch : float
            blade pitch angle, radians.
        wind: float
            wind velocity, m/s
        omega: float
            angular velocity, m/s
        args
            extra arguments
        angles : list or array
            list or array with the inflow angle (yaw, tilt), in radians.
        precone : float
            precone angle, in radians.
        elements
            slice or list of indexes of sections to consider.
            If default - slice(None) - considers all the sections,
        kwargs
            extra key arguments of the solve method.
        """
        
        if type(elements) is slice:
            sections_to_consider = self.rotor.sections[elements]
        else:
            sections_to_consider = [self.rotor.sections[i] for i in elements]
        number_sections = len(sections_to_consider)

        forces = np.zeros([number_sections,2])
        inductions = np.zeros([number_sections,2])

        # TODO: very very bad practice! replace this by obj attributes!
        pre_args = list(inspect.signature(self.pre).parameters.keys())
        for _k in ('args','kwargs'):
            if _k in pre_args: pre_args.remove(_k)
        
        _locals = locals()
        if len(pre_args) > 0:
            pre_inputs = {
                name:_locals[name] for name in pre_args if name in _locals
                }
        # loop for all sections
        for i, section in enumerate(sections_to_consider):
            velocity = tools.calculateVelocity(
                wind,omega,section.radius,azimuth,angles[0],angles[1],precone
            )
            if len(pre_args) > 0:
                # update loop evolving values
                if 'section' in pre_args: pre_inputs['section'] = section
                if 'velocity' in pre_args: pre_inputs['velocity'] = velocity
                pre_kwargs = self.pre(**pre_inputs)
                kwargs.update(pre_kwargs)
                
            fn, ft, a, aprime = self.solve(
                section,azimuth,pitch,*args,
                velocity=velocity,angles=angles,**kwargs
                )

            forces[i,:] = (fn,ft)
            inductions[i,:] = (a,aprime)
        
        return forces, inductions


    def dynamic(self,azimuth:list,pitch:list,wind:list,omega:list,*args,
               angles:list=[0.0,0.0],precone:float=0.0,elements=slice(None),
               **kwargs):
        """Dynamic solution for set of states.

        TODO: validate it!
        
        Call steady solution for set of evolving pitch and azimuths.
        Not really dynamic because the previous state is not considered?

        The arguments can either be lists or floats. If floats, they are
        considered constant and the value is copied for as many times as the
        number of steps. At least one of them must be a list or have length
        higher than 0.

        Parameters
        ----------
        azimuth: float or list
            azimuthal angle. See `bem.BaseBEM.steady`.
        pitch: float or list
            pitch angle. See `bem.BaseBEM.steady`.
        wind: float or list
            wind velocity. See `bem.BaseBEM.steady`.
        omega: float or list
            angular velocity. See `bem.BaseBEM.steady`.
        args
            reamining args
        kwargs
            remaining  key args

        """
        
        
        if type(elements) is slice:
            sections_to_consider = self.rotor.sections[elements]
        else:
            sections_to_consider = [self.rotor.sections[i] for i in elements]
        number_sections = len(sections_to_consider)
        
        # loop for all steps
        # select any input that is a list and make the others the same
        # if they are not
        inputs = dict(azimuth=azimuth,pitch=pitch,wind=wind,omega=omega)
        number_steps = 1
        for key, val in inputs.items():
            try:
                if len(val) > 1:
                    number_steps = max(len(val),number_steps)
            except TypeError:
                inputs[key] = [val]
        for key, val in inputs.items():
            if len(val) == 1: inputs[key] = val*number_steps

        forces = np.zeros((number_sections,number_steps,2))
        factors = np.zeros((number_sections,number_steps,2))
        
        # loop for all sections
        for i, section in enumerate(sections_to_consider):
            # loop for all steps
            for ii, (_azi, _pit, _wnd, _omg) in enumerate(zip(*inputs.values())):
                velocity = tools.calculateVelocity(
                    _wnd,_omg,section.radius,azimuth,
                    angles[0],angles[1],precone
                )
                forces[i,ii,0], forces[i,ii,1], factors[i,ii,0], factors[i,ii,1] = self.solve(
                    section,_azi,_pit,*args,
                    velocity=velocity,angles=angles,**kwargs
                    )
            for corr in self.corrections.__dict__.values():
                if hasattr(corr,'restart'): corr.restart()
        
        return forces, factors

    
    def cycle(self,pitch:float,wind,omega:float,*args,
              angles:list=[0.0,0.0],precone:float=0.0,
              N=1.0,dt:float=None,delta_phi:float=None,n_phi:int=36,
              **kwargs):
        """Solve for a given number of revolutions.
        
        Evaluate steady solution for several azimuths. Different from dynamic
        because the time is not considered.

        Parameters
        ----------
        pitch: float
            pitch angle. See `bem.BaseBEM.steady`.
        wind: float
            wind velocity. See `bem.BaseBEM.steady`.
        omega: float
            angular velocity. See `bem.BaseBEM.steady`.
        angles: list
            yaw and tilt angles. See `bem.BaseBEM.steady`.
        precone: float
            precone angle. See `bem.BaseBEM.steady`.
        N : float or int, optional
            number of revolutions. Default is 1.
        dt : float
            timesteps duration, s.
        delta_phi : float
            azimuthal step, rad.
        n_phi: int
            number of azimuthal steps. Optional, 36 by default.
        kwargs
            extra key arguments for call to `solve`
        
        Returns
        -------
        Arrays with evolution of forces, induction factors and azimuthal angle.
        The size of the arrays depends on the number of timesteps/azimuths and 
        number of elements.
        
        """

        if delta_phi is not None:
            pass
        elif dt is not None:
            delta_phi = dt*omega
        elif n_phi is not None:
            delta_phi = 2.0*np.pi*N/(n_phi - 1)
        else:
            raise ValueError(f'Timestep or azimuthal step must be defined')

        azimuths = np.arange(0.0,N*2*np.pi+delta_phi,delta_phi)
        n_phi = len(azimuths)

        # forces and inductions can have different size from what is defined
        # inside the steady solution (elements argument)
        forces, inductions = [], []

        # loop for all azimuths
        for i, azimuth in enumerate(azimuths):
            _forces, _inductions = self.steady(
                azimuth,pitch,wind,omega,*args,
                angles=angles,precone=precone,**kwargs
            )
            forces.append(_forces)
            inductions.append(_inductions)

        
        return np.array(forces), np.array(inductions), azimuths