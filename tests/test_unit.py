
"""Unit tests"""

import os
import numpy as np

import bemol

lib_path = os.path.abspath(os.path.dirname(bemol.__file__))


def test_airfoil():
    """Test definition of the airfoil class."""
    airfoil = bemol.airfoil.BaseAirfoil(f'{lib_path}/rotors/mexico/airfoils/RISOE.foil')
    # just check if can returns single value and array
    assert airfoil.cl(0.1).size == 1
    assert type(airfoil._cl) == np.ndarray
    assert airfoil._cl.size > 1


def test_rotor():
    """Test definition of the rotor class."""

    # test definition of rotors
    assert 'mexico' in bemol.rotor.__dict__
    assert 'iea15mw' in bemol.rotor.__dict__

    rotor = bemol.rotor.mexico
    # check iterator
    for i, section in enumerate(rotor):
        assert rotor[i] is section

    # test definition of iea15mw rotor
    iea15 = bemol.rotor.iea15mw
    assert iea15.pitchRated == 0.06499140591234773



def test_secondary_module():
    """Test definition of secondary effects."""
    dummy_rotor = bemol.rotor.mexico

    # check if all the secondary models are selected when user is not defining them
    for corrections in (None,[],{}):
        model = bemol.bem.BaseBEM(dummy_rotor,corrections=corrections)
        assert hasattr(model.corrections,'hubTipLoss')
        assert hasattr(model.corrections,'skewAngle')
        assert hasattr(model.corrections,'dynamicInflow')
        assert hasattr(model.corrections,'yawModel')
        assert hasattr(model.corrections,'turbulentWakeState')

    # check if the Dummy is selected everytime
    model = bemol.bem.BaseBEM(dummy_rotor)
    assert model.corrections.hubTipLoss.__class__ is bemol.secondary.HubTipLoss.Dummy
    assert model.corrections.skewAngle.__class__ is bemol.secondary.SkewAngle.Dummy
    assert model.corrections.dynamicInflow.__class__ is bemol.secondary.DynamicInflow.Dummy
    assert model.corrections.yawModel.__class__ is bemol.secondary.YawModel.Dummy
    assert model.corrections.turbulentWakeState.__class__ is bemol.secondary.TurbulentWakeState.Dummy

    # define a model with a few custom corrections
    model = bemol.bem.BaseBEM(
        dummy_rotor,
        corrections=[bemol.secondary.HubTipLoss.Prandtl,bemol.secondary.SkewAngle.Burton]
        )
    assert model.corrections.hubTipLoss.__class__ is bemol.secondary.HubTipLoss.Prandtl
    assert model.corrections.skewAngle.__class__ is bemol.secondary.SkewAngle.Burton
    assert model.corrections.dynamicInflow.__class__ is bemol.secondary.DynamicInflow.Dummy
    assert model.corrections.yawModel.__class__ is bemol.secondary.YawModel.Dummy
    assert model.corrections.turbulentWakeState.__class__ is bemol.secondary.TurbulentWakeState.Dummy

    # same as before, but with dict input define a model with a few custom corrections
    model = bemol.bem.BaseBEM(
        dummy_rotor,
        corrections=dict(
            hubTipLoss=bemol.secondary.HubTipLoss.Prandtl,
            skewAngle=bemol.secondary.SkewAngle.Burton
            )
        )
    assert model.corrections.hubTipLoss.__class__ is bemol.secondary.HubTipLoss.Prandtl
    assert model.corrections.skewAngle.__class__ is bemol.secondary.SkewAngle.Burton
    assert model.corrections.dynamicInflow.__class__ is bemol.secondary.DynamicInflow.Dummy
    assert model.corrections.yawModel.__class__ is bemol.secondary.YawModel.Dummy
    assert model.corrections.turbulentWakeState.__class__ is bemol.secondary.TurbulentWakeState.Dummy

    # mix class and instance
    hub_corr = bemol.secondary.HubTipLoss.Prandtl()
    model = bemol.bem.BaseBEM(
        dummy_rotor,corrections=[hub_corr,bemol.secondary.SkewAngle.Burton]
        )
    assert model.corrections.hubTipLoss is hub_corr
    assert model.corrections.skewAngle.__class__ is bemol.secondary.SkewAngle.Burton

    # check if corrections of nested solvers are the same instance
    ## TOOD: maybe this is not always the wanted behavior!
    solver_uncoupled = bemol.ning.NingUncoupled(dummy_rotor,1.0,corrections)
    solver_coupled = bemol.ning.NingCoupled(dummy_rotor,1.0,corrections)
    assert solver_coupled.corrections.dynamicInflow is not solver_uncoupled.corrections.dynamicInflow
    assert solver_coupled.corrections.dynamicInflow is solver_coupled._uncoupled.corrections.dynamicInflow