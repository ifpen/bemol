
"""

Compare if results are strictly the same

Update the folder with the reference solutions (not stored in the git).

"""

import shutil
import os
import importlib
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pytest

import bemol


lib_path = os.path.abspath(os.path.dirname(bemol.__file__))


@pytest.fixture(autouse=True,)
def change_test_dir(request,monkeypatch):
    # run tests in the tests folder, not path test was called
    monkeypatch.chdir(request.fspath.dirname)


@pytest.fixture(scope='module',autouse=True,)
def reference_folder(request):
    # define the reference folder based on current folder
    test_path = request.path.parent
    reference_path = test_path / 'ref' / 'results_ref_2025-03-25'
    return Path( reference_path ).resolve()



@pytest.fixture(scope='module')
def folder(request):
    """Remove previous results."""
    test_folder = Path(request.fspath.dirname).resolve()
    new_data_folder = test_folder / 'results'
    shutil.rmtree(new_data_folder,ignore_errors=True)
    new_data_folder.mkdir(exist_ok=True)
    return new_data_folder


@pytest.fixture(autouse=True)
def collect():
    """Close all figures."""
    plt.close('all')


@pytest.fixture
def run_aligned():
    # TODO: use import to run the case, make it less hacky
    _ = importlib.import_module('aligned')


@pytest.fixture
def run_yaw():
    # TODO: use import to run the case, make it less hacky
    _ = importlib.import_module('yaw')


@pytest.fixture
def run_pitch_maneuver():
    # TODO: use import to run the case, make it less hacky
    _ = importlib.import_module('pitch_maneuver')



def test_airfoil(folder):
    """Test airfoil with dynamic stall correction.
    
    Not really checking anything, just plot the polars.
    """

    dataFile = f'{lib_path}/rotors/mexico/airfoils/RISOE.foil'
    myFoil = bemol.airfoil.DynStallAirfoil(dataFile)
    airfoil_folder = folder / 'airfoil'
    airfoil_folder.mkdir(exist_ok=True)

    aoas_deg = np.arange(-15., 30., 0.1)
    aoas_rad = np.radians(aoas_deg)

    plt.plot(aoas_deg,myFoil.cl(aoas_rad),label='$C_L$')
    plt.plot(aoas_deg,myFoil.cl0Slope*(aoas_rad - myFoil.aoaZero),label='slope')
    plt.plot(aoas_deg,myFoil.splineFullySeparatedLift(aoas_rad),label='separated')
    plt.plot(aoas_deg,myFoil.splineStaticAttachment(aoas_rad),label='static attach.')
    plt.scatter(np.rad2deg(myFoil.aoaZero),myFoil.cl(myFoil.aoaZero),label='AoA zero')
    
    plt.xlabel('angle of attack, deg')
    plt.ylabel('lift coefficient')
    plt.legend()
    plt.grid()
    plt.ylim(-2.0,2.0)

    plt.savefig( airfoil_folder / 'graph_airfoil_dynStall_polars.png' )



@pytest.mark.parametrize('run',('run_aligned','run_yaw','run_pitch_maneuver'))
def test_compare(folder,run,request,reference_folder):
    request.getfixturevalue(run)
    name_run = run.replace('run_','')
    results_folder = folder / name_run
    for new_filename in results_folder.glob('*.csv'):
        filename = new_filename.name
        new_data =  pd.read_csv(new_filename)

        ref_filename = reference_folder / filename
        ref_data = pd.read_csv(ref_filename)

        pd.testing.assert_frame_equal(ref_data,new_data)