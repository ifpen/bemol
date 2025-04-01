"""

Sample code for aligned flow.

"""

import os

import pandas as pd
import matplotlib.pyplot as plt

## uncoment this to if ModuleNotFounError
## add repo folder to PYTHONPATH, if running from the samples folder
# import sys
# sys.path.append(os.path.abspath(f'{__file__}/../..'))

import bemol

results_folder = 'results/aligned'
os.makedirs(results_folder,exist_ok=True)

wind = 15.06
omega = 44.5163679
rho = 1.191

corrections = (
    bemol.secondary.HubTipLoss.Prandtl,
    bemol.secondary.SkewAngle.Burton,
    bemol.secondary.YawModel.PittAndPeters,
    bemol.secondary.DynamicInflow.Dummy,
    bemol.secondary.TurbulentWakeState.Buhl,
    )

mexico = bemol.rotor.mexico

# all angles are null, flow completly aligned
skewAngle = 0.0
yawAngle = 0.0
azimuthAngle = 0.0

solver_uncoupled = bemol.ning.NingUncoupled(mexico,rho,corrections)
forces, _ = solver_uncoupled.steady(
    azimuthAngle,mexico.pitch,wind,omega
)
data = pd.DataFrame({'fn':forces[:,0],'ft': forces[:,1]})
data.to_csv(f'{results_folder}/results_aligned_uncoupled.csv',index=False)


solver_coupled = bemol.ning.NingCoupled(mexico,rho,corrections)
forces_coupled, _ = solver_coupled.steady(
    azimuthAngle,mexico.pitch,wind,omega,
    angles=[yawAngle,0.0],
    skew=0.0,
)
data = pd.DataFrame({'fn':forces_coupled[:,0],'ft':forces_coupled[:,1]})
data.to_csv(f'{results_folder}/results_aligned_coupled.csv',index=False)


lib_folder = os.path.abspath(os.path.dirname(bemol.__file__))
refCastor = pd.read_csv(
    f'{lib_folder}/rotors/mexico/ref/CASTOR.dat',
    index_col=None,comment='#',sep=' ',
    )

for i, name in enumerate(['Fn','Ft']):

    # change orientation for tangent
    factor = -1 if name == 'Ft' else 1

    fig, ax = plt.subplots(1,1,constrained_layout=True)
    ax.plot(refCastor['r'],refCastor[f'{name}15'],label='CASTOR')
    ax.plot(mexico.radius,factor*forces[:,i],label='Uncoupled BEM')
    ax.plot(mexico.radius,factor*forces_coupled[:,i],'--',label='Coupled BEM')
    ax.legend()
    ax.grid()
    ax.set_xlabel('radius, m')
    ax.set_ylabel(f'{name}, N/m')
    fig.savefig(f'{results_folder}/graph_aligned_force_{name}.png')
    plt.show() # comment if you dont want to show the figure


