"""

Sample code for the calculation of forces on IEA 15MW blade, aligned flow.

"""

import os

import pandas as pd
import matplotlib.pyplot as plt

## uncoment this to if ModuleNotFounError
## add repo folder to PYTHONPATH, if running from the samples folder
# import sys
# sys.path.append(os.path.abspath(f'{__file__}/../..'))

import bemol

results_folder = 'results/iea15mw'
os.makedirs(results_folder,exist_ok=True)

turbine = bemol.rotor.iea15mw
wind = turbine.windRated
omega = turbine.omegaRated
pitch = turbine.pitchRated
precone = turbine.precone
tilt = turbine.tilt
rho = 1.225

corrections = (
    bemol.secondary.HubTipLoss.Prandtl,
    bemol.secondary.TurbulentWakeState.Buhl,
    )

# all angles are null, flow completly aligned
skewAngle = 0.0
yawAngle = 0.0
azimuthAngle = 0.0

solver = bemol.ning.NingCoupled(turbine,rho,corrections)
forces, _ = solver.steady(
    azimuthAngle,pitch,wind,omega,
    angles=[0.0,tilt],precone=precone,
    skew=0.0,
)
data = pd.DataFrame({'Fn':forces[:,0],'Ft':forces[:,1]})
data.to_csv(f'{results_folder}/results_aligned_coupled.csv',index=False)

lib_folder = os.path.abspath(os.path.dirname(bemol.__file__))
# get reference data for plot
ref = {}
for solver in ('AeroDeeP','CASTOR'):
    ref[solver] = pd.read_csv(
        f'{lib_folder}/rotors/iea15mw/ref/data_{solver}.csv',
        index_col=None,comment='#',sep=',',
        )
    

fig, axs = plt.subplots(1,2,constrained_layout=True)

for name, ax in zip(['Fn','Ft'],axs):
    
    # change orientation for tangent
    factor = -1 if name == 'Ft' else 1

    ax.plot(ref['AeroDeeP']['radius'],factor*ref['AeroDeeP'][name],
            '-ob',linewidth=1.0,markersize=3,label='AeroDeeP (BEM)')
    ax.plot(ref['CASTOR']['radius'],factor*ref['CASTOR'][name],
            '-sk',linewidth=1.0,markersize=3,label='CASTOR (FVW)')
    
    ax.plot(turbine.radius,factor*data[name],'-',label='bemol (coupled BEM)')
    ax.grid()
    ax.set_xlabel('radius, m')
    ax.set_ylabel(f'{name}, N/m')

ax.legend() # legend only on the right plot
fig.savefig(f'{results_folder}/graph_iea15mw_forces.png')
#plt.show() # uncomment if you want to show the figure


