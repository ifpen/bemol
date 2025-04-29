"""

Sample with comparison of yaw models.

"""


import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

## uncoment following lines if ModuleNotFounError
## add repo folder to PYTHONPATH, if running from the samples folder
import sys
sys.path.append(os.path.abspath(f'{__file__}/../..'))

import bemol


results_folder = 'results/yaw_models'
os.makedirs(results_folder,exist_ok=True)

wind = 15.06
omega = 44.5163679
rho = 1.191
number_revolutions = 1.0
tStep = 0.1
elements = [28]


azimuthAngle = 0.0
preconeAngle = 0.0
tiltAngle = 0.0
skewAngle = np.radians(30.)
yawAngle = skewAngle

mexico = bemol.rotor.mexico

base_corrections = [
    bemol.secondary.HubTipLoss.Prandtl,
    bemol.secondary.SkewAngle.Burton,
    bemol.secondary.TurbulentWakeState.Buhl,
]


## calculate solution for different yaw models
for yaw_model in ('Dummy','PittAndPeters','IFPEN'):
    
    corrections = base_corrections.copy()
    corrections.append(getattr(bemol.secondary.YawModel,yaw_model))

    solver = bemol.ning.NingUncoupled(mexico,rho,corrections)
    forces, _, azimuths = solver.cycle(
        mexico.pitchRated,wind,omega,angles=[yawAngle,tiltAngle],tStep=tStep,
        n_phi=180,N=number_revolutions,elements=elements,
    )

    azimuth_deg = np.degrees(azimuths)
    ## plot values
    plt.plot(
        azimuth_deg,forces[:,0,0],
        label=yaw_model
        )
    
    ## export data
    data = pd.DataFrame(
        {'azi':azimuth_deg,'fn':forces[:,0,0],'ft':forces[:,0,1]}
        )
    data.to_csv(f'{results_folder}/results_yaw_model_{yaw_model}.csv')


plt.xlabel('azimuth, deg')
plt.ylabel('normal lineic force, N/m')
plt.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=3)
plt.grid()
plt.savefig(f'{results_folder}/graph_yaw_models_azimuthal.png')
plt.show() # uncomment if you want to show the figure