
"""

Solution for 3 sections for rotor under yaw.

"""


import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import bemol


results_folder = 'results/yaw'
os.makedirs(results_folder,exist_ok=True)

wind = 15.06
omega = 44.5163679
rho = 1.191
number_revolutions = 4.0
tStep = 0.1
elements = [10, 18, 28]


azimuthAngle = 0.0
preconeAngle = 0.0
tiltAngle = 0.0
skewAngle = np.radians(30.)
yawAngle = skewAngle

mexico = bemol.rotor.mexico

corrections = (
    bemol.secondary.HubTipLoss.Prandtl,
    bemol.secondary.SkewAngle.Burton,
    bemol.secondary.TurbulentWakeState.Buhl,
    )


## solution for several revolutions, uncoupled Ning
solver_uncoupled = bemol.ning.NingUncoupled(mexico,rho,corrections)
forces_uncoupled, _, azimuths = solver_uncoupled.cycle(
    mexico.pitch,wind,omega,angles=[yawAngle,tiltAngle],tStep=tStep,
    n_phi=180,N=number_revolutions,elements=elements,
)

## solution for several revolutions, coupled Ning
solver_coupled = bemol.ning.NingCoupled(mexico,rho,corrections)
forces_coupled, _, _ = solver_coupled.cycle(
    mexico.pitch,wind,omega,angles=[yawAngle,tiltAngle],tStep=tStep,
    uInfty=wind,skew=skewAngle,
    n_phi=180,N=number_revolutions,elements=elements,
)

azimuth_deg = np.degrees(azimuths)
for i, elt in enumerate(elements):
    ## plot values
    plt.plot(
        azimuth_deg,forces_uncoupled[:,i,0],
        label=f'{elt}-th element (r = {mexico[i].radius:.4})',
        )
    plt.plot(
        azimuth_deg,forces_coupled[:,i,0],'--',
        label='Coupled BEM',
        )
    
    ## export data
    data = pd.DataFrame(
        {'azi':azimuth_deg,'fn':forces_uncoupled[:,i,0],'ft':forces_uncoupled[:,i,1]}
        )
    data.to_csv(f'{results_folder}/results_yaw_uncoupled_i{elt}.csv')
    data = pd.DataFrame(
        {'azi':azimuth_deg,'fn':forces_coupled[:,i,0],'ft':forces_coupled[:,i,1]}
        )
    data.to_csv(f'{results_folder}/results_yaw_coupled_i{elt}.csv')

plt.legend()
plt.savefig(f'{results_folder}/graph_yaw_azimuthal.png')
plt.close()