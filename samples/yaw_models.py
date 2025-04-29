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
elements = [24] # node close to the experimental data (r/R = 0.82)


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



## plot experimental data
# r/R = 0.82
# Final report of IEA-29 Phase 3, page 37, Figure 4.6 (d)
# scanned with webPlotDigitilizer
# data-points every 10 deg from 0
ref_azimuth = np.arange(0,360,10)
ref_data = [
    417.51361161524494, 413.5208711433756, 409.52813067150623, 407.7132486388384,
    406.9872958257712, 407.3502722323048, 409.52813067150623, 412.43194192377484,
    416.4246823956442, 420.05444646098, 424.7731397459164, 430.21778584392007,
    435.2994555353901, 441.10707803992733, 447.6406533575317, 453.0852994555353,
    460.3448275862068, 467.9673321234119, 475.9528130671506, 484.30127041742276,
    492.649727767695, 499.9092558983665, 504.9909255898366, 508.9836660617059,
    509.709618874773, 507.53176043557164, 503.5390199637023, 495.55353901996364,
    487.20508166969137, 477.0417422867513, 466.87840290381115, 455.9891107078039,
    445.8257713248638, 436.7513611615244, 428.76588021778576, 422.595281306715,
]
plt.plot(
    ref_azimuth,ref_data,'sk',markerfacecolor='none',
    label='exp.'
    )


plt.xlabel('azimuth, deg')
plt.ylabel('normal lineic force, N/m')
plt.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=3)
plt.grid()
plt.savefig(f'{results_folder}/graph_yaw_models_azimuthal.png')
# plt.show() # uncomment if you want to show the figure