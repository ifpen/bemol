"""

Blade pitch maneuver

"""

import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

## uncoment following lines if ModuleNotFounError
## add repo folder to PYTHONPATH, if running from the samples folder
# import sys
# sys.path.append(os.path.abspath(f'{__file__}/../..'))

import bemol


results_folder = 'results/pitch_maneuver'
os.makedirs(results_folder,exist_ok=True)

wind = 15.06
omega = 44.5163679
rho = 1.191

preconeAngle = 0.0
tiltAngle = 0.0
skewAngle = 0.0
yawAngle = skewAngle

tStep = 0.1

# solution for a single section at a given azimuth angle
azimuthAngle = 0.
iElement = 28

mexico = bemol.rotor.mexico

corrections = [
    bemol.secondary.HubTipLoss.Prandtl,
    bemol.secondary.SkewAngle.Burton,
    bemol.secondary.YawModel.Dummy,
    bemol.secondary.DynamicInflow.Knudsen,
    bemol.secondary.TurbulentWakeState.Buhl,
]

times = np.arange(0.0,50.0,tStep)
pitchs = np.zeros(len(times))
pitchs[times > 25.0] = np.radians(4.)

section = mexico[iElement]

Fns_Coupled = []
Fts_Coupled = []

Fns_Uncoupled = []
Fts_Uncoupled = []

UxPrime = wind*np.cos(preconeAngle)
UyPrime = omega*section.radius*np.cos(preconeAngle)

solver_uncoupled = bemol.ning.NingUncoupled(mexico,rho,corrections)
solver_coupled = bemol.ning.NingCoupled(mexico,rho,corrections)

UxUncoupled, UyUncoupled = bemol.tools.calculateVelocity(
    wind,omega,section.radius,azimuthAngle,yawAngle,tiltAngle,preconeAngle
    )

for (t, pitch) in zip(times, pitchs):

    ## uncoupled BEM solution
    Fn, Ft, _, _ = solver_uncoupled.solve(
        section,azimuthAngle,pitch,
        velocity=[UxUncoupled,UyUncoupled],angles=[yawAngle,0.0],
        tStep=tStep,
        )
    Fns_Uncoupled.append(Fn)
    Fts_Uncoupled.append(Ft)

    ## coupled BEM solution
    Fn, Ft, _, _ = solver_coupled.solve(
        section,azimuthAngle,pitch,
        velocity=[UxUncoupled,UyUncoupled],
        angles=[yawAngle,0.0],
        uInfty=wind,
        UxPrime=UxPrime,UyPrime=UyPrime,
        skew=skewAngle,
        tStep=tStep)
    Fns_Coupled.append(Fn)
    Fts_Coupled.append(Ft)

plt.plot(times,Fns_Uncoupled,label='uncoupled BEM')
plt.plot(times,Fns_Coupled,label='coupled BEM')
plt.xlabel('time, s')
plt.ylabel('normal lineic force, N/m')
plt.grid()
plt.legend()
plt.savefig(f'{results_folder}/graph_pitch_maneuver.png')
# plt.show() # uncomment if you want to show the figure


data = pd.DataFrame({'time':times,'fn':Fns_Uncoupled,'ft': Fts_Uncoupled})
data.to_csv(f'{results_folder}/results_pitch_maneuver_uncoupled.csv',index=False)

data = pd.DataFrame({'time':times,'fn':Fns_Coupled,'ft': Fts_Coupled})
data.to_csv(f'{results_folder}/results_pitch_maneuver_coupled.csv',index=False)