"""

Convert the OpenFAST input data to bemol format.

"""

from pathlib import Path
import re

import numpy as np
import pandas as pd



hub_radius = 3.97

# path to the IEA-15-240-RWT repo
src_folder = Path('../../../../IEA-15-240-RWT/OpenFAST/IEA-15-240-RWT')

src_blade = pd.read_csv(
    src_folder/'IEA-15-240-RWT_AeroDyn15_blade.dat',
    skiprows=[0,1,2,3,5],sep='\s+'
    )

data = {}
data['radius'] = src_blade['BlSpn'] + hub_radius
data['twist'] = np.deg2rad(src_blade['BlTwist'])
data['chord'] = src_blade['BlChord']
data['airfoil'] = [f'AF{i:02}' for i in src_blade['BlAFID']]

pd.DataFrame(data,).to_csv(
    'blade.dat',sep='\t',index=False
)


for file in Path(src_folder/'Airfoils').glob('*Polar_*.dat'):
    
    with open(file,'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if 'Alpha' in line and 'Cl' in line:
                break
    src_polar = pd.read_csv(
        file,names=('alpha','Cl','Cd','Cm'),
        skiprows=i+2,sep='\s+'
        )
    foil_number = int(re.search(r'.*Polar_(\d+)\.dat',f.name).group(1))
    foil_name = f'AF{foil_number+1:02d}'
    polar = {}
    polar['AoA'] = src_polar['alpha']
    polar['Cl'] = src_polar['Cl']
    polar['Cd'] = src_polar['Cd']

    with open(f'airfoils/{foil_name}.foil','w') as f:
        f.write('# AoA [o], Cl [.], Cd [.]\n')
        pd.DataFrame(polar,).to_csv(
            f,sep=' ',header=False,index=False
        )