# rotors

Reference rotors.

## How to add new rotors

A rotor must include (also see `rotor.py`):

- blade.dat: definition of the blade (radial distribution of twist, chord and airfoils)
- rotor.yml: file with definition of the rotor parameters (number of blades, etc...)
- airfoils: polar files (`*.foil`) with the force coefficients for each airfoil
  present in the blade.
- ref: reference result values. Optional.

The rotor is dynamically read by the bemol lib and available to use.