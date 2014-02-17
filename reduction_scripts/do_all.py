"""
Do all reductions in the appropriate order.

Must have downloaded the data and set paths in paths.py first
"""
execfile('calibrate_C.py')
execfile('calibrate_Ku.py')
execfile('makecube_Lima_C.py')
execfile('makecube_Lima_Ku.py')
execfile('make_tau_cubes.py')
