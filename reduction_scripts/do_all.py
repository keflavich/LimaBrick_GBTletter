# Define the paths for the various files
# outpath must include at least one formatting string, "%s"
# paths must include trailing slashes
global_paths ={'AGBT12B_221_path':'/Users/adam/observations/gbt/AGBT12B_221_01/',
               'AGBT14A_110_path':'/Users/adam/observations/gbt/AGBT14A_110_01/',
               'outpath':'/Users/adam/observations/gbt/%smap/',
              }

execfile('calibrate_C.py')
execfile('calibrate_Ku.py')
execfile('makecube_Lima_C.py')
execfile('makecube_Lima_Ku.py')
execfile('make_tau_cubes.py')
