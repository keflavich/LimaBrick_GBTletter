# Define the paths for the various files
# outpath must include at least one formatting string, "%s"
# paths must include trailing slashes
import os

AGBT12B_221_path = '/Users/adam/observations/gbt/AGBT12B_221_01/'
AGBT14A_110_path = '/Users/adam/observations/gbt/AGBT14A_110_01/'
AGBT14A_110_2_path = '/Users/adam/observations/gbt/AGBT14A_110_02/'
AGBT14A_110_3_path = '/Users/adam/observations/gbt/AGBT14A_110_3/'
AGBT14A_110_4_path = '/Users/adam/observations/gbt/AGBT14A_110_04/'
datapath = outpath = '/Users/adam/observations/gbt/LimaBeanmap/'
source_root = root = '/Users/adam/work/h2co/limabean/'
modelpath = os.path.join(source_root, 'models/')
figpath = os.path.join(root, 'figures/')
regpath = os.path.join(root, 'regions/')

def rpath(x, p=regpath):
    return os.path.join(p,x)

def dpath(x, p=datapath):
    return os.path.join(p,x)

def mpath(x, modelpath=modelpath):
    return os.path.join(modelpath, x)
