import os
datapath = '/Users/adam/observations/gbt/LimaBeanmap/'
outpath = '/Users/adam/observations/gbt/LimaBeanmap/'
root = '/Users/adam/work/h2co/limabean/'
figpath = os.path.join(root, 'figures/')
regpath = os.path.join(root, 'regions/')

def rpath(x, p=regpath):
    return os.path.join(p,x)

def dpath(x, p=datapath):
    return os.path.join(p,x)
