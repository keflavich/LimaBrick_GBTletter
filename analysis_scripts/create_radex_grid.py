import os
import shutil

print("RADEX grid creation requires that there is a RADEX executable on your "
      "PATH.  If there is not, it will crash in an ugly way.  The RADEX "
      "executable should be of the LVG variety.  Try using pyradex to install "
      "if you don't have it")

print("Warning: This is running a large grid of RADEX models.  "
      "It only has to run once, but it can take a long time, "
      "especially on a laptop.  The data it creates will be "
      "~100 MB")

if not os.path.isdir('radex/faure'):
    os.mkdir('radex/faure')

if not os.path.exists('radex/o-h2co-h2_faure.dat'):
    from astroquery import lamda
    data = lamda.query('oh2co-h2', return_datafile=True)
    datapath = 'radex/o-h2co-h2_faure.dat'
    with open(datapath,'w') as out:
        out.writelines([d+"\n" for d in data])

if not os.path.exists('radex/faure/1-1_2-2_XH2CO_fixed_faure.dat'):
    try:
        os.chdir('radex')
        execfile('radex_grid_faure_fixedX.py')
        shutil.move('1-1_2-2_XH2CO_fixed_faure.dat','faure/')
        shutil.move('1-1_3-3_XH2CO_fixed_faure.dat','faure/')
        shutil.move('2-2_3-3_XH2CO_fixed_faure.dat','faure/')
    except:
        raise
    finally:
        os.chdir('..')
