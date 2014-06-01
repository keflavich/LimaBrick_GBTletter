import os
import astropy.io.fits as pyfits
import itertools
import sys
from sdpy import makecube,make_off_template,calibrate_map_scans
import numpy as np
from paths import outpath
# to ignore div-by-zero errors?
np.seterr(all='ignore')

for feed in (1,2):
    cubename = os.path.join(outpath,'LimaBean_H2CO33_cube_feed{0}'.format(feed))
    # 8' x 8'
    makecube.generate_header(0.256, 0.0220, naxis1=100, naxis2=100, pixsize=10,
                             naxis3=500, cd3=1.0, clobber=True, restfreq=28.9748e9,
                             crval3=0.0, cunit3='km/s',
                             bmaj=28/3600., bmin=28/3600.)
    makecube.make_blank_images(cubename,clobber=True)

    files = [os.path.join(outpath,x) for x in
             ['14A_110_2_29to79_A9_F1.fits',
              '14A_110_2_29to79_A13_F2.fits',
              '14A_110_2_80to96_A9_F1.fits',
              '14A_110_2_80to96_A13_F2.fits',
              '14A_110_3_11to35_A9_F1.fits',
              '14A_110_3_11to35_A13_F2.fits',
              '14A_110_3_35to75_A9_F1.fits',
              '14A_110_3_35to75_A13_F2.fits',
              '14A_110_3_76to102_A9_F1.fits',
              '14A_110_3_76to102_A13_F2.fits',
              '14A_110_3_102to142_A9_F1.fits',
              '14A_110_3_102to142_A13_F2.fits',
              '14A_110_3_148to174_A9_F1.fits',
              '14A_110_3_148to174_A13_F2.fits',
              '14A_110_3_174to214_A9_F1.fits',
              '14A_110_3_174to214_A13_F2.fits',
              '14A_110_3_215to229_A9_F1.fits',
              '14A_110_3_215to229_A13_F2.fits',
              '14A_110_04_6to32_A9_F1.fits',
              '14A_110_04_6to32_A13_F2.fits',
              '14A_110_04_32to72_A9_F1.fits',
              '14A_110_04_32to72_A13_F2.fits',
              '14A_110_04_73to99_A9_F1.fits',
              '14A_110_04_73to99_A13_F2.fits',
              '14A_110_04_99to113_A9_F1.fits',
              '14A_110_04_99to113_A13_F2.fits',
              ]
             ]

    for fn in files:
        if 'F{0}'.format(feed) not in fn:
            continue
        makecube.add_file_to_cube(fn,
                                  cubename+'.fits',nhits=cubename+'_nhits.fits',
                                  chmod=True,
                                  add_with_kernel=True,
                                  kernel_fwhm=10./3600.,
                                  velocityrange=[-250,250],
                                  excludefitrange=[0,50],
                                  diagnostic_plot_name=fn.replace('.fits','_data_scrubbed.png'),
                                  smoothto=2)

    os.system(os.path.join(outpath,cubename+'_starlink.sh'))

    makecube.make_flats(cubename,vrange=[-20,60],noisevrange=[150,200])
