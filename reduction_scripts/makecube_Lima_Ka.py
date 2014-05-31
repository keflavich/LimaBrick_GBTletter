import os
import astropy.io.fits as pyfits
import itertools
import sys
from sdpy import makecube,make_off_template,calibrate_map_scans
import numpy as np
from paths import outpath
# to ignore div-by-zero errors?
np.seterr(all='ignore')

cubename = os.path.join(outpath,'LimaBean_H2CO33_cube')
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
    makecube.add_file_to_cube(fn,
                              cubename+'.fits',nhits=cubename+'_nhits.fits',
                              chmod=True,
                              add_with_kernel=True,
                              kernel_fwhm=10./3600.,
                              velocityrange=[-250,250],
                              excludefitrange=[0,50],
                              diagnostic_plot_name=fn.replace('.fits','_data_scrubbed.png'),
                              smoothto=2)

os.system(os.path.join(outpath,'LimaBean_H2CO33_cube_starlink.sh'))

makecube.make_flats(cubename,vrange=[-20,60],noisevrange=[150,200])


sampler_feeds = {'A10': 1,
                 'A14': 2,
                 'B17': 1,
                 'B18': 1,
                 'B21': 2,
                 'B22': 2,
                 }


for cubename,restfreq,samplers in (
        #('LimaBean_H113a_cube', 4.497776e9, ('D34','D38')),
        #('LimaBean_H110a_cube', 4.874157e9, ('D33','D37')),
        ('LimaBean_H213CO33_cube', 27.55567, ["A10","A14"]),
        ('LimaBean_H2C18O33_cube', 26.33016, ["B18","B22"]),
        ('LimaBean_H62a_cube', 26.93916, ["B17","B21"]),
        #('LimaBean_H109a_cube', 5.008922e9, ('B18','B22')),
        #('LimaBean_H112a_cube', 4.61879e9, ('C26','C30')),
        #('LimaBean_OHF44_cube', 5.52344e9, ('A10','A14')),
        #('LimaBean_CH3NH2_cube', 5.19543e9, ('B21','B17'))
            ):

    cubename = os.path.join(outpath,cubename)

    makecube.generate_header(0.256, 0.0220, naxis1=100, naxis2=100, pixsize=10,
                             naxis3=500, cd3=1.0, clobber=True,
                             restfreq=restfreq, crval3=0.0, cunit3='km/s',
                             bmaj=28/3600., bmin=28/3600.)

    makecube.make_blank_images(cubename,clobber=True)

    #INLINE WHOA
    files = [x for session,scan1,scan2 in zip([2, 2, 3, 3, 3,  3,  3,  3,  3,"04","04","04","04"],
                                              [29,80,11,35,76, 102,148,174,215,6, 32,73,99],
                                              [79,96,35,75,102,142,174,214,229,32,72,99,113],)
             for x in [os.path.join(outpath, '14A_110_%s_%ito%i_%s_F%i.fits' %
                                    (session,scan1,scan2,samplers[ii],sampler_feeds[samplers[ii]]))
                       for ii in xrange(len(samplers))]]

    for fn in files:
        makecube.add_file_to_cube(fn,
                                  cubename+'.fits',
                                  nhits=cubename+'_nhits.fits',
                                  add_with_kernel=True,
                                  chmod=True,
                                  kernel_fwhm=10./3600.,
                                  velocityrange=[-250,250],
                                  excludefitrange=[0,50], smoothto=2)

    os.system('chmod +x %s_starlink.sh' % cubename)
    os.system('%s_starlink.sh' % cubename)
    makecube.make_flats(cubename,vrange=[-20,60],noisevrange=[150,200])


#import FITS_tools
from astropy.io import fits
from FITS_tools.cube_regrid import gsmooth_cube,spatial_smooth_cube

for cubename in ('LimaBean_H2CO33_cube', 'LimaBean_H213CO33_cube', 'LimaBean_H2C18O33_cube'):

    cubename = os.path.join(outpath,cubename)
    cube = fits.open(cubename+"_sub.fits")
    # OUT OF DATE: for 2-2....
    # kernel = ((2.5*60)**2 -  50.**2)**0.5 / sqrt(8*log(2)) = 60 arcsec
    # 60 arcsec / 15 "/pixel = 4
    cubesm2 = gsmooth_cube(cube[0].data, [5,4,4], use_fft=True, psf_pad=False, fft_pad=False)
    cubesm = spatial_smooth_cube(cube[0].data, kernelwidth=4, interpolate_nan=True)
    cube[0].data = cubesm
    cube.writeto(cubename+"_sub_smoothtoCband.fits",clobber=True)
    cube[0].data = cubesm2
    cube.writeto(cubename+"_sub_smoothtoCband_vsmooth.fits",clobber=True)

    # etamb = 56-64% * 1.37 from GTBpg
    makecube.make_taucube(cubename,
                          cubename+"_continuum.fits",
                          etamb=1.37*0.6,
                          tex=2.0,
                          suffix="_sub_smoothtoCband_vsmooth.fits",
                          outsuffix="_smoothtoCband_vsmooth.fits",
                          TCMB=2.7315+0.4)

    makecube.make_taucube(cubename,
                          cubename+"_continuum.fits",
                          etamb=1.37*0.6,
                          tex=2.0,
                          suffix="_sub_smoothtoCband.fits",
                          outsuffix="_smoothtoCband.fits",
                          TCMB=2.7315+0.4)

    makecube.make_taucube(cubename,
                          cubename+"_continuum.fits",
                          etamb=1.37*0.6,
                          tex=2.0,
                          suffix="_sub.fits",
                          outsuffix=".fits",
                          TCMB=2.7315+0.4)

