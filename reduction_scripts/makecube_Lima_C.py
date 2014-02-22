import os
import astropy.io.fits as pyfits
from astropy.io import fits
import itertools
import sys
from sdpy import makecube
import numpy as np
np.seterr(all='ignore')
from paths import outpath

cubename=outpath+'LimaBean_H2CO11_cube'
#cubename_discrete='LimaBean_H2CO11_cube_discrete'
# 15' x 12 '
#makecube.generate_header(0.256, 0.0220, naxis1=24, naxis2=24, pixsize=60,
#        naxis3=2400, cd3=1.0, clobber=True, restfreq=4.8296594e9)
makecube.generate_header(0.256, 0.0220, naxis1=100, naxis2=100, pixsize=15,
                         naxis3=800, cd3=1.0, clobber=True, restfreq=4.8296594e9)
makecube.make_blank_images(cubename,clobber=True)
#makecube.make_blank_images(cubename_discrete,clobber=True)

files = [outpath+'12B_221_6to21_A13_F1.fits',
         outpath+'12B_221_6to21_A9_F1.fits',
         # there was one noisy scan, so the whole damned thing gets commented out
         #outpath+'12B_221_22to32_A13_F1.fits',
         #outpath+'12B_221_22to32_A9_F1.fits',
         ]

for fn in files:
    makecube.add_file_to_cube(fn,
                              cubename+'.fits', nhits=cubename+'_nhits.fits',
                              wcstype='V', 
                              chmod=True,
                              add_with_kernel=True,
                              kernel_fwhm=90./3600.,
                              velocityrange=[-400, 400],
                              excludefitrange=[-225, 250], 
                              diagnostic_plot_name=fn.replace('.fits','_data_scrubbed.png'),
                              smoothto=2)

fstr = outpath+'12B_221_{0}to{1}_{2}_F1.fits'
cn_fmt = cubename+"_{0}_sr{1}-{2}"
for sr in ([6,21],[22,32]):
    for sampler in ('A9','A13'):
        cn = cn_fmt.format(sampler,sr[0],sr[1])
        fn = fstr.format(sr[0],sr[1],sampler)
        makecube.make_blank_images(cn,clobber=True)
        makecube.add_file_to_cube(fn,
                                  cn+'.fits',nhits=cn+'_nhits.fits',wcstype='V',
                                  add_with_kernel=True,
                                  chmod=True,
                                  kernel_fwhm=90./3600.,
                                  velocityrange=[-400,400],excludefitrange=[-225,250],
                                  smoothto=2)
    c1 = fits.getdata(cn_fmt.format('A9',sr[0],sr[1])+".fits")
    c2 = fits.getdata(cn_fmt.format('A13',sr[0],sr[1])+".fits")
    hdr = fits.getheader(cn+".fits")
    poldiff = c1-c2
    hdu = fits.PrimaryHDU(data=poldiff,header=hdr)
    hdu.writeto(cn_fmt.format("poldiff",sr[0],sr[1])+".fits",clobber=True)



os.system(outpath+'LimaBean_H2CO11_cube_starlink.sh')

makecube.make_flats(cubename,vrange=[-20,60],noisevrange=[250,300])
makecube.make_taucube(cubename,cubename+"_continuum.fits",etamb=0.98)


cubename=outpath+'LimaBean_H2CO11_cube_scan2'
files = [outpath+'12B_221_22to32_A13_F1.fits',
         outpath+'12B_221_22to32_A9_F1.fits',
         ]
makecube.generate_header(0.256, 0.0220, naxis1=100, naxis2=100, pixsize=15,
                         naxis3=800, cd3=1.0, clobber=True, restfreq=4.8296594e9)
makecube.make_blank_images(cubename,clobber=True)
for fn in files:
    makecube.add_file_to_cube(fn,
                              cubename+'.fits',
                              nhits=cubename+'_nhits.fits',
                              wcstype='V',
                              add_with_kernel=True,
                              kernel_fwhm=90./3600.,
                              velocityrange=[-400,400],excludefitrange=[-225,250],
                              diagnostic_plot_name=fn.replace('.fits','_data_scrubbed.png'),
                              chmod=True,
                              smoothto=2)
os.system(outpath+'LimaBean_H2CO11_cube_scan2_starlink.sh')

makecube.make_flats(cubename,vrange=[-20,60],noisevrange=[250,300])
makecube.make_taucube(cubename,cubename+"_continuum.fits",etamb=0.98)

# SAMPLERS:
# [('A10', 5500000000.0),
#  ('A13', 4829659400.0),
#  ('A14', 5500000000.0),
#  ('A9', 4829659400.0),
#  ('B17', 5200000000.0),
#  ('B18', 5008992000.0),
#  ('B21', 5200000000.0),
#  ('B22', 5008992000.0),
#  ('C25', 4593089000.0),
#  ('C26', 4630000000.0),
#  ('C29', 4593089000.0),
#  ('C30', 4630000000.0),
#  ('D33', 4874157000.0),
#  ('D34', 4497776000.0),
#  ('D37', 4874157000.0),
#  ('D38', 4497776000.0)]
#  OH F0-1 4.66024 
#  OH F4-4 5.52344 
#  OH F4-3 5.54704 


for cubename,restfreq,samplers in (
        ('LimaBean_H113a_cube', 4.497776e9, ('D34','D38')),
        ('LimaBean_H110a_cube', 4.874157e9, ('D33','D37')),
        ('LimaBean_H213CO_cube', 4.593089e9, ('C25','C29')),
        ('LimaBean_H109a_cube', 5.008922e9, ('B18','B22')),
        ('LimaBean_H112a_cube', 4.61879e9, ('C26','C30')),
        ('LimaBean_OHF44_cube', 5.52344e9, ('A10','A14')),
        ('LimaBean_CH3NH2_cube', 5.19543e9, ('B21','B17'))
       ):

    cubename = outpath+cubename


    #makecube.generate_header(0.256, 0.0220, naxis1=24, naxis2=24, pixsize=60,
    #        naxis3=2400, cd3=1.0, clobber=True, restfreq=restfreq)
    makecube.generate_header(0.256, 0.0220, naxis1=100, naxis2=100, pixsize=15,
                             naxis3=800, cd3=1.0, clobber=True,
                             restfreq=restfreq)
    makecube.make_blank_images(cubename,clobber=True)

    files = [outpath+'12B_221_6to21_%s_F1.fits' % samplers[0],
             outpath+'12B_221_6to21_%s_F1.fits' % samplers[1],
             # there was one noisy scan, so the whole damned thing gets commented out
             #outpath+'12B_221_22to32_%s_F1.fits' % samplers[0],
             #outpath+'12B_221_22to32_%s_F1.fits' % samplers[1],]
             ]
    for fn in files:
        makecube.add_file_to_cube(fn,
                                  cubename+'.fits',
                                  nhits=cubename+'_nhits.fits', wcstype='V',
                                  add_with_kernel=True, kernel_fwhm=90./3600.,
                                  chmod=True,
                                  diagnostic_plot_name=fn.replace('.fits','_data_scrubbed.png'),
            velocityrange=[-400,400],excludefitrange=[-150,225])

    os.system('chmod +x %s_starlink.sh' % cubename)
    os.system('%s_starlink.sh' % cubename)
    makecube.make_flats(cubename,vrange=[-20,60],noisevrange=[250,300])
    if 'CO' in cubename:
        makecube.make_taucube(cubename,cubename+"_continuum.fits",etamb=0.98)



try:
    import pyspeckit
    #pyspeckit.wrappers.cube_fit('LimaBean_H2CO11_cube_smooth.fits','LimaBean_H2CO11_parfits.fits','LimaBean_H2CO11_noise.fits',fittype='formaldehyde',absorption=True,signal_cut=2)


    from pylab import *
    def some_plots(xx,yy,dolegend=False,dohalpha=False):
        c13 = h213cocube.get_spectrum(xx,yy)
        c12 = cube.get_spectrum(xx,yy)
        c12.baseline(exclude=[-225,200],order=5)
        c13.baseline(exclude=[-225,200],order=5)
        c13.plotter(label='H$_{2}$$^{13}$CO',axis=gca())
        c12.plotter(axis=c13.plotter.axis,clear=False,color='b',label="H$_{2}$CO")
        (c13*6).plotter(label='6$\\times$H$_{2}$$^{13}$CO',axis=gca(),color='r',clear=False)
        if dolegend:
            legend(loc='best')
        if dohalpha:
            halpha = h110acube.get_spectrum(xx,yy)
            halpha.baseline(exclude=[-225,200],order=5)
            (halpha*5).plotter(axis=c13.plotter.axis,clear=False,color='darkgreen')
        gca().set_xlim(-100,100)
        draw()
            

    # some plots
    if not 'do_some_plots' in locals():
        do_some_plots=False
        dohalpha=True
    if do_some_plots:
        subplots_adjust(wspace=0,hspace=0,left=0.05,right=0.95,bottom=0.05,top=0.95)
        cube = pyspeckit.Cube(outpath+'LimaBean_H2CO11_cube_sub.fits')
        h213cocube = pyspeckit.Cube(outpath+'LimaBean_H213CO_cube_sub.fits')
        h110acube = pyspeckit.Cube(outpath+'LimaBean_H110a_cube_sub.fits')

        figure(1); clf()
        nx,ny = 5,5
        for ii,(spy,spx) in enumerate(itertools.product(range(nx),range(ny))):
            sp = subplot(nx,ny,ii+1)
            sp.zorder = nx-spx+ny*spy
            some_plots(spx+9,spy+9,dolegend=(ii==24),dohalpha=True)
            if dohalpha:
                sp.set_ylim(-2.5,1)
            else:
                sp.set_ylim(-2.5,0.5)
            if spx > 0: 
                sp.set_yticks([])
            if spy < ny-1: 
                sp.set_xticks([])
            
except ImportError:
    print "Pyspeckit is required to do plots"
