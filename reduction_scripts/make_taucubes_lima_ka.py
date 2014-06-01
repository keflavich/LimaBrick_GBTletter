from sdpy import makecube
from paths import outpath
from astropy.io import fits
from FITS_tools.cube_regrid import gsmooth_cube,spatial_smooth_cube

for cubename in ('LimaBean_H2CO33_cube', 'LimaBean_H213CO33_cube',
                 'LimaBean_H2C18O33_cube'):

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
                          continuum=0.0,
                          etamb=1.37*0.6,
                          tex=0.0,
                          suffix="_sub_smoothtoCband_vsmooth.fits",
                          outsuffix="_smoothtoCband_vsmooth.fits",
                          TCMB=2.7315)

    makecube.make_taucube(cubename,
                          continuum=0.0,
                          etamb=1.37*0.6,
                          tex=0.0,
                          suffix="_sub_smoothtoCband.fits",
                          outsuffix="_smoothtoCband.fits",
                          TCMB=2.7315)

    makecube.make_taucube(cubename,
                          continuum=0.0,
                          etamb=1.37*0.6,
                          tex=0.0,
                          suffix="_sub.fits",
                          outsuffix=".fits",
                          TCMB=2.7315)

