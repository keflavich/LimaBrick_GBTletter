"""
This code is to remake the tau cubes using the C-band image from Casey Law's
survey, which contains more reliable estimates of the continuum
"""
import sdpy
import numpy as np
from astropy.io import fits
import astropy.units as u
import FITS_tools
from paths import outpath
import os
from astropy import convolution

datapath = outpath

fivecubes = ['LimaBean_H213CO_cube_sub.fits',
             'LimaBean_H2CO11_cube_sub.fits',
             ]

# 14.488 rounds to 15 in this case.
fifteencubes = ['LimaBean_H2CO22_cube_sub.fits',
                'LimaBean_H2CO22_cube_sub_smoothtoCband.fits',
                'LimaBean_H213CO22_cube_sub.fits',
                #'LimaBean_H213CO22_cube_sub_smoothtoCband.fits',
                'LimaBean_H2C18O22_cube_sub.fits',
                #'LimaBean_H2C18O22_cube_sub_smoothtoCband.fits',
                ]

# TODO: Where can these be got?
ccontfile = datapath+'GCCBand_lb.2.fits'

# First do with Casey Law's images as the background
if os.path.exists(ccontfile):
    for fn in fivecubes:
        reproj_contfile = datapath+"GCCBand_reproj.fits"

        header = FITS_tools.strip_headers.flatten_header(fits.getheader(datapath+fivecubes[0]))

        ccont = FITS_tools.project_to_header(ccontfile, header) * u.Jy

        cbfreq = 4.829 * u.GHz
        gbbeam_5ghz = 1.22 * ((cbfreq.to(u.m, u.spectral()) / (100*u.m)) * u.rad).decompose()
        fwhm = np.sqrt(8*np.log(2))
        bmaj=4.2413E-02*u.deg
        bmin=4.2413E-02*u.deg
        ktojy5ghz = (1*u.K).to(u.Jy,u.brightness_temperature((2*np.pi*(bmaj*bmin/fwhm**2)), cbfreq))

        cont_K = (ccont/ktojy5ghz).value
        outfile = fits.PrimaryHDU(data=cont_K, header=header)
        outfile.writeto(reproj_contfile,clobber=True)

        continuum = cont_K
        suffix = '_claw.fits'

    for tex in [1,1.5,2]:
        sdpy.makecube.make_taucube(datapath+fn.replace("_sub.fits",""),
                                   continuum,
                                   tex=tex,
                                   linefreq=4.82966*u.GHz,
                                   outsuffix="_tex%s%s" % (tex,suffix),
                                   etamb=0.98)

print("For our own continuum data, there is an offset.")
print("Based on analysis shown in continuum_cband_law_comparison.py, we add 5 K"
      " to the background temperature to balance out an offset")

for fn in fivecubes:
    continuum = fits.getdata(datapath+fn.replace("_sub.fits","_continuum.fits")) + 5
    suffix = '.fits'

    for tex in [1,1.5,2]:
        sdpy.makecube.make_taucube(datapath+fn.replace("_sub.fits",""),
                                   continuum,
                                   tex=tex,
                                   linefreq=4.82966*u.GHz,
                                   outsuffix="_tex%s%s" % (tex,suffix),
                                   etamb=0.98)

# Now do the 2-2 14.5 GHz ("fifteen") lines

repl = lambda x: x.replace("_sub.fits","_continuum.fits").replace("_sub_smoothtoCband.fits","_continuum_smoothtoCband.fits")

continuum22 = fits.getdata(datapath+'LimaBean_H2CO22_cube_continuum.fits') + 0.4
continuum22smooth = convolution.convolve(continuum22,convolution.Gaussian2DKernel(4))

suffix = '.fits'
for fn in fifteencubes:
    for tex in [1.5,2,2.5]:
        filepath = datapath+fn.replace("_sub.fits","")
        if 'smoothtoCband' in fn:
            filepath = filepath.replace("_sub_smoothtoCband.fits","").replace(".fits","")
            sdpy.makecube.make_taucube(filepath,
                                       continuum22smooth if 'smoothtoCband' in fn else continuum22,
                                       tex=tex,
                                       linefreq=14.488*u.GHz,
                                       suffix="_sub_smoothtoCband.fits",
                                       outsuffix="_tex%i%s%s" % (tex,'_smoothtoCband',suffix),
                                       etamb=0.886)
        else:
            sdpy.makecube.make_taucube(filepath,
                                       continuum22,
                                       tex=tex,
                                       linefreq=14.488*u.GHz,
                                       suffix='_sub.fits',
                                       outsuffix="_tex%i%s" % (tex,suffix),
                                       etamb=0.886)
