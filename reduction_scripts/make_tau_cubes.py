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

datapath = outpath

fivecubes = ['LimaBean_H213CO_cube_sub.fits',
             'LimaBean_H2CO11_cube_sub.fits',
             ]

# these are not used...
fifteencubes = ['LimaBean_H2CO22_cube_sub.fits',
                'LimaBean_H2CO22_cube_sub_smoothtoCband.fits',
                'LimaBean_H213CO22_cube_sub.fits',
                'LimaBean_H213CO22_cube_sub_smoothtoCband.fits',
                'LimaBean_H2C18O22_cube_sub.fits',
                'LimaBean_H2C18O22_cube_sub_smoothtoCband.fits',
                ]

# TODO: Where can these be got?
ccontfile = datapath+'GCCBand_lb.2.fits'
if os.path.exists(ccontfile):
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

    for fn in fivecubes:
        sdpy.makecube.make_taucube(datapath+fn.strip("_sub.fits"), 
                                   reproj_contfile,
                                   outsuffix="_claw.fits",
                                   etamb=0.98)
else:
    print "Casey Law's C-band continuum images are not available.",
    print "We'll use our own continuum instead (which is OK, but less reliable)."


    for fn in fivecubes:
        sdpy.makecube.make_taucube(datapath+fn.replace("_sub.fits",""),
                                   datapath+fn.replace("_sub.fits","_continuum.fits"),
                                   etamb=0.98)
