from paths import datapath
from astropy.io import fits
from FITS_tools import strip_headers
import numpy as np

vrange = [-20,50]

def integrate(vrange=vrange,suffix='_integrated.fits'):
    for fn in ['LimaBean_H2CO11_taucube.fits',
               'LimaBean_H2CO11_taucube_claw.fits',
               'LimaBean_H2CO22_taucube.fits',]:
        fn = datapath+fn
        hdu = fits.open(fn)
        hdr = hdu[0].header

        xarr = (np.arange(hdr['NAXIS3'])-hdr['CRPIX3']+1)*hdr['CDELT3']+hdr['CRVAL3']

        xr = [np.argmin(np.abs(xarr-x)) for x in vrange]

        integ = hdu[0].data[xr[0]:xr[1],:,:].sum(axis=0)
        integ *= hdr['CDELT3']

        fhdr = strip_headers.flatten_header(hdr)

        hdu[0].data = integ
        hdu[0].header = fhdr
        hdu.writeto(fn.replace(".fits",suffix),clobber=True)

integrate()
integrate([61,91],suffix='_70kmscloud_integrated.fits')


