from paths import datapath
import FITS_tools
import FITS_tools.cube_regrid
from astropy.io import fits
import pyregion
import numpy as np

cube22 = fits.getdata(datapath+'LimaBean_H2CO22_taucube.fits')
cube22h = fits.getheader(datapath+'LimaBean_H2CO22_taucube.fits')
# v = -9 to + 57
summed = cube22[390:456,:,:].sum(axis=0)
cutout = summed[40:60,40:60]
cmask = cutout > 0.5
mask22 = np.zeros_like(summed,dtype='bool')
mask22[40:60,40:60] = cmask

cube11f = fits.open(datapath+'LimaBean_H2CO11_taucube.fits')
cube11h = cube11f[0].header
cube11 = cube11f[0].data

summed = cube11[390:456,:,:].sum(axis=0)
cutout = summed[35:65,35:65]
cmask = cutout > 3
mask11 = np.zeros_like(summed,dtype='bool')
mask11[35:65,35:65] = cmask
flathead = FITS_tools.strip_headers.flatten_header(cube11h)
maskHDU = fits.PrimaryHDU(data=mask11.astype('int'),header=flathead)
maskHDU.writeto(datapath+'brick_mask11.fits',clobber=True)

maskHDU = fits.PrimaryHDU(data=mask22.astype('int'),header=flathead)
maskHDU.writeto(datapath+'brick_mask22.fits',clobber=True)

regions = pyregion.open(datapath+"brick_ellipse.reg")
maskreg = regions.get_mask(maskHDU)
# ellipse = [0.2570,0.0222,0.0426,0.0180,317]

# Mask out the region to interpolate into
cube11[:,maskreg] = np.nan

bgcube = FITS_tools.cube_regrid.spatial_smooth_cube(cube11, 10, interpolate_nan=True)

cube11f[0].data = bgcube
cube11f.writeto(datapath+'LimaBean_H2CO11_interpolated_background_cube.fits',clobber=True)

cube11f = fits.open(datapath+'LimaBean_H2CO11_taucube.fits')
cube11f[0].data -= bgcube
cube11f.writeto(datapath+'LimaBean_H2CO11_taucube_backgroundsubtracted.fits',clobber=True)
