"""
Create an H2CO 1-1 data cube with the "background" subtracted.  Uses an annular
"rind" mask to estimate the background, then interpolates through the central
regions (where The Brick is).
"""
from paths import datapath
import FITS_tools
import FITS_tools.cube_regrid
from astropy.io import fits
import pyregion
import numpy as np
import aplpy
import pylab as pl

cube22f = fits.open(datapath+'LimaBean_H2CO22_taucube.fits')
cube22 = cube22f[0].data
cube22h = cube22f[0].header
# v = -9 to + 57
summed22 = cube22[390:456,:,:].sum(axis=0)
cutout = summed22[40:60,40:60]
cmask = cutout > 1.5
mask22 = np.zeros_like(summed22,dtype='bool')
mask22[40:60,40:60] = cmask

cube11f = fits.open(datapath+'LimaBean_H2CO11_taucube.fits')
cube11h = cube11f[0].header
cube11 = cube11f[0].data

summed11 = cube11[390:456,:,:].sum(axis=0)
cutout = summed11[35:65,35:65]
cmask = cutout > 5.25
mask11 = np.zeros_like(summed11,dtype='bool')
mask11[35:65,35:65] = cmask
flathead = FITS_tools.strip_headers.flatten_header(cube11h)
maskHDU11 = fits.PrimaryHDU(data=mask11.astype('int'),header=flathead)
maskHDU11.writeto(datapath+'brick_mask11.fits',clobber=True)

pl.figure(1)
pl.clf()
mfig11 = aplpy.FITSFigure(fits.PrimaryHDU(data=summed11,header=flathead),
                          figure=pl.figure(1),
                          convention='calabretta', subplot=(1,2,1))
mfig11.show_grayscale()
mfig11.show_contour(maskHDU11,levels=[0.5])

maskHDU22 = fits.PrimaryHDU(data=mask22.astype('int'),header=flathead)
maskHDU22.writeto(datapath+'brick_mask22.fits',clobber=True)

mfig22 = aplpy.FITSFigure(fits.PrimaryHDU(data=summed22,header=flathead),
                          figure=pl.figure(1),
                          convention='calabretta', subplot=(1,2,2))
mfig22.axis_labels.hide_y()
mfig22.tick_labels.hide_y()
mfig22.show_grayscale()
mfig22.show_contour(maskHDU22,levels=[0.5])

for F in (mfig11,mfig22):
    F.set_tick_labels_xformat('dd.d')
    F.set_tick_labels_yformat('dd.d')
    F.add_colorbar()

regions = pyregion.open("./data/brick_ellipse.reg")
maskreg = regions.get_mask(maskHDU11)
# ellipse = [0.2570,0.0222,0.0426,0.0180,317]

# Mask out the region to interpolate into
cube11[:,maskreg] = np.nan

bgcube = FITS_tools.cube_regrid.spatial_smooth_cube(cube11, 10, interpolate_nan=True)

cube11f[0].data = bgcube
cube11f.writeto(datapath+'LimaBean_H2CO11_interpolated_background_cube.fits',clobber=True)

cube11f = fits.open(datapath+'LimaBean_H2CO11_taucube.fits')
cube11f[0].data -= bgcube
cube11f.writeto(datapath+'LimaBean_H2CO11_taucube_backgroundsubtracted.fits',clobber=True)
