"""
Integrate fit gaussian
----------------------

Fit a 2D gaussian to the integrated optical depth cube of the 2-2 line to
determine the approximate volume and expected volume density of The Brick.
"""
from paths import datapath
from astropy.io import fits
from gaussfitter import gaussfitter
from astropy import units as u
from astropy import constants
import pylab as pl
import numpy as np

cube22 = fits.getdata(datapath+'LimaBean_H2CO22_taucube.fits')
cube22h = fits.getheader(datapath+'LimaBean_H2CO22_taucube.fits')

# v = -9 to + 57

summed = cube22[390:456,:,:].sum(axis=0)
cutout = summed[40:60,40:60]
fitpars,fitimg = gaussfitter.gaussfit(cutout,returnfitimage=True)

(height, amplitude, x, y, width_x, width_y, rota) = fitpars
fwhm = np.sqrt(8*np.log(2))
print "Widths: ",width_x,width_y
width_x_deg,width_y_deg = width_x*abs(cube22h['CDELT1'])*u.deg,width_y*abs(cube22h['CDELT2'])*u.deg
print "Widths: ",width_x_deg,width_y_deg
width_x_as,width_y_as = width_x_deg.to(u.arcsec),width_y_deg.to(u.arcsec)
print "Widths: ",width_x_as,width_y_as

width_x_pc = ((width_x_as/u.rad) *8500*u.pc).to(u.pc)
width_y_pc = ((width_y_as/u.rad) *8500*u.pc).to(u.pc)
print "PC: ",width_x_pc,width_y_pc
print "Effective radius (pc): ",(width_x_pc*width_y_pc)
print "Longmore 2012 effective radius: 2.8pc"
volume_prolate = 4/3. * np.pi * width_x_pc * width_y_pc**2
volume_oblate = 4/3. * np.pi * width_x_pc**2 * width_y_pc
print "Prolate: "
print "Volume: ",volume_prolate

# From Longmore et al 2012, based on dust
mass = 1.2e5*u.M_sun
massdensity = (mass/volume_prolate).to(u.g/u.cm**3)
print "Mass Density: ",massdensity
voldensity = (massdensity / (2.8 * constants.m_p)).to(u.cm**-3)
print "Volume Density: ",voldensity

print "Oblate: "
print "Volume: ",volume_oblate
massdensity = (mass/volume_oblate).to(u.g/u.cm**3)
print "Mass Density: ",massdensity
voldensity = (massdensity / (2.8 * constants.m_p)).to(u.cm**-3)
print "Volume Density: ",voldensity

pl.figure(10)
pl.clf()
pl.subplot(1,3,1)
pl.imshow(cutout)
pl.colorbar()
pl.subplot(1,3,2)
pl.imshow(fitimg)
pl.colorbar()
pl.subplot(1,3,3)
pl.imshow(cutout-fitimg)
pl.colorbar()

# filling factor stuff
shape = [525/15.,525/15.]
xc,yc = (shape[1]-1)/2.,(shape[0]-1)/2.
sourcepars = [0, 1, xc, yc, width_x_as.value/15., width_y_as.value/15., rota]
beampars = [0, 1, xc, yc, 2.6/fwhm*60/15., 2.6/fwhm*60/15., 0]
beamimg = gaussfitter.twodgaussian(beampars, shape=shape)
sourceimg = gaussfitter.twodgaussian(sourcepars, shape=shape) 
beamimg/=beamimg.max()
sourceimg/=sourceimg.sum()
print "Filling Factor: ",(beamimg*sourceimg).sum()

pl.figure(11)
pl.clf()
pl.subplot(1,3,1)
pl.imshow(beamimg)
pl.colorbar()
pl.subplot(1,3,2)
pl.imshow(sourceimg)
pl.colorbar()
pl.subplot(1,3,3)
pl.imshow(beamimg*sourceimg)
pl.colorbar()

sourcepars = [0, 1, xc-30/15., yc+30/15., width_x_as.value, width_y_as.value, rota]
sourceimgshift = gaussfitter.twodgaussian(sourcepars, shape=shape) 
sourceimgshift/=sourceimgshift.sum()
print "Filling Factor (offset by 30,30): ",(beamimg*sourceimgshift).sum()

sourcepars = [0, 1, xc-90/15., yc+90/15., width_x_as.value, width_y_as.value, rota]
sourceimgshift = gaussfitter.twodgaussian(sourcepars, shape=shape) 
sourceimgshift/=sourceimgshift.sum()
print "Filling Factor (offset by 90,90): ",(beamimg*sourceimgshift).sum()

print "The filling factor appears to not be the correction factor for the elliptical gaussian."
from agpy import smooth # to simulate what is done to the 2 cm image...
si = sourceimg.copy()
si /= si.max()
print "The correction factor is ",smooth(si,4).max()
