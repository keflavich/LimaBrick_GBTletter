from paths import datapath,figpath
from astropy.io import fits
import FITS_tools.strip_headers
import pyspeckit
import scipy.ndimage
import pylab as pl
import aplpy

mask11 = fits.getdata(datapath+'brick_mask11.fits').astype('bool')
mask22 = fits.getdata(datapath+'brick_mask22.fits').astype('bool')
cube22 = fits.getdata(datapath+'LimaBean_H2CO22_taucube.fits')
cube22h = fits.getheader(datapath+'LimaBean_H2CO22_taucube.fits')
cube11 = fits.getdata(datapath+'LimaBean_H2CO11_taucube.fits')
cube11h = fits.getheader(datapath+'LimaBean_H2CO11_taucube.fits')

aperture22 = (scipy.ndimage.morphology.binary_dilation(mask22, iterations=6)-
              scipy.ndimage.morphology.binary_dilation(mask22, iterations=2))

spechead = FITS_tools.strip_headers.speccen_header(cube22h)
fgspec = cube22[:,mask22].mean(axis=1)
bgspec = cube22[:,aperture22].mean(axis=1)
specHDU22bg = fits.PrimaryHDU(data=bgspec,header=spechead)
specHDU22bg.writeto(datapath+"BackgroundH2CO_22.fits",clobber=True)
specHDU22 = fits.PrimaryHDU(data=fgspec,header=spechead)
specHDU22.writeto(datapath+"BrickSpectrumH2CO_22.fits",clobber=True)

cube11f = fits.open(datapath+'LimaBean_H2CO11_taucube.fits')
cube11h = cube11f[0].header
cube11 = cube11f[0].data

spechead = FITS_tools.strip_headers.speccen_header(cube11h)
fgspec = cube11[:,mask22].mean(axis=1)
bgspec = cube11[:,aperture22].mean(axis=1)
specHDU11bg = fits.PrimaryHDU(data=bgspec,header=spechead)
specHDU11bg.writeto(datapath+"BackgroundH2CO_11.fits",clobber=True)
specHDU11 = fits.PrimaryHDU(data=fgspec,header=spechead)
specHDU11.writeto(datapath+"BrickSpectrumH2CO_11.fits",clobber=True)

sp11m2 = pyspeckit.spectrum.classes.Spectrum.from_hdu(specHDU11)
bg11m2 = pyspeckit.spectrum.classes.Spectrum.from_hdu(specHDU11bg)
sp22m2 = pyspeckit.spectrum.classes.Spectrum.from_hdu(specHDU22)
bg22m2 = pyspeckit.spectrum.classes.Spectrum.from_hdu(specHDU22bg)
sp11m2.specname='H$_2$CO 1-1 Sharp Mask'
sp11m2.plotter(figure=pl.figure(12))
bg11m2.plotter(axis=sp11m2.plotter.axis, clear=False, color='b')
d11m2 = (sp11m2-bg11m2)
d11m2.plotter(axis=sp11m2.plotter.axis, clear=False, color='r',ymin=bg11m2.data.min(),ymax=sp11m2.data.max())
sp22m2.specname='H$_2$CO 2-2 Sharp Mask'
sp22m2.plotter(figure=pl.figure(13))
bg22m2.plotter(axis=sp22m2.plotter.axis, clear=False, color='b')
d22m2 = (sp22m2-bg22m2)
d22m2.plotter(axis=sp22m2.plotter.axis, clear=False, color='r',ymin=bg22m2.data.min())

sp11m2.plotter.savefig(figpath+'avgspectrum_oneone_bgsubtraction.pdf')
sp22m2.plotter.savefig(figpath+'avgspectrum_twotwo_bgsubtraction.pdf')

aperture11 = (scipy.ndimage.morphology.binary_dilation(mask11, iterations=12)-
              scipy.ndimage.morphology.binary_dilation(mask11, iterations=4))

mask11f = fits.open(datapath+'brick_mask11.fits')
mask11f[0].data = aperture11.astype('int')
mask11f.writeto(datapath+'brick_aperture11.fits',clobber=True)

# Shrink mask11 to better approximate the peak
mask11 = scipy.ndimage.morphology.binary_erosion(mask11, iterations=5)

# since we're using C-band aperture...
cube22 = fits.getdata(datapath+'LimaBean_H2CO22_taucube_smoothtoCband.fits')
cube22h = fits.getheader(datapath+'LimaBean_H2CO22_taucube_smoothtoCband.fits')

spechead = FITS_tools.strip_headers.speccen_header(cube22h)
fgspec = cube22[:,mask11].mean(axis=1)
bgspec = cube22[:,aperture11].mean(axis=1)
specHDU22bg = fits.PrimaryHDU(data=bgspec,header=spechead)
specHDU22bg.writeto(datapath+"BackgroundH2CO_22_mask11.fits",clobber=True)
specHDU22 = fits.PrimaryHDU(data=fgspec,header=spechead)
specHDU22.writeto(datapath+"BrickSpectrumH2CO_22_mask11.fits",clobber=True)

spechead = FITS_tools.strip_headers.speccen_header(cube11h)
fgspec = cube11[:,mask11].mean(axis=1)
bgspec = cube11[:,aperture11].mean(axis=1)
specHDU11bg = fits.PrimaryHDU(data=bgspec,header=spechead)
specHDU11bg.writeto(datapath+"BackgroundH2CO_11_mask11.fits",clobber=True)
specHDU11 = fits.PrimaryHDU(data=fgspec,header=spechead)
specHDU11.writeto(datapath+"BrickSpectrumH2CO_11_mask11.fits",clobber=True)

sp11m1 = pyspeckit.spectrum.classes.Spectrum.from_hdu(specHDU11)
bg11m1 = pyspeckit.spectrum.classes.Spectrum.from_hdu(specHDU11bg)
sp22m1 = pyspeckit.spectrum.classes.Spectrum.from_hdu(specHDU22)
bg22m1 = pyspeckit.spectrum.classes.Spectrum.from_hdu(specHDU22bg)
sp11m1.specname='H$_2$CO 1-1 Smooth Mask'
sp11m1.unit=r'$\tau$'
sp11m1.plotter(figure=pl.figure(14))
bg11m1.plotter(axis=sp11m1.plotter.axis, clear=False, color='b')
d11m1 = (sp11m1-bg11m1)
d11m1.plotter(axis=sp11m1.plotter.axis, clear=False, color='r',ymin=bg11m1.data.min(), ymax=sp11m1.data.max())
sp22m1.specname='H$_2$CO 2-2 Smooth Mask'
sp22m1.unit=r'$\tau$'
sp22m1.plotter(figure=pl.figure(15))
bg22m1.plotter(axis=sp22m1.plotter.axis, clear=False, color='b')
d22m1 = (sp22m1-bg22m1)
d22m1.plotter(axis=sp22m1.plotter.axis, clear=False, color='r',ymin=bg22m1.data.min())

sp11m1.plotter.savefig(figpath+'avgspectrum_oneone_bgsubtraction_smooth.pdf')
sp22m1.plotter.savefig(figpath+'avgspectrum_twotwo_bgsubtraction_smooth.pdf')

d22m1.plotter(figure=pl.figure(16), color='r')
d11m1.plotter(axis=d22m1.plotter.axis, clear=False)
d22m1.plotter.axis.set_title("H2CO 1-1 (black) and H2CO 2-2 (red) background subtracted")

d22m1.write(datapath+'Brick_H2CO22_smooth_BGsub.fits')
d11m1.write(datapath+'Brick_H2CO11_smooth_BGsub.fits')

rbg = d22m1.copy()
rbg.data = d11m1.data/d22m1.data
rbg.plotter(figure=pl.figure(17))
rbg.units='Ratio'

r2 = d22m1.copy()
r2.data = sp11m1.data/sp22m1.data
r2.plotter(xmin=-100,xmax=150,ymin=0,ymax=13,color='r',axis=rbg.plotter.axis,clear=False)
r2.units='Ratio'
rbg.plotter.axis.set_ylabel("Ratio")

rbg.write(datapath+'Brick_RatioSpectrum_BGsub.fits')
r2.write(datapath+'Brick_RatioSpectrum.fits')

rbg.plotter.savefig(figpath+'avgspectrum_ratio_bgsubtraction_smooth.pdf')

r3 = bg11m1.copy()
r3.data = bg11m1.data/bg22m1.data
r3.specname = "Background H$_2$CO 1-1/2-2"
r3.plotter(xmin=-100,xmax=150,ymin=0,ymax=13,clear=False,figure=pl.figure(19))
r3.units='Ratio'
r3.plotter.axis.set_ylabel("Ratio")
r3.plotter.savefig(figpath+'background_avgspectrum_ratio_smooth.pdf')

dd = d11m1.copy()
dd.specname = 'H$_2$CO 1-1 minus 2-2'
dd.data = d11m1.data - d22m1.data
dd.unit=r'$\tau_{1-1}-\tau_{2-2}$'
dd.plotter(figure=pl.figure(18))
d11m1.plotter(axis=dd.plotter.axis,color='b',clear=False)
dd.plotter.savefig(figpath+'avgspectrum_11minus22.pdf')
