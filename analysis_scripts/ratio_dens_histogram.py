import pyspeckit
import pylab as pl
import numpy as np

datapath = '/Users/adam/work/gc/limabean/'

r = pyspeckit.Spectrum(datapath+'Brick_RatioSpectrum_BGsub.fits')
r2 = pyspeckit.Spectrum(datapath+'Brick_RatioSpectrum.fits')

pl.figure(1)
pl.clf()
pl.hist(r2.slice(-25,55,units='km/s').data,bins=np.linspace(0,10),alpha=0.5)
pl.hist(r.slice(-25,55,units='km/s').data,bins=np.linspace(0,10),alpha=0.5)

pl.figure(2)
pl.clf()
pl.hist(r2.slice(61,91,units='km/s').data,bins=np.linspace(0,10),alpha=0.5)
pl.hist(r.slice(61,91,units='km/s').data,bins=np.linspace(0,10),alpha=0.5)
