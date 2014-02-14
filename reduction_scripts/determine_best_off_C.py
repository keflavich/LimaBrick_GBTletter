import astropy.io.fits as pyfits
import numpy as np
from gbtpy import makecube,make_off_template,calibrate_map_scans
import pylab as pl

scanrange=[6,21]
ref1 = 6
ref2 = 21
refscans = [6,21]#,22,32]
sourcename = "LimaBean"
obsmode = "DecLatMap"
mapname = 'LimaBean'
outpath = '/Users/adam/observations/gbt/%smap/' % mapname

filename = '/Users/adam/observations/gbt/AGBT12B_221_01/AGBT12B_221_01.raw.acs.fits'
filepyfits = pyfits.open(filename,memmap=True)
datapfits = filepyfits[1].data
dataarr = datapfits.DATA

from gbtpy.makecube import velo_iterator
velo = velo_iterator(datapfits,linefreq=None).next()

sampler = 'A9'
off_templates = {}
for percentile in (50,95,99):
    for ii,exclude in enumerate(([-58,103,158,193],[-127,103,158,193],[-127,43,69,103,158,193],[-127,-96,-78,38,62,96,158,193])):

        off_template,off_template_in,off_poly = make_off_template.make_off(filepyfits, scanrange=scanrange,
                exclude_velo=exclude, interp_vrange=[-350,350],
                interp_polyorder=10, sampler=sampler, return_uninterp=True, 
                return_poly=True,
                percentile=percentile, exclude_spectral_ends=10)

        off_templates[(percentile,ii)] = off_template,off_template_in,off_poly

oti_med = off_templates[(50,0)][1]

pl.figure(1)
pl.clf()
for ot in off_templates:
    off_template,off_template_in,off_poly = off_templates[ot]
    ax1 = pl.subplot(4,2,ot[1]+1)
    L, = ax1.plot(velo, off_template)
    ax1.plot(velo, off_template_in, label=ot, color=L.get_color())
    ax1.set_ylim(0.89,0.98)
    ax1.set_xlim(-250,250)
    pl.legend(loc='best')

    ax2 = pl.subplot(4,2,ot[1]+5)
    ax2.plot(velo, oti_med/off_template, color=L.get_color())
    ax2.plot(velo, oti_med/off_template_in, color=L.get_color(), label=ot)

    ax2.set_ylim(0.965,1.0)
    ax2.set_xlim(-250,250)
    #pl.legend(loc='best')

for percentile in (1,5,25,75,99.5,99.9):
    off_template,off_template_in,off_poly = make_off_template.make_off(filepyfits, scanrange=scanrange,
            exclude_velo=exclude, interp_vrange=[-350,350],
            interp_polyorder=10, sampler=sampler, return_uninterp=True, 
            return_poly=True,
            percentile=percentile, exclude_spectral_ends=10)

    off_templates[(percentile,ii)] = off_template,off_template_in,off_poly

ot = (99,3)

pl.figure(2)
pl.clf()

off_template,off_template_in,off_poly = off_templates[ot]
ax1 = pl.gca()
ax1.plot(velo, off_template_in, label="%s%%" % ot[0], linewidth=3, alpha=0.5)
L, = ax1.plot(velo, off_template, label="Interpolated", linewidth=3)
for percentile in (1,5,25,50,75,95,99,99.5,99.9):
    ax1.plot(velo,off_templates[(percentile,ii)][1], label="%s%%" % (percentile), alpha=0.5)
ax1.set_ylim(0.85,0.96)
ax1.set_xlim(-250,250)
pl.legend(loc='best')

pl.savefig('/Users/adam/work/h2co/limabean/figures/off_spectra_Cband.pdf')

pl.figure(3)
pl.clf()

ax1 = pl.gca()
ax1.hlines(1,-350,350,color='k',linestyle='--',alpha=0.5)
ax1.plot(np.array(zip(*[iter([-127,-96,-78,38,62,96,158,193])]*2)).T,[1.003,1.003], color='k')
for percentile in (50,75,95,99,99.5):
    L, = ax1.plot(velo,off_templates[(percentile,ii)][1]/off_templates[(percentile,ii)][2], label="%s%%" % (percentile), alpha=0.5)
    ax1.plot(velo,off_templates[(percentile,ii)][2], color=L.get_color(), alpha=0.5)
ax1.set_ylim(0.99,1.005)
ax1.set_xlim(-250,250)
pl.legend(loc='lower right')
