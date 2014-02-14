import astropy.io.fits as pyfits
fits = pyfits
from gbtpy import makecube,make_off_template,calibrate_map_scans
import numpy as np
import pylab as pl
import os

sourcename = "LimaBean"
mapname = 'LimaBean'
outpath = '/Users/adam/observations/gbt/%smap/' % mapname

sampler = 'A9'

if os.path.exists('/Volumes/128gbdisk/gbt/AGBT14A_110_01.samplera9.fits'):
    filename = '/Volumes/128gbdisk/gbt/AGBT14A_110_01.samplera9.fits'
    filepyfits = pyfits.open(filename,memmap=True)
    hdul = filepyfits
else:
    filename = '/Users/adam/observations/gbt/AGBT14A_110_01/AGBT14A_110_01.raw.acs.fits'
    filepyfits = pyfits.open(filename,memmap=True)
    datapfits = filepyfits[1].data
    data = datapfits[datapfits['SAMPLER'] == sampler]
    hdu = fits.BinTableHDU(data)
    phdu = fits.PrimaryHDU()
    hdul = fits.HDUList([phdu,hdu])

datapfits = filepyfits[1].data
dataarr = datapfits['DATA']

samplers = {
        0: ["A9","A13","C25","C29"],
        1: ["A10","A14","C26","C30"],
        2: ["B17","B21","D33","D37"],
        3: ["B18","B22","D34","D38"],
        }

feeds = {
        0: [1,1,2,2],
        1: [1,1,2,2],
        2: [1,1,2,2],
        3: [1,1,2,2]
        }

from gbtpy.makecube import velo_iterator
velo = velo_iterator(datapfits,linefreq=None).next()

interp_polyorder=15

scanrange = [9,54]
off_templates = {}
for percentile in (50,95,99):
    for ii,exclude in enumerate(([-58,103,158,193],[-10,70],[])):

        off_template,off_template_in,off_poly = make_off_template.make_off(hdul, scanrange=scanrange,
                exclude_velo=exclude, interp_vrange=[-150,250],
                interp_polyorder=interp_polyorder, sampler=sampler, return_uninterp=True, 
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
    ax1.set_ylim(0.915,1.0)
    ax1.set_xlim(-150,150)
    pl.legend(loc='best')

    ax2 = pl.subplot(4,2,ot[1]+5)
    ax2.plot(velo, oti_med/off_template, color=L.get_color())
    ax2.plot(velo, oti_med/off_template_in, color=L.get_color(), label=ot)

    ax2.set_ylim(0.965,1.01)
    ax2.set_xlim(-150,150)
    #pl.legend(loc='best')

for percentile in (1,5,25,75):
    off_template,off_template_in,off_poly = make_off_template.make_off(hdul, scanrange=scanrange,
            exclude_velo=exclude, interp_vrange=[-150,250],
            interp_polyorder=interp_polyorder, sampler=sampler, return_uninterp=True, 
            return_poly=True,
            percentile=percentile, exclude_spectral_ends=10)

    off_templates[(percentile,ii)] = off_template,off_template_in,off_poly

ot = (50,2)

pl.figure(2)
pl.clf()

off_template,off_template_in,off_poly = off_templates[ot]
ax1 = pl.gca()
ax1.plot(velo, off_template_in, label="%s%%" % ot[0], linewidth=3, alpha=0.5)
L, = ax1.plot(velo, off_template, label="Interpolated", linewidth=3)
for percentile in (1,5,25,50,75,95,99):
    ax1.plot(velo,off_templates[(percentile,ii)][1], label="%s%%" % (percentile), alpha=0.5)
ax1.set_ylim(0.915,1.0)
ax1.set_xlim(-150,150)
pl.legend(loc='best')

TrueOffs = {}
for scanN in (8,61,99,101,107,141):
    TrueOffs[scanN] = make_off_template.make_off(hdul, scanrange=[scanN,scanN],
                                                 interp_vrange=[-150,250],
                                                 sampler=sampler, return_uninterp=False, 
                                                 interp_polyorder=interp_polyorder,
                                                 percentile=50, exclude_spectral_ends=10)

meantrueoff = np.mean(TrueOffs.values(),axis=0)/np.mean(TrueOffs.values())*off_templates[ot][0].mean()
pl.plot(velo, meantrueoff, 'k', alpha=0.5, linewidth=2)

pl.figure(3)
pl.clf()
pl.plot(velo, np.array(TrueOffs.values()).T, alpha=0.5, linewidth=2)
ax = pl.gca()
ax.set_ylim(0.915,1.0)
ax.set_xlim(-150,150)


pl.figure(4)
pl.clf()

ax1 = pl.gca()
ax1.hlines(1,-350,350,color='k',linestyle='--',alpha=0.5)
#ax1.plot(np.array(zip(*[iter([-127,-96,-78,38,62,96,158,193])]*2)).T,[1.003,1.003], color='k')
for percentile in (25,50,75,95,99):
    L, = ax1.plot(velo,off_templates[(percentile,ii)][1]/off_templates[(percentile,ii)][2], label="%s%%" % (percentile), alpha=0.5)
    #ax1.plot(velo,off_templates[(percentile,ii)][2], color=L.get_color(), alpha=0.5)
ax1.set_ylim(0.998,1.002)
ax1.set_xlim(-150,200)
pl.legend(loc='lower right')

pl.figure(5)
pl.clf()
ax1 = pl.gca()
for percentile in (25,50,75,95,99):
    M = off_templates[(percentile,ii)][1]/meantrueoff
    L, = ax1.plot(velo, M/M[200:-200].mean(), label="%s%%" % (percentile), alpha=0.5)
ax1.set_ylim(0.99,1.005)
ax1.set_xlim(-150,200)
pl.legend(loc='lower right')


pl.figure(6)
pl.clf()

off_template,off_template_in,off_poly = off_templates[ot]
ax1 = pl.gca()
#ax1.plot(velo, off_template_in, label="%s%%" % ot[0], linewidth=3, alpha=0.5)
#L, = ax1.plot(velo, off_template, label="Interpolated", linewidth=3)
for percentile in (1,5,25,50,75,95,99):
    M = off_templates[(percentile,ii)][1]
    L, = ax1.plot(velo, M/M[400:-400].mean(), label="%s%%" % (percentile), alpha=0.5)
    I = off_templates[(percentile,ii)][2]
    ax1.plot(velo,0.002+I/M[400:-400].mean(), color=L.get_color(), alpha=0.5)
ax1.set_ylim(0.915,1.0)
ax1.set_xlim(-150,150)
pl.legend(loc='best')
