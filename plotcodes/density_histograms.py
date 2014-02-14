from astropy.io import fits
import numpy as np
import pylab as pl
import matplotlib
from astropy import units as u
from astropy import constants

matplotlib.rc_file('ggplotrc')

datapath = '/Users/adam/work/gc/limabean/'
figpath = '/Users/adam/work/h2co/limabean/figures/'

vrange = [-20,50]

densf = fits.open(datapath+'LimaBean_H2CO11_to_22_logdensity.fits',clobber=True)
mask11 = fits.getdata(datapath+'brick_mask11.fits').astype('bool')
hdr = densf[0].header
dens = densf[0].data

xarr = (np.arange(hdr['NAXIS3'])-hdr['CRPIX3']+1)*hdr['CDELT3']+hdr['CRVAL3']

xr = [np.argmin(np.abs(xarr-x)) for x in vrange]

dens_pts = dens[xr[0]:xr[1],mask11]
dens_pts = dens_pts[dens_pts==dens_pts]

dens_max = np.nanmax(dens[xr[0]:xr[1],mask11],axis=0)
dens_mean = np.nanmean(dens[xr[0]:xr[1],mask11],axis=0)

bins = np.linspace(1.9,4,100)


# Formatting and ticking junk
def my_formatter_fun(x, p):
    return "$%0.1f\\times10^{%i}$" % (10**(x % 1), np.floor(x))


pl.figure(1)
pl.clf()
ax = pl.gca()
pl.hist(10**dens_pts, bins=10**bins, alpha=0.8, histtype='stepfilled',
        normed=True, weights=10**dens_pts)
ylim = ax.get_ylim()
pl.hist(10**dens_max.ravel(), bins=10**bins[::3], alpha=0.5, histtype='stepfilled',
        normed=True, weights=10**dens_max.ravel())
ax.set_xscale('log')
ax.set_xlim(10**bins.min(),10**bins.max())
#ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(my_formatter_fun))
ax.xaxis.set_major_locator(matplotlib.ticker.LogLocator())
ax.xaxis.set_tick_params('major',size=10)
ax.xaxis.set_tick_params('minor',size=6)
ax.set_yticks([])

pl.xlabel("Density $n(H_2)$ cm$^{-3}$")
pl.ylabel("Fraction of Voxels")

pl.savefig(figpath+"voxel_density_histogram.pdf",bbox_inches='tight')

ax.set_xlim(10**bins.min(),1e5)
#ax.vlines(7.3e4,ylim[0],ylim[1],linewidth=3,alpha=0.4,color='k')
rho = (((1.3e5*u.M_sun)/(4/3.*np.pi*(2.8*u.pc)**3))/ (constants.m_p*2.8)).to(u.cm**-3).value
ax.vlines(rho,ylim[0],ylim[1],linewidth=3,alpha=0.4,color='k')
ax.set_ylim(*ylim)
pl.savefig(figpath+"voxel_density_histogram_withSteves.pdf",bbox_inches='tight')





pl.figure(2)
pl.clf()
pl.hist(10**dens_max.ravel(), bins=10**bins, alpha=0.8, histtype='stepfilled')
pl.hist(10**dens_mean.ravel(), bins=10**bins, alpha=0.8, histtype='stepfilled')
ax = pl.gca()
ax.set_xscale('log')
ax.set_xlim(10**bins.min(),10**bins.max())
ax.xaxis.set_tick_params('major',size=10)
ax.xaxis.set_tick_params('minor',size=6)
#ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(my_formatter_fun))
pl.xlabel("Density $n(H_2)$ cm$^{-3}$")
pl.ylabel("Number of Pixels")

# these are quicklooks that are now made in tau_ratio_cube
pl.figure(3)
pl.clf()
ax1 = pl.subplot(1,2,1)
im1 = ax1.imshow(np.nanmax(dens[xr[0]:xr[1],:,:],axis=0))
pl.colorbar(im1,ax=ax1)
ax2 = pl.subplot(1,2,2)
im2 = ax2.imshow(np.nanmean(dens[xr[0]:xr[1],:,:],axis=0))
pl.colorbar(im2,ax=ax2)

pl.show()
