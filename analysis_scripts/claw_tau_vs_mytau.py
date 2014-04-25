from paths import datapath, figpath
from astropy.io import fits
import numpy as np
from mpl_plot_templates import adaptive_param_plot
import pylab as pl

h2co11_my = fits.getdata(datapath+'LimaBean_H2CO11_taucube_tex1.5.fits')
h2co11_law = fits.getdata(datapath+'LimaBean_H2CO11_taucube_claw_tex1.5.fits')

noise11 = h2co11_my[:50,:,:].std(axis=0)

sn11 = h2co11_my/noise11

mask = (sn11 > 2)

brick_mask = fits.getdata(datapath+'brick_mask11.fits').astype('bool')
mask_nan = np.ones(mask.shape,dtype='float')
mask_nan[np.logical_not(mask)] = np.nan


fig = pl.figure(0)
fig.clf()
ax = pl.gca()
pl.plot([0,0.7],[0,0.7],'k--',linewidth=3,alpha=0.5,zorder=0)
pl.plot([0,0.7],[0,0.7*1.65],'g--',linewidth=3,alpha=0.5,zorder=0)
adaptive_param_plot(h2co11_law[mask], h2co11_my[mask], bins=30, threshold=100,
                    marker='.', alpha=0.5, zorder=5,
                    markersize=2)
pl.plot((h2co11_law*mask_nan)[:,brick_mask],
        (h2co11_my*mask_nan)[:,brick_mask], marker='.', alpha=0.5, color='r',
        markersize=2, linestyle='none')
ax.set_xlabel(r"$\tau_{1-1}$ Law")
ax.set_ylabel(r"$\tau_{1-1}$ Ours")
ax.set_xlim(0,h2co11_law[mask].max())
ax.set_ylim(0,h2co11_my[mask].max())
pl.show()

pl.savefig(figpath+'claw_tau_vs_mytau.pdf',bbox_inches='tight')
