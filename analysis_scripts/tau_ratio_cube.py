from paths import datapath
from astropy.io import fits
import numpy as np
import h2co_modeling

h2co11 = fits.getdata(datapath+'LimaBean_H2CO11_taucube.fits')
h2co22 = fits.getdata(datapath+'LimaBean_H2CO22_taucube_smoothtoCband.fits')

noise11 = h2co11[:50,:,:].std(axis=0)
noise22 = h2co22[:50,:,:].std(axis=0)

sn11 = h2co11/noise11
sn22 = h2co22/noise22

mask = (sn11 > 2) & (sn22 > 2)
ratio = h2co11/h2co22
ratio[True-mask] = np.nan

ratioF = fits.open(datapath+'LimaBean_H2CO11_taucube.fits')
ratioF[0].data = ratio
ratioF[0].header['BUNIT'] = 'none'
ratioF.writeto(datapath+'LimaBean_H2CO11_to_22_tau_ratio.fits',clobber=True)

from astropy.convolution import convolve

filt = np.ones([3,3,3],dtype='bool')
filt[1,1,1] = 0
nneighbors = convolve(np.isfinite(ratio), filt)
ratio[(nneighbors<7) + (True-np.isfinite(nneighbors))] = np.nan

ratioF[0].data = ratio
ratioF.writeto(datapath+'LimaBean_H2CO11_to_22_tau_ratio_neighbors.fits',clobber=True)

nfin = np.isfinite(ratio).sum()
nok = (np.isfinite(sn11)*np.isfinite(sn22)).sum()
print "Unfiltered voxels: ",nfin," of ", nok," or ",nfin/float(nok)*100,"%"
nfinpix = np.isfinite(ratio).max(axis=0).sum()
nokpix = (np.isfinite(sn11)*np.isfinite(sn22)).max(axis=0).sum()
print "Unfiltered pixels ",nfinpix," of ", nokpix," or ",nfinpix/float(nokpix)*100,"%"

radexdatapath = './radex/'
faur = h2co_modeling.SmoothtauModels(datafile=radexdatapath+'faure/1-1_2-2_XH2CO_fixed_faure.dat')

abund = -8.5
tau1,tau2,dens,col = faur.select_data(abundance=abund, opr=0.1, temperature=50)
tau,vtau,vtau_ratio = faur.generate_tau_functions(abundance=abund, opr=0.1, temperature=50)

tauratio = vtau_ratio(dens, line1=tau1, line2=tau2, sigma=1.0)

ok = np.arange(tauratio.size) > np.argmax(tauratio)

def ratio_to_dens(ratio):
    inds = np.argsort(tauratio[ok])
    return np.interp(ratio, tauratio[ok][inds], dens[ok][inds], np.nan, np.nan)

dcube = ratio_to_dens(ratio)

ratioF = fits.open(datapath+'LimaBean_H2CO11_taucube.fits')
ratioF[0].data = dcube
ratioF[0].header['BUNIT'] = 'log volume density'
ratioF.writeto(datapath+'LimaBean_H2CO11_to_22_logdensity.fits',clobber=True)


hdr = ratioF[0].header
xarr = (np.arange(hdr['NAXIS3'])-hdr['CRPIX3']+1)*hdr['CDELT3']+hdr['CRVAL3']

flatfile = fits.open(datapath+'LimaBean_H2CO11_cube_continuum.fits')
flatfile[0].header['BMAJ'] = 2.5/60./np.sqrt(8*np.log(2))
flatfile[0].header['BMIN'] = 2.5/60./np.sqrt(8*np.log(2))

for vrange in ([-20,50],[50,70]):
    xr = [np.argmin(np.abs(xarr-x)) for x in vrange]
    dens_peak = np.nanmax(dcube[xr[0]:xr[1],:,:],axis=0)
    dens_mean = np.nanmean(dcube[xr[0]:xr[1],:,:],axis=0)
    flatfile[0].data = dens_peak
    flatfile.writeto(datapath+'LimaBean_H2CO11_peak_density_vr_%ito%i.fits' % tuple(vrange),clobber=True)
    flatfile[0].data = dens_mean
    flatfile.writeto(datapath+'LimaBean_H2CO11_mean_density_vr_%ito%i.fits' % tuple(vrange),clobber=True)
