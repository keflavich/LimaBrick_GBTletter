"""
This is analysis code to determine the atmospheric optical depth and examine
data quality.  It is not part of the automatic data pipeline.

Analysis of this data revealed that the first and last of the 3 maps made
during this session were severely affected by excess noise at low elevations
(airmass > 3).  This changed the overall calibration, not just adding noise but
damaging signal, so it has to be excluded.
"""
from paths import AGBT14A_110_path,outpath
import astropy.io.fits as pyfits
import numpy as np
import sdpy
import pylab as pl
pl.matplotlib.rc_file('/Users/adam/.matplotlib/ggplotrc')

outpath = outpath % 'LimaBean'

filename = AGBT14A_110_path+'AGBT14A_110_01.raw.acs.fits'
filepyfits = pyfits.open(filename,memmap=True)
datapfits = filepyfits[1].data
dataarr = datapfits.DATA
data = datapfits

exslice=slice(4096*0.1,4096*0.9)
sdpy.calibrate_map_scans.compute_tsys(data, tsysmethod='perint',
                                      exslice=exslice)

isfin = np.isfinite(data['DATA'].sum(axis=1))

for sampler in ('C25','A9','A13','C29',):
    feednum = 1 if 'A' in sampler else 2
    OK = data['SAMPLER'] == sampler
    OK *= (data['FEED'] == feednum) * isfin

    OKsource = OK.copy()
    CalOn  = (data['CAL']=='T')*isfin
    CalOff = (data['CAL']=='F')*isfin

    on_data = dataarr[CalOn*OKsource,exslice]
    off_data = dataarr[CalOff*OKsource,exslice]
    tcal = data['TCAL'][CalOn*OKsource]
    elev = data['ELEVATIO'][CalOn*OKsource]
    az = data['AZIMUTH'][CalOn*OKsource]

    # terrible, disgusting, absolutely terrifying hack
    if sampler == 'C29':
        sz = min([on_data.shape[0],off_data.shape[0]])
        on_data = on_data[:sz,:]
        off_data = off_data[:sz,:]
        tcal = tcal[:sz]
        elev = elev[:sz]
        az = az[:sz]

    offmean = np.mean(off_data,axis=1)
    onmean  = np.mean(on_data,axis=1)
    diffmean = onmean-offmean

    tsys = ( offmean / diffmean * tcal + tcal/2.0 )
    airmass = 1/np.cos((90-elev)/180*np.pi)

    OK2 = data['SAMPLER'] == sampler
    OK2 *= data['FEED'] == feednum
    OK2 *= np.isfinite(data['DATA'].sum(axis=1))

    slope, trec = np.polyfit(airmass,tsys,1)

    # same, but for just the off positions
    el1,ts1 = (data['ELEVATIO'][OK*CalOn*(data['OBJECT']=='LimaBeanOff')],data['TSYS'][OK*CalOn*(data['OBJECT']=='LimaBeanOff')])
    am1 = 1/np.cos((90-el1)/180*np.pi)

    slope, trec = np.polyfit(am1,ts1,1)

    # from the logs
    tatm = 259
    tau  = -1*np.log(1-(tsys-trec)/tatm)
    tauz = tau / airmass

    pl.figure(1)
    pl.clf()
    pl.plot(airmass,tsys,'.', alpha=0.7, zorder=1)
    pl.plot(am1,ts1,'.', alpha=0.9, zorder=3)
    pl.plot(np.linspace(2.5,4.5),np.linspace(2.5,4.5)*slope+trec, color='k', alpha=0.5, linewidth=3, linestyle='--',zorder=5)
    pl.xlabel("Airmass")
    pl.ylabel("$T_{sys}$")
    pl.savefig(outpath+'TSYS_vs_airmass_BrickKU_%s.pdf' % sampler, bbox_inches='tight')
    pl.quiver(airmass[:-200:200],tsys[:-200:200],airmass[200::200]-airmass[:-200:200],tsys[200::200]-tsys[:-200:200], scale_units='xy', angles='xy', scale=1, alpha=0.5, linewidth=3, color='r', zorder=10, edgecolor='none')
    pl.savefig(outpath+'TSYS_vs_airmass_BrickKU_arrows_%s.pdf' % sampler, bbox_inches='tight')

    pl.figure(2)
    pl.clf()
    az1 = data['AZIMUTH'][OK*CalOn*(data['OBJECT']=='LimaBeanOff')]
    pl.plot(az,tsys,'.')
    pl.plot(az1,ts1,'.')
    pl.xlabel('Azimuth')
    pl.ylabel('$T_{sys}$')
    pl.savefig(outpath+"TSYS_vs_azimuth_BrickKU_%s.pdf" % sampler,bbox_inches='tight')
    #plot(az,tsys/(airmass/2.5),'.')
    #pl.quiver(az[:-1000:1000],tsys[:-1000:1000],az[1000::1000]-az[:-1000:1000],tsys[1000::1000]-tsys[:-1000:1000], scale_units='xy', angles='xy', scale=1, alpha=0.5, linewidth=3, color='r', zorder=10, edgecolor='none')

    pl.figure(3)
    pl.clf()
    pl.plot(onmean,label='On')
    pl.plot(offmean,label='Off')
    pl.plot(diffmean+1.5,label='Diff+1.5')
    pl.xlabel("Integration #")
    pl.ylabel("Instrument Counts")
    pl.legend(loc='best')
    pl.savefig(outpath+"OnOffvsTime_BrickKU_%s.pdf" % sampler,bbox_inches='tight')

    pl.figure(4)
    pl.clf()
    pl.plot(data['SCAN'][CalOn*OKsource])
    pl.xlabel("Integration #")
    pl.ylabel("Scan #")
