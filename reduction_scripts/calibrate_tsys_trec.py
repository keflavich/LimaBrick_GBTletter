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
from sdpy import calibrate_map_scans as cms
import pylab as pl
pl.matplotlib.rc_file('/Users/adam/.matplotlib/ggplotrc')

filename = AGBT14A_110_path+'AGBT14A_110_01.raw.acs.fits'
filepyfits = pyfits.open(filename,memmap=True)
datapfits = filepyfits[1].data
dataarr = datapfits.DATA
data = datapfits

exslice=slice(4096*0.1,4096*0.9)
sdpy.calibrate_map_scans.compute_tsys(data, tsysmethod='perint',
                                      exslice=exslice)

isfin = np.isfinite(data['DATA'].sum(axis=1))

for sampler in ('A9',):# 'C25','A9','A13',):
    feednum = 1 if 'A' in sampler else 2
    OK = data['SAMPLER'] == sampler
    OK *= (data['FEED'] == feednum) * isfin

    OKsource = OK.copy()
    CalOn  = (data['CAL']=='T')*isfin
    CalOff = (data['CAL']=='F')*isfin

    # only do this for the first one
    refscans=[9,54]
    temp_ref = cms.get_reference(data, refscans, CalOn=CalOn, CalOff=CalOff,
                                 exslice=exslice, OK=OK)
    LSTrefs, refarray, ref_cntstoK, tsysref = temp_ref
    OKobs1 = OK * (refscans[0] < data['SCAN'])*(data['SCAN'] < refscans[1])
    ref_scale,ref_airmass = cms.get_min_scale_reference(data,
                                                        1,
                                                        OKsource=OKobs1,
                                                        CalOn=CalOn,
                                                        CalOff=CalOff,
                                                        exslice=exslice,
                                                        airmass_method='maddalena')
    specRef = (refarray[1,:]+refarray[0,:])/2.0
    specRef = specRef/specRef[exslice].mean() * ref_scale

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

    cntstok = tcal / diffmean

    tsys = ( offmean / diffmean * tcal + tcal/2.0 )
    airmass = sdpy.calibrate_map_scans.elev_to_airmass(elev)

    OK2 = data['SAMPLER'] == sampler
    OK2 *= data['FEED'] == feednum
    OK2 *= np.isfinite(data['DATA'].sum(axis=1))

    slope, trec = np.polyfit(airmass,tsys,1)

    # same, but for just the off positions
    el1, ts1 = (data['ELEVATIO'][OK*CalOn*(data['OBJECT']=='LimaBeanOff')],
                data['TSYS'][OK*CalOn*(data['OBJECT']=='LimaBeanOff')])
    am1 = sdpy.calibrate_map_scans.elev_to_airmass(el1)

    ok = (am1==am1)*(ts1==ts1)
    slope, trec = np.polyfit(am1[ok],ts1[ok],1)

    # from the logs
    tatm = 259
    tau  = -1*np.log(1-(tsys-trec)/tatm)
    tauz = np.median(tau / airmass)
    tauz = 0.0110

    tsys_corrected = tsys - (np.exp(tauz*airmass)-1)*tatm
    ref_correction = tatm/ref_cntstoK*(np.exp(-tauz*ref_airmass)-np.exp(-tauz*airmass))

    pl.figure(1)
    pl.clf()
    pl.plot(airmass,tsys,'.', alpha=0.7, zorder=1)
    pl.plot(airmass,tsys_corrected,'.', alpha=0.5, zorder=0, color='r')
    pl.plot(am1,ts1,'.', alpha=0.9, zorder=3)
    pl.plot(np.linspace(2.5, 4.5), np.linspace(2.5, 4.5)*slope+trec,
            color='k',  alpha=0.5,  linewidth=3,  linestyle='--', zorder=5)
    pl.xlabel("Airmass")
    pl.ylabel("$T_{sys}$")
    pl.savefig(outpath+'TSYS_vs_airmass_BrickKU_%s.pdf' % sampler, bbox_inches='tight')
    pl.quiver(airmass[:-200:200], tsys[:-200:200],
              airmass[200::200]-airmass[:-200:200],
              tsys[200::200]-tsys[:-200:200],  scale_units='xy',  angles='xy',
              scale=1,  alpha=0.5,  linewidth=3,  color='r',  zorder=10,
              edgecolor='none')
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
    #pl.quiver(az[:-1000:1000], tsys[:-1000:1000],
    #          az[1000::1000]-az[:-1000:1000], tsys[1000::1000]-tsys[:-1000:1000],
    #          scale_units='xy',  angles='xy',  scale=1,  alpha=0.5,  linewidth=3,
    #          color='r',  zorder=10,  edgecolor='none')

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

    pl.figure(5)
    pl.clf()
    pl.plot(cntstok)
    pl.xlabel("Integration #")
    pl.ylabel("K/count")

    pl.figure(6)
    pl.clf()
    pl.plot(ref_correction)
    pl.xlabel("Integration #")
    pl.ylabel("Reference Level Correction")

    pl.figure(7)
    pl.clf()
    pl.plot(onmean[OKobs1[CalOn*OKsource]]  - ref_correction[OKobs1[CalOn*OKsource]])
    pl.plot(offmean[OKobs1[CalOn*OKsource]] - ref_correction[OKobs1[CalOn*OKsource]])
    pl.xlabel("Integration #")
    pl.ylabel("Instrument Counts, de-corrected for airmass")
