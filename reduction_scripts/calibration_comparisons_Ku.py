"""
Not part of the standard data pipeline:
    This is a set of experiments to examine whether there may be issues with
    the continuum estimated with the standard reduction technique.
"""
from astropy.io import fits
pyfits = fits
from sdpy import makecube,make_off_template,calibrate_map_scans
from sdpy.make_off_template import make_off
import numpy as np
from astropy import units as u
from paths import AGBT14A_110_path,outpath

sourcename = "LimaBean"
mapname = 'LimaBean'
outpath = outpath % mapname

filename = AGBT14A_110_path+'AGBT14A_110_01.raw.acs.fits'
filepyfits = pyfits.open(filename,memmap=True)
datapfits = filepyfits[1].data
dataarr = datapfits.DATA

obsmode = 'DecLatMap'
refscans = [9,54]
scanrange = [9,54]
ref1,ref2 = refscans

ifnum = 0
sampler = 'A9'
feednum = 1


for min_scale_reference in (False,1,10):
    for tsysmethod in ('perint','perscan'):

        # determine_best_off_Ku reveals that there is no need to interpolate the offs;
        # even with a standard median there is no obvious signal
        # (even at the 1% level, no signal!  This is surprising.)
        savefile = AGBT14A_110_path+"AGBT14A_110_01_{0}_fd{1}_if{2}_sr{3}-{4}".format(sampler,feednum,ifnum,ref1,ref2)
        if False: # skip this because we don't want to override the "real" version
            off_template,off_template_in = \
                    make_off(filename, scanrange=scanrange,
                             #exclude_velo=[-10,70],
                             interp_vrange=[-150,250],
                             interp_polyorder=10,
                             sampler=sampler,
                             return_uninterp=True,
                             feednum=feednum, percentile=50,
                             sourcename=sourcename,
                             savefile=savefile, clobber=True,
                             #debug=True,
                             # needed to get velo right...
                             linefreq=14.48848e9,
                             extension=1,
                             exclude_spectral_ends=10)
        else:
            off_template, off_template_in = fits.getdata(savefile+'_offspectra.fits')

        experiment_suffix = "%s_msr%i" % (tsysmethod,min_scale_reference)
        outfn = outpath+'14A_110_%ito%i_%s_F%i_%s.fits' % (ref1,ref2,sampler,feednum,experiment_suffix)
        calibrate_map_scans.calibrate_cube_data(filename,
                                                outfn,
                                                scanrange=scanrange,
                                                min_scale_reference=min_scale_reference,
                                                feednum=feednum,
                                                refscans=refscans,
                                                sampler=sampler,
                                                filepyfits=filepyfits,
                                                datapfits=datapfits,
                                                tau=0.0116,
                                                dataarr=dataarr,
                                                obsmode=obsmode,
                                                sourcename=sourcename,
                                                tsysmethod=tsysmethod,
                                                off_template=off_template)

        cubename = os.path.join(outpath,'LimaBean_H2CO22_cube_experiment_'+experiment_suffix)
        makecube.generate_header(0.256, 0.0220, naxis1=100, naxis2=100, pixsize=15,
                naxis3=800, cd3=1.0, clobber=True, restfreq=14.48848e9)
        makecube.make_blank_images(cubename,clobber=True)

        fn = outfn
        makecube.add_file_to_cube(fn,
                                  cubename+'.fits',
                                  nhits=cubename+'_nhits.fits', wcstype='V',
                                  chmod=True,
                                  add_with_kernel=True,
                                  kernel_fwhm=20./3600.,
                                  velocityrange=[-400,400],excludefitrange=[-125,250],
                                  diagnostic_plot_name=fn.replace('.fits','_data_scrubbed.png'),
                                  smoothto=2)

        os.system(cubename+'_starlink.sh')

        makecube.make_flats(cubename,vrange=[-20,60],noisevrange=[250,300])
