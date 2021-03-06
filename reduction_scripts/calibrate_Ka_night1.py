import astropy.io.fits as pyfits
from sdpy import makecube,make_off_template,calibrate_map_scans
import numpy as np
from astropy import units as u
from paths import AGBT14A_110_2_path,outpath

sourcename = "LimaBean"
mapname = 'LimaBean'

filename = AGBT14A_110_2_path+'AGBT14A_110_02.raw.acs.fits'
filepyfits = pyfits.open(filename,memmap=True)
datapfits = filepyfits[1].data
dataarr = datapfits.DATA

samplers = {
        0: ["A9","A13", ],
        1: ["A10","A14",],
        2: ["B17","B21",],
        3: ["B18","B22",],
        }

feeds = {
        0: [1,2],
        1: [1,2],
        2: [1,2],
        3: [1,2]
        }

#for obsmode,refscans,scanrange in zip(('DecLatMap','RALongMap','DecLatMap'),([9,54],[62,98],[108,140]),([9,54],[62,98],[108,140])):
for obsmode,refscans,scanrange in zip(('DecLatMap','RALongMap'),([29,79],[80,96]),([30,78],[81,96])):

    ref1,ref2 = refscans

    for ifnum in samplers:
        for sampler,feednum in zip(samplers[ifnum],feeds[ifnum]):

            savefile = AGBT14A_110_2_path+"AGBT14A_110_02_{0}_fd{1}_if{2}_sr{3}-{4}".format(sampler,feednum,ifnum,ref1,ref2)

            # determine_best_off_Ku reveals that there is no need to interpolate the offs;
            # even with a standard median there is no obvious signal
            # (even at the 1% level, no signal!  This is surprising.)
            if sampler in ('A9','A13','C25','C29'):
                off_template,off_template_in = make_off_template.make_off(filename, scanrange=scanrange,
                        #exclude_velo=[-10,70], interp_vrange=[-150,250],
                        interp_polyorder=10, sampler=sampler, return_uninterp=True,
                        feednum=feednum,
                        percentile=50,
                        sourcename=sourcename,
                        savefile=savefile,
                        clobber=True,
                        #debug=True,
                        linefreq=28.9748e9, # needed to get velo right...
                        extension=1, exclude_spectral_ends=10)
            # there is some RFI or weird ringing (recvr related?) from 70-80 km/s, about
            # This section of the data should be flagged out, actually
            # (note that this section was never used: if the MEDIAN shows this, it
            # should be treated as "real off")
            #elif sampler in ('C25','C29'):
            #    off_template,off_template_in = make_off_template.make_off(filename, scanrange=scanrange,
            #            #exclude_velo=[-10,70], interp_vrange=[-150,250],
            #            interp_polyorder=10, sampler=sampler, return_uninterp=True, 
            #            feednum=feednum,
            #            percentile=50,
            #            sourcename=sourcename,
            #            savefile=savefile,
            #            clobber=True,
            #            #debug=True,
            #            exclude_velo=[72,80],
            #            interp_vrange=[50,100],
            #            linefreq=14.48848e9, # needed to get velo right...
            #            extension=1, exclude_spectral_ends=10)

            elif sampler in ["B17","B21","D33","D37"]: # h2 13CO
                off_template = make_off_template.make_off(filename, scanrange=scanrange,
                        #exclude_velo=[-10,70], interp_vrange=[-150,250],
                        interp_polyorder=10, sampler=sampler,
                        feednum=feednum,
                        sourcename=sourcename,
                        savefile=savefile,
                        clobber=True,
                        extension=1, exclude_spectral_ends=10)
            elif sampler in ["B18","B22","D34","D38"]: # h2c18o
                off_template = make_off_template.make_off(filename, scanrange=scanrange,
                        #exclude_velo=[-10,70], interp_vrange=[-150,250],
                        interp_polyorder=10, sampler=sampler,
                        feednum=feednum,
                        sourcename=sourcename,
                        savefile=savefile,
                        clobber=True,
                        extension=1, exclude_spectral_ends=10)
            else:
                off_template = None

            outfn = outpath+'14A_110_2_%ito%i_%s_F%i.fits' % (ref1,ref2,sampler,feednum)
            calibrate_map_scans.calibrate_cube_data(filename,
                                                    outfn,
                                                    scanrange=scanrange,
                                                    #min_scale_reference=10,
                                                    feednum=feednum,
                                                    refscans=refscans,
                                                    sampler=sampler,
                                                    filepyfits=filepyfits,
                                                    datapfits=datapfits,
                                                    tau=0.03,
                                                    dataarr=dataarr,
                                                    obsmode=obsmode,
                                                    sourcename=sourcename,
                                                    off_template=off_template)


