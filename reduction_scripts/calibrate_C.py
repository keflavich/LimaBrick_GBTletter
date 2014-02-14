import astropy.io.fits as pyfits
import numpy as np
from sdpy import makecube,make_off_template,calibrate_map_scans
from paths import AGBT12B_221_path,outpath

scanrange=[6,21]
ref1 = 6
ref2 = 21
refscans = [6,21]#,22,32]
sourcename = "LimaBean"
obsmode = "DecLatMap"
mapname = 'LimaBean'
outpath = outpath % mapname

obsmode = 'RALongMap'
scanrange=[22,32]
ref1 = 22
ref2 = 32
refscans = [22,32]#,22,32]

filename = AGBT12B_221_path+'AGBT12B_221_01.raw.acs.fits'
filepyfits = pyfits.open(filename,memmap=True)
datapfits = filepyfits[1].data
dataarr = datapfits.DATA

for scanrange, obsmode in zip(([6,21],[22,32]), ['DecLatMap','RALongMap']):
    for sampler in ['A10', 'A13', 'A14', 'A9', 'B17', 'B18', 'B21', 'B22', 'C25', 'C26',
                    'C29', 'C30', 'D33', 'D34', 'D37', 'D38']:

        ref1,ref2 = refscans = scanrange

        savefile = AGBT12B_221_path+'AGBT12B_221_01_{0}_sr{1}-{2}'.format(sampler,ref1,ref2)

        if sampler in ('A9','A13'):
            off_template,off_template_in = make_off_template.make_off(filepyfits,
                                                                      scanrange=scanrange,
                                                                      obsmode=obsmode,
                                                                      #exclude_velo=[-58,103,158,193],
                                                                      #exclude_velo=[-127,-96,-78,38,62,96,158,193],
                                                                      exclude_velo=[-59,-42,-33,-25,-12,8,62,96,158,193],
                                                                      interp_vrange=[-300,300],
                                                                      interp_polyorder=10,
                                                                      sampler=sampler,
                                                                      return_uninterp=True,
                                                                      percentile=99,
                                                                      savefile=savefile,
                                                                      debug=True,
                                                                      clobber=True,
                                                                      exclude_spectral_ends=10)
                    # previously used median with no exclusion; I think this is much better.
                    # 95% preserves the overall shape very well, but does a far superior job of excluding
                    # lines
                    # -127 to -58 clearly shows some structure at the 50 and 95 level, but maybe not 99
        elif sampler in ('C25','C29'): # h2 13CO
            off_template = make_off_template.make_off(filename,
                                                      scanrange=scanrange,
                                                      obsmode=obsmode,
                                                      exclude_velo=[-10,70],
                                                      interp_vrange=[-100,150],
                                                      interp_polyorder=10,
                                                      sampler=sampler,
                                                      savefile=savefile,
                                                      clobber=True,
                                                      exclude_spectral_ends=10)
        elif sampler in ('D33','D37'): # H110a
            off_template = make_off_template.make_off(filename,
                                                      scanrange=scanrange,
                                                      obsmode=obsmode,
                                                      exclude_velo=[-50,150],
                                                      interp_vrange=[-150,250],
                                                      interp_polyorder=10,
                                                      sampler=sampler,
                                                      savefile=savefile,
                                                      clobber=True,
                                                      exclude_spectral_ends=10)
        else:
            off_template = None

        print "Sanity check.  Off_template: ",off_template

        calibrate_map_scans.calibrate_cube_data(filename,
                outpath+'12B_221_%ito%i_%s_F1.fits' %
                (ref1,ref2,sampler),scanrange=scanrange,refscan1=ref1,refscan2=ref2,
                feednum=1, refscans=refscans, sampler=sampler, filepyfits=filepyfits,
                datapfits=datapfits, tau=0, dataarr=dataarr, obsmode=obsmode,
                sourcename=sourcename, off_template=off_template)
