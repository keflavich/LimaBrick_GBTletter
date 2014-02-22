from paths import outpath
import os
import aplpy
from aplpy_pars import set_aplpy_pars
import pylab as pl

obsmode = 'DecLatMap'
refscans = [9,54]
scanrange = [9,54]
ref1,ref2 = refscans

ifnum = 0
sampler = 'A9'
feednum = 1

mappars = zip(('DecLatMap','RALongMap','DecLatMap'),
              ([9,54],[62,98],[108,140]),
              ([9,54],[62,98],[108,140]))

fignum = 1
for min_scale_reference in (False,1,10):
    for tsysmethod in ('perint','perscan'):
        FF = [1,2,3]
        for ii,(obsmode,refscans,scanrange) in enumerate(mappars):


            experiment_suffix = ("%s_msr%i_s%i-%i" %
                                 (tsysmethod, min_scale_reference,
                                  scanrange[0], scanrange[1]))
            #fn = (outpath+'14A_110_%ito%i_%s_F%i_%s.fits' %
            #      (ref1, ref2, sampler, feednum, experiment_suffix))

            #fn = (AGBT14A_110_path+
            #      "AGBT14A_110_01_{0}_fd{1}_if{2}_sr{3}-{4}_cube_continuum.fits".
            #      format(sampler, feednum, ifnum, ref1, ref2))
            fn = os.path.join(outpath,
                              'LimaBean_H2CO22_cube_experiment_'+
                              experiment_suffix+
                              "_continuum.fits")

            if ii == 0:
                fig = pl.figure(fignum)
                fig.clf()
                fignum += 1

            FF[ii] = aplpy.FITSFigure(fn, subplot=(1,3,ii+1),
                                      convention='calabretta',
                                      figure=fig)
            FF[ii].show_colorscale(vmin=-0.1,vmax=8, vmid=1, stretch='arcsinh', cmap=pl.cm.hot)
            FF[ii].show_contour(fn, levels=[2.5], colors=['b'])

            set_aplpy_pars(FF[ii], width=0.3, height=0.3)
            if ii>0:
                FF[ii].axis_labels.hide_y()
                FF[ii].tick_labels.hide_y()
            if ii != 1:
                FF[ii].axis_labels.hide_x()
            FF[ii]._ax1.set_title('Obs %i' % (ii+1))

        #FF[1]._ax1.set_title("%s, min=%i" % (tsysmethod, min_scale_reference))
        FF[2].add_colorbar()

        figfilename = (outpath+"14A_110_scan_comparison_%s_msr%i.pdf" %
                       (tsysmethod, min_scale_reference))
        FF[ii].save(figfilename)


pl.show()
