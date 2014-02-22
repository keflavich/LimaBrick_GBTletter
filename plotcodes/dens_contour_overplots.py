from paths import datapath,figpath
import aplpy
import pylab as pl
import numpy as np

datapath = '/Users/adam/work/gc/limabean/'
figpath = '/Users/adam/work/h2co/limabean/figures/'

pl.figure(10)
pl.clf()
F = aplpy.FITSFigure(datapath+'limabean_MIPS_24_crop.fits',
                     figure=pl.figure(10),
                     convention='calabretta')
F.set_auto_refresh(False)
vm = [(50,None,None),]
vmin,vmid,vmax = vm[0]
F.show_grayscale(stretch='log',vmin=vmin,vmid=vmid,vmax=vmax)

#datapath+'LimaBean_H2CO11_mean_density_vr_-20to50.fits'
#datapath+'LimaBean_H2CO11_mean_density_vr_50to70.fits'
#datapath+'LimaBean_H2CO11_peak_density_vr_50to70.fits'
F.show_contour(datapath+'LimaBean_H2CO11_peak_density_vr_-20to50.fits',
               levels=np.log10([5e2,1e3,5e3,1e4,5e4,1e5,5e5]),
               colors=['c','b','g','orange','red','magenta'],
               filled=True,
               alpha=0.25)
F.recenter(0.26, 0.035, height=0.2, width=0.2)
F.set_tick_labels_xformat('dd.dd')
F.set_tick_labels_yformat('dd.dd')
B = aplpy.Beam(F)
B.show(major=2.6/60/2.35,minor=2.6/60/2.35,angle=0)
B.set(facecolor='orange',edgecolor='orange',alpha=0.5)
B.set_corner('top right')

F.refresh()

F.save(figpath+'peakdensity_contours_on_24um.pdf')

