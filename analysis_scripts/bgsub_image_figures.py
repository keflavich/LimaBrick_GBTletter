import aplpy
import pylab as pl

datapath = '/Users/adam/work/gc/limabean/'

f = pl.figure(9,figsize=[12,4])
if f.get_figheight() != 4 or f.get_figwidth != 12:
    pl.close(9)
    pl.figure(9,figsize=[12,4])
else:
    pl.clf()

Figures = [
    aplpy.FITSFigure('/Users/adam/work/gc/limabean/limabean_MIPS_24_crop.fits',
                     figure=pl.figure(9),convention='calabretta', subplot=[0.1,0.1,0.3,0.8]),
    aplpy.FITSFigure('/Users/adam/work/gc/limabean/LimaBean_H2CO11_taucube_integrated.fits',
                     figure=pl.figure(9),convention='calabretta', subplot=[0.4,0.1,0.3,0.8]),
    aplpy.FITSFigure('/Users/adam/work/gc/limabean/LimaBean_H2CO22_taucube_integrated.fits',
                     figure=pl.figure(9),convention='calabretta', subplot=[0.7,0.1,0.3,0.8]),
    ]

vm = [('log',50,None,None),
      ('linear',None,None,None),#(0,-0.1,7.5),
      ('linear',None,None,None),#(0,-0.1,7.5),
      ]

for F,(stretch,vmin,vmid,vmax) in zip(Figures,vm)[::-1]:
    F.set_auto_refresh(False)
    F.set_tick_labels_xformat('dd.d')
    F.set_tick_labels_yformat('dd.d')
    F.set_tick_labels_size(20)
    F.axis_labels.set_font(size=20)
    F.show_grayscale(stretch=stretch,vmin=vmin,vmid=vmid,vmax=vmax)
    F.recenter(0.253, 0.016, height=0.2, width=0.2)
    F.show_contour(datapath+'brick_aperture11.fits', alpha=0.2, levels=[0.5,1.5], colors=['r']*3, filled=True)

for F in Figures[1:]:
    F.ticks.hide_y()
    F.tick_labels.hide_y()
    F.axis_labels.hide_y()

Figures[0]._ax1.set_title('24 $\mu$m')
Figures[1]._ax1.set_title('H$_2$CO 1-1')
Figures[2]._ax1.set_title('H$_2$CO 2-2')

for F in Figures:
    F.refresh()


Figures[0].save('/Users/adam/work/h2co/limabean/figures/bgsub_aperture_mask.pdf')

pl.show()
