from paths import figpath,datapath
import aplpy
import pylab as pl

pl.figure(9)
pl.clf()

Figures = [aplpy.FITSFigure(datapath+'limabean_MIPS_24_crop.fits',
                            figure=pl.figure(9), convention='calabretta',
                            subplot=[0.1,0.55,0.35,0.35]),
           aplpy.FITSFigure(datapath+'GCCBand_Brick_zoom_crop.fits',
                            figure=pl.figure(9), convention='calabretta',
                            subplot=[0.55,0.55,0.35,0.35]),
           aplpy.FITSFigure(datapath+'GCXBand_Brick_zoom_crop.fits',
                            figure=pl.figure(9), convention='calabretta',
                            subplot=[0.1,0.1,0.35,0.35]),
           aplpy.FITSFigure(datapath+'LimaBean_H2CO22_cube_continuum.fits',
                            figure=pl.figure(9), convention='calabretta',
                            subplot=[0.55,0.1,0.35,0.35]),
           ]

vm = [(50,None,None),
      (0,-0.5,15),
      (0,-0.1,10),
      (-0.7,-0.8,10),
      ]

for F,(vmin,vmid,vmax) in zip(Figures,vm):
    F.set_auto_refresh(False)
    F.set_tick_labels_xformat('dd.d')
    F.set_tick_labels_yformat('dd.d')
    F.show_grayscale(stretch='log',vmin=vmin,vmid=vmid,vmax=vmax)
    F.recenter(0.26, 0.035, height=0.2, width=0.2)
    F.refresh()
    

Figures[3].show_regions('./brick_circles_nooffs.reg')
Figures[3].refresh()

F.save(figpath+'continuum_comparison.pdf', dpi=85)

pl.figure(10)
pl.clf()
F = aplpy.FITSFigure(datapath+'limabean_MIPS_24_crop.fits',
                     figure=pl.figure(10),
                     convention='calabretta')
vmin,vmid,vmax = vm[0]
F.show_grayscale(stretch='log',vmin=vmin,vmid=vmid,vmax=vmax)
F.show_regions('./brick_circles_nooffs.reg')
F.refresh()
F.save(figpath+'apertures_on_24um.pdf')
