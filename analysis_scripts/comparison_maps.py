from paths import figpath,datapath
import aplpy
import pylab as pl
#import montage

pl.close(9)
pl.figure(9)
pl.clf()

#http://irsa.ipac.caltech.edu/data/SPITZER/GLIMPSE/images/II/1.2_mosaics_v3.5/GLM_00000+0000_mosaic_I4.fits
#http://irsa.ipac.caltech.edu/data/SPITZER/MIPSGAL/images/mosaics24/MG0000p005_024.fits
#http://irsa.ipac.caltech.edu/data/SPITZER/MIPSGAL/images/mosaics24/MG0000n015_024.fits
#for fullfn,cropfn in [('GCCBand_lb.2.fits','GCCBand_Brick_zoom_crop.fits'),
#                      ('GCXBand2_lb.2.fits','GCXBand_Brick_zoom_crop.fits'),
#                      TODO - THIS IS TEH FAILS
#                     ]

Figures = [aplpy.FITSFigure(datapath+'limabean_MIPS_24_crop.fits',
                            figure=pl.figure(9), convention='calabretta',
                            #subplot=(2,2,1)),
                            subplot=[0.1,0.55,0.45,0.45]),
           aplpy.FITSFigure(datapath+'GCCBand_Brick_zoom_crop.fits',
                            figure=pl.figure(9), convention='calabretta',
                            #subplot=(2,2,2)),
                            subplot=[0.55,0.55,0.45,0.45]),
           aplpy.FITSFigure(datapath+'GCXBand_Brick_zoom_crop.fits',
                            figure=pl.figure(9), convention='calabretta',
                            #subplot=(2,2,3)),
                            subplot=[0.1,0.1,0.45,0.45]),
           aplpy.FITSFigure(datapath+'LimaBean_H2CO22_cube_continuum.fits',
                            figure=pl.figure(9), convention='calabretta',
                            #subplot=(2,2,4)),
                            subplot=[0.55,0.1,0.45,0.45]),
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
    
for ii in (1,3):
    Figures[ii].axis_labels.hide_y()
    Figures[ii].tick_labels.hide_y()
for ii in (0,1):
    Figures[ii].axis_labels.hide_x()
    Figures[ii].tick_labels.hide_x()
pl.subplots_adjust(hspace=0,wspace=0)

Figures[3].show_regions('./data/brick_circles_nooffs.reg')
Figures[3].refresh()

F.save(figpath+'continuum_comparison.pdf', dpi=85)

pl.figure(10)
pl.clf()
F = aplpy.FITSFigure(datapath+'limabean_MIPS_24_crop.fits',
                     figure=pl.figure(10),
                     convention='calabretta')
vmin,vmid,vmax = vm[0]
F.show_grayscale(stretch='log',vmin=vmin,vmid=vmid,vmax=vmax)
F.show_regions('./data/brick_circles_nooffs.reg')
F.refresh()
F.save(figpath+'apertures_on_24um.pdf')
