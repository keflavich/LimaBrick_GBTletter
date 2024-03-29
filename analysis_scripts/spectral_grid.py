from paths import datapath,figpath
import pyspeckit
import itertools
import astropy.io.fits as pyfits
from astropy import wcs
import sys
import pylab as pl
import numpy as np
from FITS_tools import strip_headers
from sdpy import makecube
do_some_plots=True
np.seterr(all='ignore')

def spectral_grid(cube11=pyspeckit.Cube(datapath+'LimaBean_H2CO11_cube_sub.fits'),
                  h213co11cube=pyspeckit.Cube(datapath+'LimaBean_H213CO_cube_sub.fits'),
                  cube22=pyspeckit.Cube(datapath+'LimaBean_H2CO22_cube_sub_smoothtoCband.fits'),
                  cube33=pyspeckit.Cube(datapath+'LimaBean_H2CO33_cube_sub_smoothtoCband.fits'),
                  figure=pl.figure(1,figsize=(10,10)),
                  yrange=(-2.1,0.5),
                  nx=6,
                  ny=6,
                  dx=2,
                  dy=2,
                  xc=48,
                  yc=45,
                  xlabel="Velocity",
                  ylabel="$T_A^*$ [K]$",
                  ratio=False):

    for c in [cube11,h213co11cube,cube22,cube33]:
        c.xarr.convert_to_unit('km/s')

    def some_plots(xx,yy,dolegend=False):
        c13 = h213co11cube.get_spectrum(xx,yy)
        c12 = cube11.get_spectrum(xx,yy)
        c22 = cube22.get_spectrum(xx,yy)
        c33 = cube33.get_spectrum(xx,yy)

        for c in (c12,c13,c22,c33):
            c.plotter.autorefresh=False

        c12.baseline(exclude=[-225,200],order=5)
        c13.baseline(exclude=[-225,200],order=5)
        c13.plotter(label='H$_{2}$$^{13}$CO 1-1', axis=pl.gca(), color='b',
                    alpha=0.5, refresh=False)
        c12.plotter(axis=c13.plotter.axis, clear=False, color='k',
                    label="H$_{2}$CO 1-1", alpha=0.8, linewidth=1.5,
                    refresh=False)
        #(c13*6).plotter(label='6$\\times$H$_{2}$$^{13}$CO',axis=pl.gca(),color='r',clear=False)
        c22.plotter(axis=c13.plotter.axis,clear=False,color='r',linewidth=2,alpha=0.5,
                    label='H$_2$CO 2-2', refresh=False)
        c33.plotter(axis=c13.plotter.axis,clear=False,color='orange',linewidth=2,alpha=0.5,
                    label='H$_2$CO 3-3', refresh=False)

        if ratio:
            r = c12.copy()
            r.data = c22.data/c12.data
            r.data[(r.data>1)] = np.nan
            r.data[(r.data<1/13.)] = np.nan
            r.plotter(axis=c13.plotter.axis,clear=False,color='r', refresh=False)
        if dolegend:
            pl.legend(loc='best')
        pl.gca().set_xlim(-100,100)
      

    pl.clf()

    pl.subplots_adjust(wspace=0,hspace=0,left=0.05,right=0.95,bottom=0.05,top=0.95)

    #nx,ny = 6,6
    #xc,yc = 52,48 # center
    #xc,yc = 45,42 # bottom left
    xl = xc - (nx/2.)*dx
    yl = yc - (ny/2.)*dy
    
    w = wcs.WCS(strip_headers.flatten_header(cube11.header))
    for ii,(spy,spx) in enumerate(itertools.product(range(nx),range(ny))):
        sp = pl.subplot(nx,ny,ii+1)
        sp.zorder = nx-spx+ny*spy
        some_plots(spx*dx+xl,(ny-spy-1)*dy+yl,dolegend=False) # ,dolegend=(ii==24))
        sp.set_ylim(*yrange)
        if spx > 0:
            sp.set_yticks([])
        else:
            sp.set_yticks(sp.get_yticks()[1:])
        if spy < ny-1:
            sp.set_xticks([])
        else:
            sp.set_xticks(sp.get_xticks()[1:])
        sp.set_xlabel("")
        sp.set_ylabel("")

        (l,b), = w.wcs_pix2world([[spx*dx+xl,(ny-spy-1)*dy+yl]],0)
        sp.annotate('%0.3f %+0.3f' % (l,b),(0.5,0.9),xycoords='axes fraction',fontsize=12,
                    horizontalalignment='center')
        sp.xaxis.set_tick_params(labelsize=14)
        sp.yaxis.set_tick_params(labelsize=14)

    pl.figtext(0.5,0.01,xlabel)
    pl.figtext(-0.01,0.5,ylabel,rotation='vertical')

    pl.draw()

spectral_grid()
pl.savefig(figpath+'spectralgrid_absorption.pdf')

spectral_grid(cube11=pyspeckit.Cube(datapath+'LimaBean_H2CO11_taucube_claw.fits'),
              cube22=pyspeckit.Cube(datapath+'LimaBean_H2CO22_taucube_smoothtoCband.fits'),
              cube33=pyspeckit.Cube(datapath+'LimaBean_H2CO33_taucube_smoothtoCband.fits'),
              h213co11cube=pyspeckit.Cube(datapath+'LimaBean_H213CO_taucube_claw.fits'),
              figure=pl.figure(2,figsize=(10,10)),
              yrange=(-0.05,0.30),
              ylabel=r'$\tau_{observed}$',
              xlabel=r'$V_{LSR}$ km s$^{-1}$',
              ratio=False)
pl.savefig(figpath+'spectralgrid_optdepth.pdf')

spectral_grid(cube11=pyspeckit.Cube(datapath+'LimaBean_H2CO11_taucube_claw.fits'),
              cube22=pyspeckit.Cube(datapath+'LimaBean_H2CO22_taucube_smoothtoCband.fits'),
              cube33=pyspeckit.Cube(datapath+'LimaBean_H2CO33_taucube_smoothtoCband.fits'),
              h213co11cube=pyspeckit.Cube(datapath+'LimaBean_H213CO_taucube_claw.fits'),
              figure=pl.figure(3,figsize=(10,10)),
              yrange=(-0.05,0.30),
              dx=4,dy=4,
              xc=51,
              yc=49,
              ylabel=r'$\tau_{observed}$',
              xlabel=r'$V_{LSR}$ km s$^{-1}$',
              ratio=False)
pl.savefig(figpath+'spectralgrid_optdepth_wide.pdf')

spectral_grid(cube11=pyspeckit.Cube(datapath+'LimaBean_H2CO11_taucube_claw.fits'),
              cube22=pyspeckit.Cube(datapath+'LimaBean_H2CO22_taucube_smoothtoCband.fits'),
              cube33=pyspeckit.Cube(datapath+'LimaBean_H2CO33_taucube_smoothtoCband.fits'),
              h213co11cube=pyspeckit.Cube(datapath+'LimaBean_H213CO_taucube_claw.fits'),
              figure=pl.figure(4,figsize=(10,10)),
              yrange=(-0.05,0.30),
              dx=6,dy=6,
              xc=51,
              yc=49,
              ylabel=r'$\tau_{observed}$',
              xlabel=r'$V_{LSR}$ km s$^{-1}$',
              ratio=False)
pl.savefig(figpath+'spectralgrid_optdepth_vwide.pdf')
