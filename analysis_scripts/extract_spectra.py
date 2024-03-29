from paths import datapath,outpath
from paths import figpath as figurepath
import os
import pyspeckit
import pyregion
import pylab as pl

scube11 = pyspeckit.Cube(datapath+'LimaBean_H2CO11_cube_sub.fits')
taucube11 = pyspeckit.Cube(datapath+'LimaBean_H2CO11_taucube_claw.fits')
scube11.xarr.convert_to_unit('km/s')
taucube11_13 = pyspeckit.Cube(datapath+'LimaBean_H213CO_taucube_claw.fits')
scube11_13 = pyspeckit.Cube(datapath+'LimaBean_H213CO_cube_sub.fits')
scube11_13.xarr.convert_to_unit('km/s')
scube22 = pyspeckit.Cube(datapath+'LimaBean_H2CO22_cube_sub_smoothtoCband.fits')
taucube22 = pyspeckit.Cube(datapath+'LimaBean_H2CO22_taucube_smoothtoCband.fits')
scube22.xarr.convert_to_unit('km/s')

scube33 = pyspeckit.Cube(datapath+'LimaBean_H2CO33_cube_sub_smoothtoCband.fits')
taucube33 = pyspeckit.Cube(datapath+'LimaBean_H2CO33_taucube_smoothtoCband.fits')
scube33.xarr.convert_to_unit('km/s')

if os.path.exists(datapath+'Brick_HOPS_MOPRA_NH3_22_cube.fits'):
    scubenh3 = pyspeckit.Cube(datapath+'Brick_HOPS_MOPRA_NH3_22_cube.fits')
else:
    scubenh3 = pyspeckit.Cube(datapath+'NH3_rotcrop.fits')
scubenh3.xarr.convert_to_unit('km/s')

regfn = './data/brick_circles.reg'
regions = pyregion.open(regfn)

spectra11 = {}
spectra11_13 = {}
spectra22 = {}
spectra33 = {}
spectranh3 = {}

# Loading
for ii, reg in enumerate(regions):
    sp11 = scube11.get_apspec(reg.coord_list, coordsys='galactic', wunit='degree')
    sp22 = scube22.get_apspec(reg.coord_list, coordsys='galactic', wunit='degree')
    sp33 = scube33.get_apspec(reg.coord_list, coordsys='galactic', wunit='degree')
    sp11_13 = scube11_13.get_apspec(reg.coord_list, coordsys='galactic', wunit='degree')
    taup11 = taucube11.get_apspec(reg.coord_list, coordsys='galactic', wunit='degree')
    taup22 = taucube22.get_apspec(reg.coord_list, coordsys='galactic', wunit='degree')
    taup33 = taucube33.get_apspec(reg.coord_list, coordsys='galactic', wunit='degree')
    taup11_13 = taucube11_13.get_apspec(reg.coord_list, coordsys='galactic', wunit='degree')

    spnh3 = scubenh3.get_apspec(reg.coord_list, coordsys='galactic', wunit='degree')
    spnh3.smooth(3,downsample=False)
    spnh3.data *= sp11.data.max()/spnh3.data.max()
    spnh3.units = 'Normalized'
    sp11.specname = reg.attr[1]['text']
    sp22.specname = reg.attr[1]['text']
    sp33.specname = reg.attr[1]['text']
    spnh3.specname = reg.attr[1]['text']
    taup11.specname = reg.attr[1]['text']
    taup22.specname = reg.attr[1]['text']
    taup33.specname = reg.attr[1]['text']
    spectra11[ii] = taup11
    spectra11_13[ii] = taup11_13
    spectra22[ii] = taup22
    spectra33[ii] = taup33
    spectranh3[ii] = spnh3
    sp11.write(outpath+'H2CO_11_%s.fits' % sp11.specname.replace(" ","_"))
    sp22.write(outpath+'H2CO_22_%s.fits' % sp22.specname.replace(" ","_"))
    sp33.write(outpath+'H2CO_33_%s.fits' % sp33.specname.replace(" ","_"))
    sp11_13.write(outpath+'H213CO_11_%s.fits' % sp11.specname.replace(" ","_"))

nh3offset = 0.1
twotwooffset = -0.025
thirtoffset = -0.05

# Plotting
for ii in spectra11:
    sp11,sp11_13,sp22,sp33,spnh3 = [s[ii] for s in (spectra11,spectra11_13,spectra22,spectra33,spectranh3)]

    fig = pl.figure(ii)
    pl.clf()
    sp11.plotter(figure=fig, clear=True,zorder=50,xmin=-75,xmax=150, linewidth=2,
                 alpha=0.6)
    sp11_13.plotter(axis=sp11.plotter.axis, clear=False, color='b',
                    alpha=0.5, offset=thirtoffset, use_window_limits=True,zorder=20,
                    linewidth=2)
    sp22.plotter(axis=sp11.plotter.axis, clear=False, color='r',
                 offset=twotwooffset, use_window_limits=True,zorder=25,
                 linewidth=2)
    sp33.plotter(axis=sp11.plotter.axis, clear=False, color='orange',
                 offset=twotwooffset, use_window_limits=True,zorder=25,
                 linewidth=2)
    spnh3.plotter(axis=sp11.plotter.axis, clear=False, color='g',
                  offset=nh3offset, use_window_limits=True,zorder=10,
                  linewidth=2)
    sp11.plotter.axis.hlines([twotwooffset,0.0,nh3offset],
                             -100,100,color='purple',linestyle='--',alpha=0.5,zorder=5,
                             linewidth=2)
    sp11.plotter.axis.set_ylim(min([sp22.data.min()+twotwooffset,sp11_13.data.min()+thirtoffset]),
                               max([spnh3.data.max()+nh3offset,sp11.data.max()]))

    sp11.plotter.refresh()

    ratio = sp11.copy()
    ratio.data = sp11.data/sp22.data
    ratio.units=r'$\tau$ ratio'
    ratio2 = sp11.copy()
    ratio2.data = sp11.data/sp11_13.data
    ratio2.units=r'$\tau$ ratio'
    ratio3 = sp11.copy()
    ratio3.data = sp11_13.data*25 / sp22.data
    ratio3.units=r'$\tau$ ratio'
    ratio4 = sp11.copy()
    sp33i = pyspeckit.interpolation.interp_on_axes(sp33, sp22.xarr)
    ratio4.data = sp22.data / sp33i.data
    ratio4.units=r'$\tau$ ratio'

    inset = pl.axes([0.65,0.65,0.25,0.25])
    xmin,xmax = -20,60
    ratio.plotter(axis=inset, xmin=xmin, xmax=xmax, ymin=0,ymax=20, 
                  linewidth=2, alpha=0.6)
    ratio2.plotter(axis=inset, xmin=xmin, xmax=xmax, ymin=0,ymax=20,
                   clear=False, color='b', alpha=0.5, linewidth=2)
    ratio3.plotter(axis=inset, xmin=xmin, xmax=xmax, ymin=0,ymax=20,
                   clear=False, color='r', alpha=0.5, linewidth=2)
    ratio4.plotter(axis=inset, xmin=xmin, xmax=xmax, ymin=0,ymax=20,
                   clear=False, color='orange', alpha=0.5, linewidth=2)

    if 'off' not in sp11.specname.lower():
        sp11.plotter.savefig(figurepath+sp11.specname.replace(" ","_")+"_11_22_nh3_spectra.pdf")

