from paths import figpath,datapath
from astropy.io import fits
import FITS_tools
import os
import aplpy
import numpy as np
from astropy import units as u
import pylab as pl


#lawfile = os.path.join(datapath,'GCCBand_Brick_zoom_crop.fits')
lawfile = os.path.join(datapath,'GCCBand_lb.2.fits')
lawdata = fits.getdata(lawfile).squeeze()
law_hdr = FITS_tools.strip_headers.flatten_header(fits.getheader(lawfile))

fwhm = np.sqrt(8*np.log(2))
beam = (2*np.pi*u.deg**2*(4.241274E-02)**2/fwhm**2)
kperjy = (1*u.Jy).to(u.K, equivalencies=u.brightness_temperature(beam,4.829*u.GHz))
lawdata = lawdata * kperjy.value

ccontfile = os.path.join(datapath,'LimaBean_H2CO11_cube_continuum.fits')
print "TODO: REMOVE THE +1 HERE IT'S JUST A TEST"
ccontdata = fits.getdata(ccontfile).squeeze()
#ccontdata -= np.nanmin(ccontdata) # HACK
ccont_hdr = fits.getheader(ccontfile)

brick_mask = fits.getdata(datapath+'brick_mask11.fits').astype('bool')

okimg = np.isfinite(ccontdata).astype('float')

resampled_law = FITS_tools.hcongrid.hcongrid(lawdata, law_hdr, ccont_hdr)

xok = (np.isfinite(ccontdata) * (ccontdata != 0)).max(axis=0)
yok = (np.isfinite(ccontdata) * (ccontdata != 0)).max(axis=1)

croprange_x = np.argmax(xok),np.argmax(xok)+xok.sum()
croprange_y = np.argmax(yok),np.argmax(yok)+yok.sum()
slice_x = slice(*croprange_x)
slice_y = slice(*croprange_y)

cropped_ccont = ccontdata[slice_y,slice_x]
cropped_mask = brick_mask[slice_y,slice_x]
cropped_law = resampled_law[slice_y,slice_x]

crop_hdr = ccont_hdr.copy()
crop_hdr['CRPIX1'] -= croprange_x[0]
crop_hdr['CRPIX2'] -= croprange_y[0]

law_hdu = fits.PrimaryHDU(cropped_law, crop_hdr)
ccont_hdu = fits.PrimaryHDU(cropped_ccont, crop_hdr)

fig = pl.figure(1)
sp1 = pl.subplot(2,2,1)
sp1b = sp1.bbox._bbox._points[0].tolist() + (sp1.bbox._bbox._points[1]-sp1.bbox._bbox._points[0]).tolist()
sp2 = pl.subplot(2,2,2)
sp2b = sp2.bbox._bbox._points[0].tolist() + (sp2.bbox._bbox._points[1]-sp2.bbox._bbox._points[0]).tolist()
pl.clf()
F1 = aplpy.FITSFigure(ccont_hdu, convention='calabretta', subplot=sp1b, figure=fig)
F1._ax1.set_title("Our Data")
F2 = aplpy.FITSFigure(law_hdu, convention='calabretta', subplot=sp2b, figure=fig)
F2._ax1.set_title("Law 2008")
vmin = -1
vmax = 37
for F in (F1,F2):
    F.show_grayscale(vmin=vmin, vmax=vmax, stretch='log', vmid=vmin-0.5, invert=True)
    F.set_tick_labels_xformat('dd.d')
    F.set_tick_labels_yformat('dd.d')
    F.set_tick_xspacing(0.1)
    F.recenter(0.253,0.010,0.125)
    #F.show_contour(fits.PrimaryHDU(okimg,ccont_hdr), convention='calabretta', levels=[0.5], colors=['r'])
    F.show_contour(fits.PrimaryHDU(brick_mask.astype('int'), ccont_hdr), convention='calabretta', levels=[0.5], colors=['b'])
    F.add_colorbar()
F1.remove_colorbar() # hack to preserve figure size
F2.axis_labels.hide_y()
F2.tick_labels.hide_y()
F2.colorbar.set_ticks([0.1,0.5,1,2.5,5,10,20])


pl.subplot(2,2,3)
ok = np.isfinite(cropped_law) # ccont is all "ok"
pl.plot(np.linspace(vmin,vmax),np.linspace(vmin,vmax),'k--',linewidth=2,alpha=0.5)
for low,high in ([0,10],[10,20],[20,40]):
    cutok = ok*(cropped_law>low)*(cropped_law<high)
    pl.plot(cropped_law[cutok],cropped_ccont[cutok],',',color='r')
pl.plot(cropped_law[ok*cropped_mask],cropped_ccont[ok*cropped_mask],',',color='b')
#pl.plot(cropped_ccont[np.round(cropped_rsok).astype('bool')],cropped_law[np.round(cropped_rsok).astype('bool')],'.',color='r')
#mpl_plot_templates.adaptive_param_plot(cropped_ccont[ok],cropped_law[ok],bins=30,threshold=10,fill=True)
pl.ylabel(r'$T_B(K)$ Ours')
pl.xlabel(r'$T_B(K)$ Law')
pl.axis([vmin,vmax,vmin,vmax])

pl.subplot(2,2,4)
pl.hist(cropped_ccont[ok]/cropped_law[ok], histtype='step',
        bins=np.linspace(0.0,1.2), alpha=0.5, edgecolor='k', linewidth=4,
        zorder=10)
pl.hist(cropped_ccont[cropped_mask*ok]/cropped_law[cropped_mask*ok], 
        histtype='stepfilled', bins=np.linspace(0.0,1.2), alpha=0.5,
        edgecolor='none', facecolor='b')
#for low,high in ([0,10],[10,20],[20,40]):
#    cutok = ok*(cropped_law>low)*(cropped_law<high)
#    pl.hist(cropped_ccont[cutok]/cropped_law[cutok], histtype='stepfilled',
#            bins=np.linspace(0.0,1.2), alpha=0.5, edgecolor='none')
pl.xlabel(r'$T_B(K)$ Ours/$T_B(K)$ Law')

pl.savefig(figpath+'comparison_Law_to_Cband.pdf')

pl.show()
