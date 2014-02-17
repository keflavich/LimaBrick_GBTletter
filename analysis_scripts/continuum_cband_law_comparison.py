from paths import figpath,datapath
from astropy.io import fits
import FITS_tools
import os
import aplpy


#lawfile = os.path.join(datapath,'GCCBand_Brick_zoom_crop.fits')
lawfile = os.path.join(datapath,'GCCBand_lb.2.fits')
lawdata = fits.getdata(lawfile).squeeze()
law_hdr = FITS_tools.strip_headers.flatten_header(fits.getheader(lawfile))

ccontfile = os.path.join(datapath,'LimaBean_H2CO11_cube_continuum.fits')
ccontdata = fits.getdata(ccontfile).squeeze()
ccont_hdr = fits.getheader(ccontfile)

okimg = np.isfinite(ccontdata).astype('float')

resampled_law = FITS_tools.hcongrid.hcongrid(lawdata, law_hdr, ccont_hdr)

xok = (np.isfinite(resampled_law) * (resampled_law != 0)).max(axis=0)
yok = (np.isfinite(resampled_law) * (resampled_law != 0)).max(axis=1)

croprange_x = np.argmax(xok),np.argmax(xok)+xok.sum()
croprange_y = np.argmax(yok),np.argmax(yok)+yok.sum()
slice_x = slice(*croprange_x)
slice_y = slice(*croprange_y)

cropped_ccont = ccontdata[slice_y,slice_x]
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
for F in (F1,F2):
    F.show_grayscale(vmin=0.0, vmax=20, stretch='log', vmid=-0.5, invert=True)
    F.set_tick_labels_xformat('dd.d')
    F.set_tick_labels_yformat('dd.d')
    F.set_tick_xspacing(0.3)
    F.show_contour(fits.PrimaryHDU(okimg,law_hdr), convention='calabretta', levels=[0.5], colors=['r'])
    F.add_colorbar()
F1.remove_colorbar() # hack to preserve figure size
F2.axis_labels.hide_y()
F2.tick_labels.hide_y()
F2.colorbar.set_ticks([0.1,0.5,1,2.5,5,10,20])


pl.subplot(2,1,2)
ok = np.isfinite(cropped_law) # ccont is all "ok"
pl.plot(np.linspace(0,20),np.linspace(0,20),'k--',linewidth=2,alpha=0.5)
pl.plot(cropped_ccont[ok],cropped_law[ok],',')
#pl.plot(cropped_ccont[np.round(cropped_rsok).astype('bool')],cropped_law[np.round(cropped_rsok).astype('bool')],'.',color='r')
#mpl_plot_templates.adaptive_param_plot(cropped_ccont[ok],cropped_law[ok],bins=30,threshold=10,fill=True)
pl.xlabel(r'$T_B(K)$ Ours')
pl.ylabel(r'$T_B(K)$ Law')
pl.axis([0,20,0,20])
pl.savefig(figpath+'comparison_Law_to_Cband.pdf')

pl.show()
