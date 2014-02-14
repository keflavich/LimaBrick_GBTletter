import matplotlib
import pyspeckit
from astropy.io import fits
import pylab as pl

datadir_C = '/Users/adam/observations/gbt/AGBT12B_221_01/'
datadir_U = '/Users/adam/observations/gbt/AGBT14A_110_01/'
figpath = '/Users/adam/work/h2co/limabean/figures/'

matplotlib.rc_file('/Users/adam/.matplotlib/ggplotrc')

colors = ['k','r','g','b', 'm', "#348ABD", "#7A68A6"]

def load_off_spectra(prefix='AGBT12B_221_01_A9',datadir=datadir_C):
    sp0 = pyspeckit.Spectrum(datadir+prefix+"_offspectra.fits",specnum=0)
    sp0_nom = pyspeckit.Spectrum(datadir+prefix+"_offspectra.fits",specnum=1)

    spectra = [sp0,sp0_nom]

    return spectra

def make_off_plots(spectra, dobaseline=True,
                   vcut1=[-70,-10],
                   vcut2=[-30,110]):

    ax1 = pl.subplot(2,1,1)
    ax2 = pl.subplot(2,2,3)
    ax3 = pl.subplot(2,2,4)

    xmin,xmax = -150,200

    if dobaseline:
        blspectra = []
        for s,c in zip(spectra,colors):
            #s.baseline.selectregion(xmin=-100,xmax=200, reset=True)
            s.plotter.autorefresh=False
            s.plotter.axis = None # just to make damned sure
            # then fit the baseline so we can see a flattened version
            s.baseline(xmin=xmin,xmax=xmax, reset=True, exclude=[-50,100],
                       annotate=False, order=5, subtract=False)
            s2 = pyspeckit.Spectrum(xarr=s.xarr, data=s.data/s.baseline.basespec-1)
            # Set the background level to zero...
            s.baseline(xmin=xmin,xmax=xmax, reset=True, exclude=[-50,100],
                       annotate=False, order=0, subtract=True)
            s.unit = ""
            s2.unit = ""
            blspectra.append(s2)
            #s.baseline.highlight_fitregion()

    if dobaseline:
        for s2,c2 in zip(blspectra, colors[::-1]):
            s2.plotter(axis=ax1, clear=False, color=c2, xmin=xmin,
                       xmax=xmax, alpha=0.5, linewidth=2)
            s2.plotter(axis=ax2, clear=False, color=c2,
                       xmin=vcut1[0], xmax=vcut1[1],
                       alpha=0.5, linewidth=2)
            s2.plotter(axis=ax3, clear=False, color=c2,
                       xmin=vcut2[0], xmax=vcut2[1],
                       alpha=0.5, linewidth=2)

    for s,c in zip(spectra, colors):
        s.plotter(axis=ax1, clear=False, color=c, xmin=xmin,
                  xmax=xmax, alpha=0.5, linewidth=2)
        s.plotter(axis=ax2, clear=False, color=c,
                  xmin=vcut1[0], xmax=vcut1[1],
                  alpha=0.5, linewidth=2)
        s.plotter(axis=ax3, clear=False, color=c,
                  xmin=vcut2[0], xmax=vcut2[1],
                  alpha=0.5, linewidth=2)

    for ax in (ax1,ax2,ax3):
        ax.set_ylabel("")

def do_C():
    all_spectra = {}
    for ref1,ref2 in ([6,21],[22,32]):
        for ii,sampler in enumerate(('A9','A13')):
            pl.figure(ii)
            pl.clf()
            all_spectra['C_'+sampler] = load_off_spectra('AGBT12B_221_01_{0}_sr{1}-{2}'.format(sampler,ref1,ref2),
                                                         datadir=datadir_C)
            make_off_plots(all_spectra['C_'+sampler], dobaseline=False)
            all_spectra['C_'+sampler][0].plotter.savefig(figpath+"offspectra_interpolation_C_{0}_sr{1}-{2}.pdf".format(sampler,ref1,ref2))

def do_U():
    all_spectra = {}
    for ref1,ref2 in ([9,54],[62,98],[108,140]):
        for ii,sampler in enumerate(('A9_fd1_if0','A13_fd1_if0','C25_fd2_if0','C29_fd2_if0')):
            pl.figure(ii+2)
            pl.clf()
            all_spectra['U_'+sampler] = load_off_spectra('AGBT14A_110_01_{0}_sr{1}-{2}'.format(sampler,ref1,ref2)
                                                         ,datadir=datadir_U)
            make_off_plots(all_spectra['U_'+sampler], dobaseline=True)
            all_spectra['U_'+sampler][0].plotter.savefig(figpath+"offspectra_interpolation_U_{0}_sr{1}-{2}.pdf".format(sampler,ref1,ref2))

    return all_spectra

if __name__ == "__main__":
    all_spectra_C = do_C()
    #all_spectra_U = do_U()
