import pyspeckit
import numpy as np
import sys
import h2co_modeling
import pylab as pl
import itertools

# Multiplication factor for the optical depth to account for the larger
# velocity gradient (8 km/s/pc) in The Brick
# See pyradex/notes.rst to confirm that tau propto 1/deltav
# the OBSERVED optical depth should also be scaled up by 1/fillingfactor
taumultfact = 1/8.
fillingfactor = 0.58
vrange = [-20,50]

pl.rcParams['font.size'] = 20

master_linestyles = itertools.cycle(['-','--','-.',':'])
master_linestyles = ['-']*10
master_colors = (['#880000', '#008888', '#9933CC']+
                 ["#"+x for x in '348ABD, CCCCCC, A60628, 467821, CF4457, 188487, E24A33'.split(', ')])

radexdatapath = '/Users/adam/work/h2co/modeling_paper/radex_data/' 

faur = h2co_modeling.SmoothtauModels(datafile=radexdatapath+'faure/1-1_2-2_XH2CO_fixed_faure.dat')

datapath = '/Users/adam/work/gc/limabean/'

rbg = pyspeckit.Spectrum(datapath+'Brick_RatioSpectrum_BGsub.fits')
rraw = pyspeckit.Spectrum(datapath+'Brick_RatioSpectrum.fits')

h2co11 = pyspeckit.Spectrum(datapath+'Brick_H2CO11_smooth_BGsub.fits')
h2co22 = pyspeckit.Spectrum(datapath+'Brick_H2CO22_smooth_BGsub.fits')


abund = -9
for abund in (-8,-8.5,-9,-10):

    pl.figure(1)
    pl.clf()
    pl.title("$X=10^{{{}}}$".format(abund))
    pl.figure(2)
    pl.clf()
    pl.title("$X=10^{{{}}}$".format(abund))

    colors = iter(master_colors)
    linestyles = iter(master_linestyles)
    # Need T=50 models?  Have now!
    for sigma in (1.0, 2.0):
        # OPR = 0.1 corresponds to 38K:
        # OPR ~ 9 exp (-170.5/T) [Faure priv. comm]
        tau1,tau2,dens,col = faur.select_data(abund, temperature=50, opr=0.1)
        tau,vtau,vtau_ratio = faur.generate_tau_functions(abundance=abund, temperature=50, opr=0.1)

        tauratio = vtau_ratio(dens, line1=tau1, line2=tau2, sigma=sigma)
        tauA = vtau(dens, line=tau1, sigma=sigma)
        tauB = vtau(dens, line=tau2, sigma=sigma)

        C = color = colors.next()

        #ok = np.arange(tauratio.size) > np.argmax(tauratio)

        #def ratio_to_dens(ratio):
        #    inds = np.argsort(tauratio[ok])
        #    return np.interp(ratio, tauratio[ok][inds], dens[ok][inds], np.nan, np.nan)

        pl.figure(1)
        pl.plot(dens,tauratio,label='$\sigma=%0.1f$' % (sigma), linewidth=3, alpha=0.7,
                #linestyle=linestyles.next(),
                color=color)

        pl.figure(2)

        pl.plot(dens,tauA*taumultfact,label='$\sigma=%0.1f$' % (sigma), linewidth=3, alpha=0.7,
                #linestyle=linestyles.next(),
                color=C)
        pl.plot(dens,tauB*taumultfact, linewidth=3, alpha=0.7,
                linestyle='--',
                color=C)

    C = colors.next()

    pl.figure(1)
    pl.plot(dens,tau1/tau2,label='$\delta$', linewidth=3, alpha=0.7,
                linestyle=linestyles.next(),
                color=C)

    pl.figure(2)
    pl.plot(dens,tau1*taumultfact,label='$\delta$', linewidth=3, alpha=0.7,
                linestyle=linestyles.next(),
                color=C)
    pl.plot(dens,tau2*taumultfact, linewidth=3, alpha=0.7,
                linestyle='--',
                color=C)


    pl.figure(1)
    pl.grid()

    pl.axis([0,7,0,15])
    pl.legend(loc='best')
    pl.xlabel(r'Volume-averaged density $\log(n(H_2))$')
    pl.ylabel(r'Ratio $\tau_{1-1}/\tau_{2-2}$')

    pl.savefig('/Users/adam/work/h2co/limabean/figures/tau_ratio_vs_density_thinlimit_sigmavary_Xm{}.pdf'.format(abs(abund)),bbox_inches='tight')

    rvals_bg = rbg.slice(vrange[0],vrange[1],units='km/s').data
    lc = pl.hlines(rvals_bg,dens.min(),dens.max(),color='b',alpha=0.1,linewidth=10)
    pl.savefig('/Users/adam/work/h2co/limabean/figures/tau_ratio_vs_density_thinlimit_sigmavary_Xm{}_withBrickBGsubdata.pdf'.format(abs(abund)),bbox_inches='tight')

    lc.set_visible(False)
    rvals_raw = rraw.slice(vrange[0],vrange[1],units='km/s').data
    pl.hlines(rvals_raw,dens.min(),dens.max(),color='g',alpha=0.1,linewidth=10)
    pl.savefig('/Users/adam/work/h2co/limabean/figures/tau_ratio_vs_density_thinlimit_sigmavary_Xm{}_withBrickRawdata.pdf'.format(abs(abund)),bbox_inches='tight')


    pl.figure(2)
    pl.grid()

    pl.axis([0,7,0.001,10])
    pl.gca().set_yscale('log')
    pl.legend(loc='best')
    pl.xlabel(r'Volume-averaged density $\log(n(H_2))$')
    pl.ylabel(r'$\tau$')

    pl.savefig('/Users/adam/work/h2co/limabean/figures/tau_vs_density_thinlimit_sigmavary_Xm{}.pdf'.format(abs(abund)),bbox_inches='tight')

    h11vals = h2co11.slice(vrange[0],vrange[1],units='km/s').data / fillingfactor
    h22vals = h2co22.slice(vrange[0],vrange[1],units='km/s').data / fillingfactor
    lc11 = pl.hlines(h11vals,dens.min(),dens.max(),color='k',alpha=0.05,linewidth=5,linestyle='-')
    lc22 = pl.hlines(h22vals,dens.min(),dens.max(),color='r',alpha=0.1,linewidth=5,linestyle='--')
    pl.savefig('/Users/adam/work/h2co/limabean/figures/tau_vs_density_thinlimit_sigmavary_Xm{}_withBrickBGsubdata.pdf'.format(abs(abund)),bbox_inches='tight')

pl.show()

