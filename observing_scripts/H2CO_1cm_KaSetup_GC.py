# Ka-band H2CO 1 cm (and friends)
# keywords defined at www.gb.nrao.edu/~fghigo/gbtdoc/config/configparams_rev1.doc
#
# CONFIGURATION:
# Bandwidth Polarization Level Number of Number of Lags - Approximate Resolution
#   (MHz) Cross-Products Sampling Spectral Beams Low Medium High
# BW Pol Lev  Windows Beams Channels / resolution
# 50 No  9    4       2     4096 - 12.2070 kHz 4096 - 12.2070 kHz 4096 - 12.2070 kHz

receiver  = 'Rcvr26_40'               # select Ka-band receiver
beam      = 'B1'                      # use two beams: NO! Memo 255 says only B1
# https://safe.nrao.edu/wiki/pub/GB/Knowledge/GBTMemos/GBT_Memo_255.pdf
obstype   = 'Spectroscopy'
backend   = 'Spectrometer'
nwin      = 4                         # eight spectral windows (b/c 1 pol)
# H2CO, H213CO, H2C18O, H62a; max separation 4000 MHz
# HCOOH (Formic acid), SO, ccs, HC3N
# CCS 2-1: 29.4777
# HCCN: 28.65117 through 28.75931	
# CCCS 5-4: 28.90369	
restfreq  = 28974.8,27555.67,26330.14,26939.16#,28.08636,30.00154,29.4777,27.29429
deltafreq = 0,0,0,0#,0,0,0,0           # DO NOT MENTION IF3FREQ IN SETUP!
bandwidth = 50.0                      # MHz Moderate-resolution mode (0.25 km/s)
swmode    = "tp"                      # set switching scheme (tp(total power with cal), tp_nocal, sp(switched power with cal), sp_nocal )
swtype    = "none"                    # for frequency switching; not used for tp mode
swper     = 1.0                       # one second cycle for switching
swfreq    = 0.0, 0.0                  # for freq switching
tint      = 1.2                       # integration time (for 4 quadrants, 1.2-40 sec.  Important to avoid smearing)
vlow      = 0
vhigh     = 0
vframe    = "lsrk"                    # LSR - kinematic is the "normal" definition (don't use dynamic)
vdef      = "Radio"                   # radio (optical is also acceptable, but not the norm for Galactic observations)
noisecal  = "lo"
pol       = "Circular"
nchan     = "high"                     # 4096 channels over 2 beams and 4 windows
        # spectrometer guide says 12.5 MHz = 237 km/s bandwidth = 8192 channels,
        # 50 MHz = 535 km/s BW
        # 4096 channels, 12.2 KHz = 0.12 km/s
        # numsamplers=8, nwin=4, chanwidth = 1.526 KHz -> resolution = .06 km/s (2 channels)... nchan = high
spect.levels = 9                      # nine level sampling
#iftarget = 0.25                       # IF target is now always 1
