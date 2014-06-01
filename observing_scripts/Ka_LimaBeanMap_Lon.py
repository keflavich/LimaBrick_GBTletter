"""
Create a 12x15' map of the Lima Bean at Ka band

*This observation will take 2.4 hours on-source, *PLUS* overheads, so about 2.6 hours*

Per Glen Langston's recommendations, I should:
    1. not use refs for each scan (because the pipeline doesn't use them)
    2. use "Track" to observe offs
    3. do pointing scans as normal; break up map into 2 sub-maps

Following pipeline recommendations:
https://safe.nrao.edu/wiki/bin/view/Kbandfpa/ObserverGuide?sortcol=table;up=#Reduction_Execute_Pipeline_with
"""

Break("Make sure you run ConfigureFocusKa.py before beginning this observation")

cat = Catalog("/users/aginsbur/GBT12B-221/limabean.astrid")
Configure("/users/aginsbur/GBT12B-221/H2CO_1cm_KaSetup_GC.py")


Slew("LimaBean")
Balance()

amintodeg = 1/60.
# samplerate = 30/minute = 0.5/s
# beam ~ 25"
# 4 samples/beam
# 8 seconds / arcminute
scanrate = 7.5 # arcmin/min = arcsec/sec
scanheight = 12. # arcmin
scanwidth = 15. # arcmin
# nscans = 6 * 15 = 90

# horizontal scans
RALongMap('LimaBean',     #center of map
    hLength = Offset("Galactic",scanwidth*amintodeg,0.0,cosv=True), 
    vLength = Offset("Galactic",0.0,scanheight*amintodeg,cosv=True), 
    vDelta  = Offset("Galactic",0.0,(1./6.)*amintodeg,cosv=True), 
    scanDuration = scanwidth/scanrate * 60,
    beamName="1")
