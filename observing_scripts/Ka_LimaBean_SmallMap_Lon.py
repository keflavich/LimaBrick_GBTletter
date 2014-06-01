"""
Create a 8x8' map of the Lima Bean at Ka band

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
Slew("LimaBeanOff")
Track("LimaBeanOff",None,60)

amintodeg = 1/60.
arcsectodeg = 1/3600.
# samplerate = 30/minute = 0.5/s
# beam ~ 25"
# 4 samples/beam
# 8 seconds / arcminute
scanrate = 9. # arcmin/min
scanheight = 8. # arcmin
scanwidth = 8. # arcmin
vdelta = 10.0 * arcsectodeg
# nscans = 6 * 15 = 90

# horizontal scans
RALongMap('LimaBean',     #center of map
    hLength = Offset("Galactic",scanwidth*amintodeg,0.0,cosv=True),
    vLength = Offset("Galactic",0.0,scanheight*amintodeg,cosv=True),
    vDelta  = Offset("Galactic",0.0,vdelta,cosv=True),
    scanDuration = scanwidth/scanrate * 60,
    beamName="1")

Slew("LimaBeanOff")
Track("LimaBeanOff",None,60)
