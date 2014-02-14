"""
Create a 12x15' map of the Lima Bean 

Per Glen Langston's recommendations, I should:
    1. not use refs for each scan (because the pipeline doesn't use them)
    2. use "Track" to observe offs
    3. do pointing scans as normal; break up map into 2 sub-maps

Following pipeline recommendations:
https://safe.nrao.edu/wiki/bin/view/Kbandfpa/ObserverGuide?sortcol=table;up=#Reduction_Execute_Pipeline_with
"""


cat = Catalog("/users/aginsbur/GBT12B-221/limabean.astrid")

Configure("/users/aginsbur/GBT12B-221/H2CO_6cm_CSetup.py")

Slew("LimaBean")
AutoPeakFocus( frequency=4829., beamName="1" )
Break("Check pointing/focus")
Configure("/users/aginsbur/GBT12B-221/H2CO_6cm_CSetup.py")

Slew("LimaBean")
Balance()

amintodeg = 1/60.
scanrate = 3. # arcmin/min

# horizontal scans
RALongMap('LimaBean',     #center of map
    hLength = Offset("Galactic",15*amintodeg,0.0,cosv=True), 
    vLength = Offset("Galactic",0.0,12*amintodeg,cosv=True), 
    vDelta  = Offset("Galactic",0.0,1.*amintodeg,cosv=True), 
    scanDuration = 15/scanrate * 60,
    beamName="1")

