# SCR_Sign_Search
An SCR model for sign searches (area or transect) where continuous sign locations are recorded.

This is an SCR model for sign searches where:
1) search effort varies across the landscape (and is recorded)
2) the continuous location of a sign is recorded when they are detected
3) an individuals sign can be detected in multiple locations per occasion (the model currently only considers one occasion)

The status quo for sign searches is to discretize the effort into "effective detectors", and snap the continuous
sign locations to the centroid of the effective detector in which they were found. When doing this, if the effective detectors
are spaced too far apart relative to sigma, the sigma estimates are positively biased. The model here allows the 
continuous locations to be used so this does not occur.

This model is modified from the Royle-Young model for area searches found here:
https://github.com/benaug/RoyleYoung

Within home range resource selection can be modeled. This feature is currently turned off, but instructions in test script to
turn it back on. 

Individual expected detection rate at a trap is assumed to be proportional to the individual use probability of the cell containing trap.
(perfect "compensatory heterogeneity").
lam[i,j] <- lambda.detect[j]*use.dist[i,trap.to.cell[j]] #trap.to.cell maps traps to cells

The model is also set up for supplemental telemetry data (not the same individuals beign detected). This can be turned off,
or telemetry from the individuals being detected can be included.

I have not tested this code, yet, but will get to it. It appears to be working correctly and is built on the 
well tested Royle-Young code. There may be more efficiency gains to be made, will look at that.