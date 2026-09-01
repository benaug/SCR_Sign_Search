# SCR_Sign_Search
An SCR model for sign searches (area or transect) where continuous sign locations are recorded.

This is an SCR model for sign searches where:
1) search effort varies across the landscape (and is recorded)
2) the continuous location of a sign is recorded when they are detected
3) an individuals sign can be detected in multiple locations per occasion (the model currently only considers one occasion)

The status quo* for sign searches is to discretize the effort into "effective detectors", and snap the continuous
sign locations to the centroid of the effective detector in which they were found. When doing this, if the effective detectors
are spaced too far apart relative to sigma, the sigma estimates are positively biased. The model here allows the 
continuous locations to be used so this does not occur.

*Efford has a similar model in secr and Zhang et al. is similar, too. The main difference I see is that my approach
includes within home range resource selection and capitalizes on the separability of the x and y dimensions
with the bivariate normal distribution to reduce computations, store reusable quantities, and speed up the MCMC. However,
this only works for a BVN availability distribution.
https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.3887

This model is modified from the Royle-Young model for area searches of *individuals* found here:
https://github.com/benaug/RoyleYoung

The data augmentation approach is explained here:
https://github.com/benaug/SCR-N-Prior-Data-Augmentation

One difference, though, is that in this model, I am gating the activity center and within home range space use model by z, so that no
likelihood is contributed when z_i=0 and activity centers are jointly turned on/off with z_i. When z_i=0, I set the s values and 
within home range space use objects to 0. This makes the MCMC more efficient.

In this model, individual expected detection rate at a trap is assumed to be proportional to the individual use probability of the cell containing trap.
(perfect "compensatory heterogeneity").

lam[i,j] <- lambda.detect[j]*use.dist[i,trap.to.cell[j]] #trap.to.cell maps traps to cells

where lambda.detect[j] is a function of (log-transformed) effort

log(lambda.detect[j]) <- beta0.lam + beta1.lam*E[j]

The within home range space use model is the same as the between primary period relocation model used in the Jolly-Seber N-Prior Data Augmentation
repository. It is an RSF weighted BVN movement model. The factored representation used to fit the model here is the same as used in 
Jolly-Seber and is described in the movement model section of the readme for that repository:
https://github.com/benaug/Jolly-Seber-N-Prior-DA

There are two sets of model files. Files with "1K" are set up for 1 occasion. Files with "multiK" are set up for 2+ occasions.


The model is also set up for supplemental telemetry data (not the same individuals being detected). This can be turned off,
or telemetry from the individuals being detected can be included with some modification. I expect telemetry will often be required
to estimate within home range space use with precision.

The test script does not simulate a transect search to define traps. I just allocate some equally spaced cells to search. This is 
more efficient than a real transect search for the same amount of effort. Simulate transect searches for more realistic scenarios.