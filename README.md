# SCR Sign Search

Spatial capture-recapture models for transect or area searches in which individually identifiable sign, such as scat, is detected and the continuous location of each sign is recorded.

The models allow:

1. search effort to vary across the landscape and among sampling occasions
2. the exact continuous location of each detected sign to be retained rather than snapping detections to detector centroids
3. multiple signs from the same individual to be detected at multiple locations during the same sampling occasion
4. spatial density covariates
5. within home range resource selection

Files with `1K` are set up for a single sampling occasion. Files with `multiK` are set up for two or more occasions.

## Related models

A common approach for sign searches is to discretize the searched area into effective detectors and assign each sign to the centroid of the detector in which it was found. When these spatial units are coarse relative to the scale of space use, replacing the observed sign coordinates with detector centroids can bias estimation of the spatial scale parameter. The models here instead use the exact continuous sign locations while still representing search effort and habitat covariates on a spatial grid.

Related approaches include the SCR model for area searches of Efford (2011), which is implemented for polygon and transect searches in `secr`, and the continuous space point process SCR framework of Zhang et al. (2023). The implementation in this repository combines a continuous sign location likelihood with within home range resource selection and a factored representation of the bivariate normal availability distribution that reduces repeated calculations during MCMC.

The model is also closely related to the Royle-Young model for transect or area searches of individuals:

https://github.com/benaug/RoyleYoung

The main observation model distinction is that the Royle-Young model represents one realized animal location on each sampling occasion, so an individual can occur in only one location on that occasion. In the sign search model, an individual can produce multiple detected signs during the same occasion and those signs can occur in multiple searched cells. The same underlying within home range space use distribution is used to allocate sign intensity across space.

The N-prior data augmentation approach used here is described in:

https://github.com/benaug/SCR-N-Prior-Data-Augmentation

## Within home range space use

Let $\mathbf{s}_i$ denote the activity center of individual $i$. For grid cell $c$, define the resource selection weight

$$
r_c = \exp(\mathbf{x}_c^\top\boldsymbol{\beta}^{\mathrm{RSF}}).
$$

The availability of cell $c$ under the bivariate normal home range kernel is

$$
a_{ic} = P(\mathbf{u} \in c \mid \mathbf{s}_i,\sigma),
$$

where $\sigma$ controls the spatial scale of within home range use and $\mathbf{u}$ denotes a possible continuous location.

Resource selection reweights this availability distribution. Define

$$
D_i = \sum_{h=1}^{C} r_h a_{ih}.
$$

The probability that space use by individual $i$ occurs in cell $c$ is then

$$
q_{ic} =
\frac{r_c a_{ic}}
{D_i}.
$$

Therefore, proximity to the activity center determines availability through the bivariate normal kernel, while habitat covariates modify use through the cell-specific RSF weights.

This is the same within home range space use model described in more detail in the Royle-Young repository.

It is also the same normalized RSF-weighted bivariate normal structure used for between primary period activity center relocation in the Jolly-Seber N-Prior Data Augmentation repository:

https://github.com/benaug/Jolly-Seber-N-Prior-DA

## Sign detection model

Let $y_{ikj}$ be the number of signs from individual $i$ detected during sampling occasion $k$ in search unit $j$, and let $c(j)$ denote the grid cell containing that search unit. The expected number of detected signs is proportional to the individual's use probability for that cell,

$$
y_{ikj} \sim \mathrm{Poisson}(\mu_{ikj}),
$$

with

$$
\mu_{ikj} = \lambda_{kj} q_{i,c(j)}.
$$

The search specific detection intensity is modeled as

$$
\lambda_{kj} = S_{kj}\exp(\beta_{0,k}^{\lambda}+\beta_1^{\lambda}E_{kj}),
$$

where $S_{kj}$ indicates whether search unit $j$ was surveyed on occasion $k$ and $E_{kj}$ is the effort covariate. In the supplied examples, effort is log transformed before entering the model. The `multiK` implementation allows the detection intercept $\beta_{0,k}^{\lambda}$ to vary among occasions. The `1K` model is the corresponding single occasion version.

This construction assumes that sign intensity in a searched cell is proportional to the individual's use probability of that cell. In the terminology of Efford and Mowat (2014), this corresponds to a compensatory relationship between the spatial scale of use and local encounter intensity. If search intensity is constant and the entire state space is searched, changing the spatial distribution of use redistributes where an individual's signs are expected without changing the individual's total expected number of signs. With incomplete spatial coverage or spatially varying effort, expected detections additionally depend on how the individual's use distribution overlaps the searched area and effort surface.

Unlike the Royle-Young observation model, this is a Poisson sign process. An individual can therefore contribute zero, one, or multiple detected signs during an occasion, including signs in multiple cells.

## Continuous sign locations

The count model above determines how many signs are detected in each searched cell. The exact location of each detected sign is then modeled continuously within its observed cell.

For an observed sign location $\mathbf{u} _{ik\ell}$ in cell $c$,
the conditional within-cell density is

$$
f(\mathbf{u} _{ik\ell}\mid\mathbf{u} _{ik\ell}\in c,\mathbf{s}_i,\sigma) =
\frac{
f_{\mathrm{BVN}}(\mathbf{u} _{ik\ell}\mid\mathbf{s}_i,\sigma)
}{
a_{ic}
}.
$$

The RSF weight does not appear in this conditional density because it is constant within a grid cell. Combining the cell use probability with the conditional within-cell density gives

$$
q_{ic}
f(\mathbf{u} _{ik\ell}\mid\mathbf{u} _{ik\ell}\in c,\mathbf{s}_i,\sigma) =
\frac{
r_c f_{\mathrm{BVN}}(\mathbf{u} _{ik\ell}\mid\mathbf{s}_i,\sigma)
}{
D_i
}.
$$

Consequently, the count plus location formulation can also be viewed as a continuous Poisson point process. Within searched cell $c(j)$, the point intensity for signs from individual $i$ on occasion $k$ is

$$
\rho_{ikj}(\mathbf{u}) =
\lambda_{kj}
\frac{
r_{c(j)} f_{\mathrm{BVN}}(\mathbf{u}\mid\mathbf{s}_i,\sigma)
}{
D_i
},
\qquad \mathbf{u}\in c(j).
$$

Integrating this intensity over cell $c(j)$ gives $\lambda_{kj}q_{i,c(j)}$, which is the mean of the Poisson count model above. Thus, the implementation uses the grid to represent search effort and cell level resource selection while retaining the exact continuous coordinates of detected signs.

## Factored representation of space use

For the isotropic bivariate normal kernel used here, the $x$ and $y$ dimensions are independent. On a rectangular grid, the bivariate normal probability mass of cell $c$ therefore factors exactly. If cell $c$ corresponds to x-cell $j(c)$ and y-cell $k(c)$,

$$
a_{ic} =
a^x_{i,j(c)}
a^y_{i,k(c)}.
$$

The normalizing constant can then be calculated as

$$
D_i =
\sum_{c=1}^{C}
r_c
a^x_{i,j(c)}
a^y_{i,k(c)}.
$$

The fitted models store only the two one-dimensional availability vectors, `avail.x` and `avail.y`, together with the RSF normalizing constant `use.denom`. When an activity center or $\sigma$ changes, the one-dimensional availability distributions and normalizing constant are updated. When the RSF coefficient changes, the availability distributions do not change, so only the normalizing constant needs to be recomputed.

This avoids storing and repeatedly recalculating a full two-dimensional availability distribution for every individual. The factorization is exact for the independent coordinate bivariate normal model on the rectangular grid used here. The specific implementation would need to be modified for a nonseparable availability distribution.

A more detailed description of the same factorization is provided in the Royle-Young and Jolly-Seber repositories linked above.

## Gating activity centers and space use by population membership

The activity center and within home range space use components are gated by the data augmentation indicator $z_i$. When $z_i=0$, the individual's activity center, one-dimensional availability vectors, and use normalizing constant are set to zero and contribute no likelihood. When an augmented individual is turned on, its activity center and associated space use quantities are proposed jointly with $z_i$.

This avoids maintaining and updating latent activity centers and space use distributions for augmented individuals that are not currently part of the population. In these models this reduces unnecessary calculations and improves MCMC efficiency.

## Model versions

There are two sets of model files:

1. Files containing `1K` are set up for one sampling occasion.
2. Files containing `multiK` are set up for two or more sampling occasions.

The model for multiple occasions retains a single activity center and within home range space use distribution for each individual while allowing the search effort and detection intercept to vary among occasions.

## Simulation note

The supplied test scripts do not simulate a literal transect layout. Instead, a set of regularly spaced grid cells is designated as searched. This produces a relatively efficient spatial sampling design for a given amount of effort. For evaluating performance under realistic transect or area search designs, the search geometry should be simulated accordingly.

## References

Efford, M. G. 2011. Estimation of population density by spatially explicit capture-recapture analysis of data from area searches. *Ecology* 92:2202-2207.

https://doi.org/10.1890/11-0332.1

Efford, M. G., and G. Mowat. 2014. Compensatory heterogeneity in spatially explicit capture-recapture data. *Ecology* 95:1341-1348.

https://doi.org/10.1890/13-1497.1

Royle, J. A., and K. V. Young. 2008. A hierarchical model for spatial capture-recapture data. *Ecology* 89:2281-2289.

https://doi.org/10.1890/07-0601.1

Zhang, W., J. D. Chipperfield, J. B. Illian, P. Dupont, C. Milleret, P. de Valpine, and R. Bischof. 2023. A flexible and efficient Bayesian implementation of point process models for spatial capture-recapture data. *Ecology* 104:e3887.

https://doi.org/10.1002/ecy.3887