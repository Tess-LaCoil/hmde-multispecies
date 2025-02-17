# hmde-multispecies
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

This code implements an analysis of growth behaviour across 17 tropical tree species from Barro Colorado Island.

## Methods
We use the statistical method of a hierarchical Bayesian longitudinal model as 
introduced by O'Brien et al. (2024), implemented through the [hmde](https://github.com/traitecoevo/hmde) package. 

The growth function we fit was first used in Canham et al. (2004), and is given by

$$g \left( S \left( t \right), g_{max}, S_{max}, k \right) = \frac{dS}{dt} = g_{max} \exp \Bigg( -\frac{1}{2} \bigg( \frac{ \ln \left( S \left( t \right) / S_{max} \right) }{k} \bigg)^2 \Bigg), $$ 

where $g_{max}$ is the maximum growth rate, $S_{max}$ is the size at which that maximum occurs, and $k$ controls how narrow or wide the peak is.

The hierarchical model estimates $g_{max}$, $S_{max}$ and $k$ for each individual, 
with a species-level distribution sitting over the top of the individual level.
We implement log-normal priors at the species level for each individual-level parameter. 
Our error model is
$$s_{ij} = S_i(t_j) + \mathcal{N}(0, \sigma_e),$$
where we fit $\sigma_e$ as a half-Cauchy distributed error parameter. 

### Data
We take data from the Barro Colorado Island (BCI) long term forest plot (Condit et al. 2019),
from 17 species of tropical trees. From those 17 species we select individuals with
6 observations (all censuses since a change of precision in 1990). Our intended sample 
size is 300, but for species with fewer than 300 individuals we take the entire sample.
For species with more than 300 individuals with 6 observations, we take a simple random
sample without replacement of size 300.

The table below gives the species names, species code from the BCI data, sample size, 
and runtime for the hierarchical model fit to that species.

|Species | Species code | Sample size | Walltime (hrs)|
| :---------- | :----------: | :------: | :------: |
| *Alseis blackiana* | alsebl | 300 | 153 |
| *Beilschmiedia tovarensis* | beilpe | 300 | 147 |
| *Cordia bicolor* | cordbi | 300 | 28 |
| *Faramea occidentalis* | faraoc | 300 | 21 |
| *Garcinia recondita* | gar2in | 300 | 32 |
| *Hirtella triandra* | hirttr | 300 | 36 |
| *Jacaranda copaia* | jac1co | 127 | 4 |
| *Prioria copaifera* | pri2co | 300 | 41 |
| *Protium panamense* | protpa | 300 | 75 |
| *Protium tenuifolium* | protte | 300 | 108 |
| *Quararibea asterolepis* | quaras | 300 | 50 |
| *Swartzia simplex var. 1* | swars1 | 300 | 143 |
| *Swartzia simplex var. 2* | swars2 | 300 | 185 |
| *Simarouba amara* | simaam | 234 | 26 |
| *Tachigali panamensis* | tachve | 300 | 64 |
| *Protium stevensonii* | tet2pa | 300 | 84 |
| *Trichilia tuberculata* | tri2tu | 300 | 135 |

## References
Canham, C. D., Papaik, M. J., Uriarte, M., McWilliams, W. H., Jenkins, J. C., & Twery, M. J. (2006). Neighborhood analyses of canopy tree competition along environmental gradients in New England forests. Ecological applications, 16(2), 540-554.

Condit, R., Pérez, R., Aguilar, S., Lao, S., Foster, R., & Hubbell, S. (2019). Complete data from the Barro Colorado 50-ha plot: 423617 trees, 35 years. https://datadryad.org/stash/dataset/doi:10.15146/5xcp-0d46

O'Brien, T., Warton, D., & Falster, D. (2024). Yes, they're all individuals: Hierarchical models for repeat survey data improve estimates of tree growth and size. Methods in Ecology and Evolution, 16(1), 183-196.

O’Brien, T. A., Kar, F., Warton, D., and Falster, D. S. (2025). hmde: Hierarchical methods for differential equations. bioRxiv.
