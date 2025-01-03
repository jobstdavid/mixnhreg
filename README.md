
# Mixtures of Non-Homogeneous Linear Regression Models

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/mixnhreg)](https://CRAN.R-project.org/package=mixnhreg)
[![R-CMD-check](https://github.com/jobstdavid/mixnhreg/workflows/R-CMD-check/badge.svg)](https://github.com/jobstdavid/mixnhreg/actions)
[![version](https://img.shields.io/badge/version-0.1.0-green.svg?style=flat)](https://github.com/jobstdavid/mixnhreg)

<!-- badges: end -->

## Overview

The R package `mixnhreg` allows to estimate mixtures of non-homogenous
linear regression models. Linear predictors are separately employed for
the location, scale and weight parameter of each mixture component
allowing for heteroscedastic mixture regression models.

The infrastructure offers:

- predefined mixture distributions, such as e.g. mixture normal
  distribution.
- the possibility for the user to extending the zoo of usable mixture
  distributions by its own.

The model estimation supports:

- numerical optimization algorithms such as e.g. BFGS or Nelder-Mead.
- non-cyclic gradient boosting and therefore a variable selection
  procedure.
- different loss functions for the model estimation, i.e. negative
  log-likelihood or continuous ranked probability score.

## Installation

You can install the development version from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("jobstdavid/mixnhreg")
```

## Example

### Model estimation

``` r
library(mixnhreg)

# load data
data("station")

# fit mixture normal distribution with two components via BFGS
(fit_optim <- mixnhreg(formula = obs ~ sin1 + cos1 + temp_mean | temp_ctrl,
                       scale.formula = ~ sin1 + cos1 + log(temp_sd) | 1,
                       weight.formula = ~ sin1 + cos1 | 1,
                       data = station,
                       control = control_optim()))
#> Family: Mixture Normal Distribution with 2 components
#> CRPS: 1654.3006, logLik: -3353.8597, AIC: 6737.7195, BIC: 6820.3266

# fit mixture normal distribution with two components via gradient-boosting
(fit_boost <- mixnhreg(formula = obs ~ sin1 + cos1 + temp_mean | temp_ctrl,
                       scale.formula = ~ sin1 + cos1 + log(temp_sd) | 1,
                       weight.formula = ~ sin1 + cos1 | 1,
                       data = station,
                       control = control_boost(mstop = "cv")))
#> Family: Mixture Normal Distribution with 2 components
#> CRPS: 1707.8399, logLik: -3405.8353, AIC: 6837.6705, BIC: 6909.2634 
#> Boosting iterations: 5998, stopping criterion: cv
```

### Model prediction

``` r
# parameter predictions
par_optim <- predict(fit_optim, type = "parameter")
par_boost <- predict(fit_boost, type = "parameter")

# get observations and mean forecasts
obs <- na.omit(station$obs)
mean_optim <- rowSums(par_optim$location * par_optim$weight)
mean_boost <- rowSums(par_boost$location * par_boost$weight)
```

### Residual plot

``` r
plot(obs - mean_optim, 
     xlab = "Index", 
     ylab = "Observation - Mean", 
     pch = 19, 
     col = adjustcolor("red", alpha = 0.5))
points(obs - mean_boost, 
       pch = 19, 
       col = adjustcolor("steelblue", alpha = 0.25))
legend("bottomright", 
       legend = c("optim", "boost"),
       pch = 19,
       col = c("red", "steelblue"))
grid()
```

<img src="man/figures/README-example3-1.png" width="100%" style="display: block; margin: auto;" />

## Contact

Feel free to contact <jobstd@uni-hildesheim.de> if you have any
questions or suggestions.

## References

- Grimit E.P., Gneiting T., Berrocal V., Johnson N.A. (2006). The
  continuous ranked probability score for circular variables and its
  application to mesoscale forecast ensemble verification. Quarterly
  Journal of the Royal Meteorological Society, 132, pp. 2925–2942. doi:
  <https://doi.org/10.1256/qj.05.235>.
- Hepp, T., J. Zierk, and E. Bergherr (2023). “Component-wise boosting
  for mixture distributional regression models”. In: Proceedings of the
  37th International Workshop on Statistical Modelling. url:
  <https://iwsm2023.statistik.tu-dortmund.de/>.
- Jobst, D. (2024). Gradient-Boosted Mixture Regression Models for
  Postprocessing Ensemble Weather Forecasts. doi:
  <https://doi.org/10.48550/arXiv.2412.09583>.
- Messner, J. W., G. J. Mayr, and A. Zeileis (2017). “Nonhomogeneous
  Boosting for Predictor Selection in Ensemble Postprocessing”. In:
  Monthly Weather Review 145.1, pp. 137–147. doi:
  <https://doi.org/10.1175/mwr-d-16-0088.1>.
- Thomas, J. et al. (2017). “Gradient boosting for distributional
  regression: faster tuning and improved variable selection via
  noncyclical updates”. In: Statistics and Computing 28.3, pp. 673–687.
  doi: <https://doi.org/10.1007/s11222-017-9754-6>.
