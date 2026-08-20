# RtForecastR

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/RtForecastR)](https://CRAN.R-project.org/package=RtForecastR)
[![R-CMD-check](https://github.com/rajsubediresearch/RtForecastR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rajsubediresearch/RtForecastR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

A lightweight R package for real-time effective reproduction number
(R_t) estimation and one-step-ahead forecasting from case-count time
series. Disease-agnostic: works for any pathogen with a known
generation interval.

## Installation

```r
install.packages("RtForecastR")
```

The development version can be installed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("rajsubediresearch/RtForecastR")
```

## Usage

```r
library(RtForecastR)

data(measles_cdmx)
fit <- rt_forecast(measles_cdmx$time, measles_cdmx$cases,
                    mean_GI = 11/7, var_GI = (4/7)^2)
fit
plot(fit, which = "Rt")
plot(fit, which = "forecast")
```

See `vignette("rtforecastr-walkthrough")` for the full walkthrough,
and `?rt_forecast`, `?wis`, `?score_batches` for calibration/forecast
scoring (mean absolute error, mean squared error, root mean squared
error, empirical coverage, weighted interval score).

## What it does

- Filtered (real-time/causal) and smoothed (retrospective) R_t
  estimates, each with 95% credible intervals
- One-step-ahead predicted case counts (in-sample), for checking model
  adequacy
- A genuine out-of-sample forecast one time-step beyond your most
  recent data point, with 50%/95% predictive intervals
- Elimination probability, P(R_t < 1), at every time point
- Forecast calibration metrics: mean absolute error, mean squared
  error, root mean squared error, empirical coverage, and the
  weighted interval score (WIS)

## Attribution

This package builds directly on existing, published methods. Please
cite accordingly if you use it:

**Core R_t estimation algorithm** - `R/epiFilter.R`, `R/epiSmoother.R`,
and the base logic in `R/recursPredict.R` are ports of the EpiFilter
algorithm:
> Parag KV. (2021) "Improved estimation of time-varying reproduction
> numbers at low case incidence and between epidemic waves." *PLOS
> Computational Biology* 17(9): e1009347.
> <https://github.com/kpzoo/EpiFilter>

**Generation-interval discretization** - `compute_lambda()` follows
the convention used in `get_Rt.m`, part of the GrowthPredict toolbox:
> Chowell G, Bleichrodt A, Dahal S, Tariq A, Roosa K, Hyman JM, Luo R.
> (2024) "GrowthPredict: A toolbox and tutorial-based primer for
> fitting and forecasting growth trajectories using phenomenological
> growth models." *Scientific Reports* 14, 1630.
> <https://github.com/gchowell/forecasting_growthmodels>

**Weighted interval score** - `wis()` implements:
> Bracher J, Ray EL, Gneiting T, Reich NG. (2021) "Evaluating epidemic
> forecasts in an interval format." *PLOS Computational Biology* 17(2):
> e1008618.

**Everything else** - the R package interface, the configurable
prediction grid, the out-of-sample forecast extension, the multi-level
quantile prediction (for WIS), and the calibration/scoring functions -
is original work by Raj Subedi, built on top of the above.

## History

This package supersedes an earlier script-based version of the same
tool (`run_rtforecast.R` + `core/*.R`). That version's archived weekly
output is preserved on the `pre-package-v0` tag for reference.

## License

GPL-3. Copyright holders and their specific contributions are declared
in `Authors@R` (see `DESCRIPTION`) - Kris V. Parag holds copyright on
the unmodified/lightly modified EpiFilter files (`R/epiFilter.R`,
`R/epiSmoother.R`, `R/recursPredict.R`'s original logic); Raj Subedi
holds copyright on everything else.
