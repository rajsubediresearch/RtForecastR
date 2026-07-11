# RtForecastR

A lightweight R tool for real-time effective reproduction number (R_t)
estimation and one-step-ahead forecasting from case-count time series.
Disease-agnostic: works for any pathogen with a known generation
interval, not just the measles example shipped here.

**Author:** Raj Subedi

## What it does

- Filtered (real-time/causal) and smoothed (retrospective) R_t
  estimates, each with 95% credible intervals
- One-step-ahead predicted case counts (in-sample), for checking model
  adequacy against what actually happened
- A genuine out-of-sample forecast one time-step beyond your most
  recent data point, with its own 50%/95% predictive interval
- Elimination probability, P(R_t < 1), at every time point
- Three report-ready plots per run, each generated with and without an
  optional "last observed value" reference line
- Output automatically organized by dataset and batch, so re-running
  weekly as new data arrives never overwrites previous results

## Attribution

This tool builds directly on existing, published methods rather than
introducing a new one. Please cite accordingly if you use it:

**Core R_t estimation algorithm** - `core/epiFilter.R`, `epiSmoother.R`,
and `recursPredict.R` are unmodified files from the EpiFilter package:

> Parag KV. (2021) "Improved estimation of time-varying reproduction
> numbers at low case incidence and between epidemic waves." *PLOS
> Computational Biology* 17(9): e1009347.
> https://github.com/kpzoo/EpiFilter

`core/recursPredictExt.R` is a modified version of `recursPredict.R`
(same source above): the original hardcodes its internal prediction
grid to 800 and throws `"Epidemic size too large"` on any outbreak
whose predicted incidence exceeds ~720. This version takes the grid's
upper bound as a parameter instead, computed automatically from your
data's peak incidence in `run_rtforecast.R`.

**Generation-interval discretization** - `core/computeLambda.R`
computes the renewal-equation "total infectiousness" term EpiFilter
needs, using a Gamma(mean, variance) discretization. This follows the
same convention used in `get_Rt.m`, part of the GrowthPredict toolbox:

> Chowell G, Bleichrodt A, Dahal S, Tariq A, Roosa K, Hyman JM, Luo R.
> (2024) "GrowthPredict: A toolbox and tutorial-based primer for
> fitting and forecasting growth trajectories using phenomenological
> growth models." *Scientific Reports* 14, 1630.
> https://github.com/gchowell/forecasting_growthmodels

**Everything else** (the out-of-sample forecast extension, the
three-plot reporting layout, dataset/batch output organization) is
original work by Raj Subedi, built on top of the above.

## Requirements

```r
install.packages("caTools")
```

## Usage

1. Prepare your data as a two-column, no-header text file: time index,
   case count (see `examples/measles_cdmx_example.txt` for the format).
2. In `run_rtforecast.R`, set:
   - `DATA_FILE` to your data's path
   - `mean_GI` / `var_GI` to your disease's generation interval (mean
     and SD, converted to variance) - the shipped values are a measles
     example, **not** a default assumption; replace them for your data
3. Run:
   ```
   Rscript run_rtforecast.R
   ```
   or source it in RStudio with the working directory set to the repo
   root.

## Output

Results land in `./output/<dataset_name>/asof_week<N>/`, where N is the
last observed time point in that run - so re-running with a new
week/period of data lands in its own folder automatically.

Each batch folder contains:

- `epifilter_1_default_Rt.png` - filtered vs smoothed R_t
- `epifilter_2_current_situation.png` / `_refline.png` - zoomed,
  recent-period view + forecast
- `epifilter_3_observed_vs_predicted.png` / `_refline.png` - full
  period, in-sample one-step-ahead predictions + forecast
- `epifilter_results.csv` - R_t (filtered + smoothed) with 95% CIs and
  elimination probability, per time point
- `epifilter_predictions.csv` - in-sample one-step-ahead predictions
- `epifilter_forecast_next_week.csv` - the out-of-sample forecast
  record, meant to be compared against the next period's actual count
  once available, for tracking real forecast performance over time

## Repository structure

```
RtForecastR/
├── README.md
├── LICENSE                          (inherited from EpiFilter)
├── run_rtforecast.R                 (main script)
├── examples/
│   ├── measles_cdmx_example.txt     (worked example dataset)
│   └── jalisco_example.txt          (larger-outbreak example dataset)
└── core/
    ├── epiFilter.R                  (Parag et al., unmodified)
    ├── epiSmoother.R                (Parag et al., unmodified)
    ├── recursPredict.R              (Parag et al., unmodified)
    ├── recursPredictExt.R           (Subedi, modified from recursPredict.R
    │                                  - configurable grid size)
    └── computeLambda.R              (Subedi, GI-convolution helper)
```

## A note on methodology

Every prediction in the observed-vs-predicted plots is a genuine
one-step-ahead prediction: the predicted value at time *t* is computed
using only the R_t posterior estimated from data through time *t-1*,
then propagated forward via the Poisson renewal equation. The model
never has access to time *t*'s own outcome when generating its
prediction for time *t* - this is what makes the observed-vs-predicted
comparison a meaningful check of forecast performance rather than a
retrospective curve fit. The out-of-sample forecast at the end of each
plot applies the identical logic one additional step past the edge of
the observed data.

The first time point in any series has no prediction, since a
one-step-ahead prediction requires a prior R_t estimate that does not
yet exist at the start of a series (marked with a hollow circle in the
observed-vs-predicted plot).
