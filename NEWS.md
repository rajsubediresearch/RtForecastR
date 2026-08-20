# RtForecastR 0.1.1

* Bug fix: `rt_forecast()`'s `forecast_quantiles` element (used for
  weighted interval score calculations via `wis()`) was computed from
  a disconnected in-sample calculation that did not correspond to the
  genuine out-of-sample forecast reported in `fit$forecast`. It now
  derives directly from the same forecast distribution, so e.g.
  `forecast_quantiles[["0.025"]]` is guaranteed identical to
  `c(fit$forecast$lo95, fit$forecast$hi95)`. Also fixes a duplicate
  `"0.25"` entry that could appear in `forecast_quantiles` when the
  default `quantile_levels` was used.

# RtForecastR 0.1.0

* Initial CRAN release.
* `rt_forecast()`: filtered (real-time) and smoothed (retrospective)
  R_t estimation, one-step-ahead in-sample predictions, a genuine
  out-of-sample forecast, and elimination probability P(R_t < 1).
* `plot.rtforecast()`: R_t, forecast, and observed-vs-predicted plots.
* Forecast calibration metrics: `mae()`, `mse()`, `rmse()`,
  `coverage()`, `interval_score()`, `wis()`.
* `score_batches()`: score archived weekly forecast output against
  realized case counts.
* Two bundled example datasets: `measles_cdmx`, `jalisco`.
