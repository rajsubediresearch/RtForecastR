#' Estimate and forecast the time-varying effective reproduction number
#'
#' Fits filtered (real-time/causal) and smoothed (retrospective) R_t
#' estimates to a case-count time series using the EpiFilter algorithm
#' (Parag 2021), computes one-step-ahead in-sample predictions (a model
#' adequacy check), a genuine out-of-sample one-step forecast with
#' predictive intervals, and the elimination probability P(R_t < 1) at
#' every time point. Disease-agnostic: supply the generation interval
#' for your own pathogen via `mean_GI`/`var_GI`.
#'
#' Unlike a script, this returns an object rather than writing files -
#' use [plot.rtforecast()] if you want the plots `run_rtforecast.R`
#' used to write directly.
#'
#' @param time numeric vector of time indices (e.g. epi weeks)
#' @param cases numeric vector of case counts, same length as `time`
#' @param mean_GI mean generation interval, in the same time units as
#'   `time` (e.g. weeks if your data is weekly)
#' @param var_GI variance of the generation interval, same time units
#' @param Rmin,Rmax bounds of the grid searched over for R_t (default 0.01, 10)
#' @param grid_size number of points in the R_t grid (default 200)
#' @param eta diffusion (state) noise controlling smoothness of R_t (default 0.1)
#' @param ci_level tail probability defining the reported main credible
#'   interval (default 0.025, i.e. a 95% interval)
#' @param quantile_levels additional quantile levels (each in (0, 0.5))
#'   to compute for the out-of-sample forecast, used by [wis()]. Default
#'   adds the levels a standard 11-interval WIS wants; set to `NULL` to
#'   skip and only compute the 50%/95% forecast interval.
#' @param maxI upper bound of the internal prediction grid; if `NULL`
#'   (default) it's set automatically to `3 * max(cases)` (minimum 2000)
#' @return an object of class `"rtforecast"`: a list with elements
#'   `results` (per-time R_t estimates and elimination probability),
#'   `predictions` (in-sample one-step-ahead predictions), `forecast`
#'   (the one-step-ahead out-of-sample forecast, with a `quantiles`
#'   element if `quantile_levels` was set), and the raw filter/smoother
#'   objects for advanced use.
#' @examples
#' data(measles_cdmx)
#' fit <- rt_forecast(measles_cdmx$time, measles_cdmx$cases,
#'                     mean_GI = 11/7, var_GI = (4/7)^2)
#' head(fit$results)
#' fit$forecast
#' @export
rt_forecast <- function(time, cases, mean_GI, var_GI,
                         Rmin = 0.01, Rmax = 10, grid_size = 200, eta = 0.1,
                         ci_level = 0.025,
                         quantile_levels = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45),
                         maxI = NULL){

  if (length(time) != length(cases)) stop("time and cases must be the same length")
  n <- length(cases)
  if (n < 3) stop("need at least 3 time points to fit R_t")
  if (is.null(maxI)) maxI <- max(50, round(3 * max(cases, na.rm = TRUE)))

  Lam <- compute_lambda(cases, mean_GI, var_GI)

  m <- grid_size
  pR0 <- (1/m) * rep(1, m)
  Rgrid <- seq(Rmin, Rmax, length.out = m)
  a <- ci_level
  tvec <- 2:n

  Rfilt <- epi_filter(Rgrid, m, eta, pR0, length(tvec), Lam[tvec], cases[tvec], a)
  Ifilt <- recurs_predict(Rgrid, Rfilt$pR, Lam[tvec], Rfilt$Rmean, a, maxI = maxI)
  Rsmooth <- epi_smoother(Rgrid, m, Rfilt$pR, Rfilt$pRup, length(tvec), Rfilt$pstate, a)
  Ismooth <- recurs_predict(Rgrid, Rsmooth$qR, Lam[tvec], Rsmooth$Rmean, a, maxI = maxI)

  # Genuine out-of-sample forecast: one step beyond the data
  cases_extended <- c(cases, NA)
  Lam_next <- compute_lambda(cases_extended, mean_GI, var_GI)
  Lam_forecast <- utils::tail(Lam_next, 1)
  pR_last <- Rfilt$pR[nrow(Rfilt$pR), ]

  Igrid <- 0:maxI
  rate <- Lam_forecast * Rgrid
  pI <- sapply(Igrid, function(k) sum(stats::dpois(k, rate) * pR_last))
  Fpred <- cumsum(pI) / sum(pI)

  forecast_mean <- Lam_forecast * Rfilt$Rmean[length(Rfilt$Rmean)]
  forecast_lo95 <- Igrid[which(Fpred >= a)[1]]
  forecast_hi95 <- Igrid[which(Fpred >= 1 - a)[1]]
  forecast_lo50 <- Igrid[which(Fpred >= 0.25)[1]]
  forecast_hi50 <- Igrid[which(Fpred >= 0.75)[1]]
  forecast_time <- time[n] + 1

  forecast <- data.frame(
    time = forecast_time, pred_mean = round(forecast_mean, 1),
    lo95 = forecast_lo95, hi95 = forecast_hi95,
    lo50 = forecast_lo50, hi50 = forecast_hi50,
    R_used = round(Rfilt$Rmean[length(Rfilt$Rmean)], 3)
  )

  forecast_quantiles <- NULL
  if (!is.null(quantile_levels) && length(quantile_levels) > 0){
    # Derived directly from Fpred/Igrid above - the exact same forecast
    # distribution that produced forecast_mean/lo95/hi95/lo50/hi50, so
    # forecast_quantiles[["0.025"]] is guaranteed identical to
    # c(forecast$lo95, forecast$hi95), etc. Previously this used a
    # separate in-sample calculation that did not correspond to the
    # true one-step-ahead forecast point; fixed in 0.1.1.
    all_levels <- sort(unique(c(a, 0.25, quantile_levels)))
    forecast_quantiles <- lapply(all_levels, function(lev){
      lo <- Igrid[which(Fpred >= lev)[1]]
      hi <- Igrid[which(Fpred >= 1 - lev)[1]]
      c(lo, hi)
    })
    names(forecast_quantiles) <- as.character(all_levels)
  }

  # Elimination probability: P(R_t < 1) at each time point
  below1 <- Rgrid < 1
  Pelim_filt <- apply(Rfilt$pR, 1, function(p) sum(p[below1]))
  Pelim_smooth <- apply(Rsmooth$qR, 1, function(p) sum(p[below1]))

  results <- data.frame(
    time          = time[tvec],
    cases         = cases[tvec],
    R_filtered    = round(Rfilt$Rmean, 3),
    R_filt_lo95   = round(Rfilt$Rhat[1, ], 3),
    R_filt_hi95   = round(Rfilt$Rhat[2, ], 3),
    R_smoothed    = round(Rsmooth$Rmean, 3),
    R_smooth_lo95 = round(Rsmooth$Rhat[1, ], 3),
    R_smooth_hi95 = round(Rsmooth$Rhat[2, ], 3),
    P_elim_filt   = round(Pelim_filt, 3),
    P_elim_smooth = round(Pelim_smooth, 3)
  )

  tvec_pred <- tvec[-1]
  predictions <- data.frame(
    time      = time[tvec_pred],
    cases     = cases[tvec_pred],
    pred_next = round(Ifilt$pred, 1),
    pred_lo95 = round(Ifilt$predInt[1, ], 1),
    pred_hi95 = round(Ifilt$predInt[2, ], 1),
    pred_lo50 = round(Ifilt$predInt[3, ], 1),
    pred_hi50 = round(Ifilt$predInt[4, ], 1)
  )

  out <- list(
    results = results,
    predictions = predictions,
    forecast = forecast,
    forecast_quantiles = forecast_quantiles,
    call_args = list(mean_GI = mean_GI, var_GI = var_GI, Rmin = Rmin, Rmax = Rmax,
                      grid_size = grid_size, eta = eta, ci_level = ci_level, maxI = maxI),
    raw = list(Rgrid = Rgrid, time = time, cases = cases, tvec = tvec, tvec_pred = tvec_pred,
               Rfilt = Rfilt, Rsmooth = Rsmooth, Ifilt = Ifilt, Ismooth = Ismooth)
  )
  class(out) <- "rtforecast"
  out
}

#' @export
print.rtforecast <- function(x, ...){
  cat("<rtforecast>", nrow(x$results), "time points\n")
  cat("Latest R_t (filtered):", x$results$R_filtered[nrow(x$results)],
      sprintf("(%.3f-%.3f)", x$results$R_filt_lo95[nrow(x$results)], x$results$R_filt_hi95[nrow(x$results)]), "\n")
  cat("One-step-ahead forecast for time", x$forecast$time, ":\n")
  cat(sprintf("  %.1f cases (95%% CI: %d-%d, 50%% CI: %d-%d)\n",
              x$forecast$pred_mean, x$forecast$lo95, x$forecast$hi95,
              x$forecast$lo50, x$forecast$hi50))
  invisible(x)
}
