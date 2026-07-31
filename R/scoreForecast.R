######################################################################
## RtForecastR - forecast calibration / performance metrics
## Author: Raj Subedi
##
## Standard forecast-evaluation metrics, implemented from their usual
## definitions (not derived from any other package's code):
##   - MAE, MSE, RMSE: standard point-forecast accuracy metrics.
##   - Empirical coverage: standard predictive-interval calibration
##     check (does the realized value fall in the stated interval the
##     stated fraction of the time?).
##   - Weighted interval score (WIS): as defined in
##     Bracher J, Ray EL, Gneiting T, Reich NG. (2021) "Evaluating
##     epidemic forecasts in an interval format." PLOS Computational
##     Biology 17(2): e1008618. WIS is the interval-format scoring rule
##     used by the CDC/ECDC forecast hubs.
######################################################################

# ---- Point-forecast accuracy -----------------------------------------

#' Mean absolute error
#'
#' @param observed numeric vector of realized values
#' @param predicted numeric vector of point predictions (same length/order)
#' @return a single numeric value
#' @examples
#' mae(c(10, 20, 30), c(12, 18, 33))
#' @export
mae <- function(observed, predicted){
  if (length(observed) != length(predicted)) stop("observed and predicted must be the same length")
  mean(abs(observed - predicted), na.rm = TRUE)
}

#' Mean squared error
#'
#' @inheritParams mae
#' @return a single numeric value
#' @examples
#' mse(c(10, 20, 30), c(12, 18, 33))
#' @export
mse <- function(observed, predicted){
  if (length(observed) != length(predicted)) stop("observed and predicted must be the same length")
  mean((observed - predicted)^2, na.rm = TRUE)
}

#' Root mean squared error
#'
#' @inheritParams mae
#' @return a single numeric value
#' @examples
#' rmse(c(10, 20, 30), c(12, 18, 33))
#' @export
rmse <- function(observed, predicted){
  sqrt(mse(observed, predicted))
}

# ---- Interval calibration ---------------------------------------------

#' Empirical coverage of a prediction interval
#'
#' @param observed numeric vector of realized values
#' @param lo,hi numeric vectors giving the lower/upper interval bounds
#' @return the fraction of observations falling within \verb{[lo, hi]}
#' @examples
#' coverage(c(5, 15, 25), c(0, 10, 20), c(10, 20, 30))
#' @export
coverage <- function(observed, lo, hi){
  if (length(observed) != length(lo) || length(observed) != length(hi)){
    stop("observed, lo, and hi must be the same length")
  }
  mean(observed >= lo & observed <= hi, na.rm = TRUE)
}

# ---- Weighted interval score --------------------------------------------

#' Interval score for a single central prediction interval
#'
#' The building block of the weighted interval score
#' (Gneiting & Raftery 2007).
#'
#' @param observed realized value
#' @param lo,hi interval bounds
#' @param level the interval's alpha (e.g. 0.025 for a 95\% interval)
#' @return a single numeric value (lower is better)
#' @examples
#' interval_score(observed = 5, lo = 2, hi = 8, level = 0.05)
#' @export
interval_score <- function(observed, lo, hi, level){
  width <- hi - lo
  penalty_lo <- (2/level) * pmax(lo - observed, 0)
  penalty_hi <- (2/level) * pmax(observed - hi, 0)
  width + penalty_lo + penalty_hi
}

#' Weighted interval score (WIS)
#'
#' Bracher J, Ray EL, Gneiting T, Reich NG. (2021) "Evaluating epidemic
#' forecasts in an interval format." PLOS Computational Biology 17(2):
#' e1008618. \doi{10.1371/journal.pcbi.1008618}
#'
#' @param observed numeric vector of realized values, length T
#' @param median_pred numeric vector of median (or mean, if no median is
#'   available) point predictions, length T
#' @param quantiles a named list, as returned by
#'   \code{rt_forecast()$forecast_quantiles} or directly from
#'   [recurs_predict_quantiles()]$quantiles: each element is a 2 x T
#'   matrix with rows lo/hi for a (level, 1-level) interval.
#' @return a data.frame with per-time-point WIS and its components,
#'   plus the overall mean WIS as an attribute (\code{"mean_wis"})
#' @examples
#' q <- list("0.025" = matrix(c(2, 8, 15, 25), nrow = 2),
#'           "0.25"  = matrix(c(4, 6, 18, 22), nrow = 2))
#' result <- wis(observed = c(5, 20), median_pred = c(5, 21), quantiles = q)
#' attr(result, "mean_wis")
#' @export
wis <- function(observed, median_pred, quantiles){
  levels <- as.numeric(names(quantiles))
  K <- length(levels)
  Tn <- length(observed)
  for (lev_name in names(quantiles)){
    if (ncol(quantiles[[lev_name]]) != Tn) stop("each quantile matrix must have length(observed) columns")
  }

  score_matrix <- matrix(0, Tn, K)
  for (k in seq_len(K)){
    lev <- levels[k]
    lo <- quantiles[[k]][1, ]; hi <- quantiles[[k]][2, ]
    score_matrix[, k] <- interval_score(observed, lo, hi, lev)
  }

  # WIS = (1/(K+0.5)) * [0.5*|observed - median| + sum_k (level_k/2)*IS_k]
  abs_err <- abs(observed - median_pred)
  weighted_interval_sum <- rowSums(sweep(score_matrix, 2, levels/2, `*`))
  wis_per_t <- (0.5*abs_err + weighted_interval_sum) / (K + 0.5)

  out <- data.frame(observed = observed, median_pred = median_pred,
                     abs_error = abs_err, wis = wis_per_t)
  attr(out, "mean_wis") <- mean(wis_per_t, na.rm = TRUE)
  out
}

# ---- Batch scoring across asof_week folders ------------------------------

#' Score realized-vs-forecast performance across weekly output batches
#'
#' Recursively finds every \verb{asof_week*/epifilter_forecast_next_week.csv}
#' under \code{root_dir} and matches each batch's one-week-ahead forecast to
#' the actual case count once it appears as an observed row in a *later*
#' batch's \code{epifilter_results.csv} - i.e. this scores genuine
#' prospective (not in-sample) forecast performance over time.
#'
#' Surveillance case counts are commonly revised between runs as reporting
#' catches up (backfill), so the same (location, time) can legitimately
#' have different case counts across different asof_week batches. This
#' function makes that choice explicit via \code{revision}: \code{"first"}
#' uses the count as it was first reported (matches classic prospective
#' forecast evaluation); \code{"latest"} uses the most complete/least
#' revised count available.
#'
#' This function is a filesystem helper for scoring output produced by the
#' standalone \code{run_rtforecast.R} script's weekly batch folders. If
#' you're using [rt_forecast()] directly in R (recommended for new code),
#' accumulate a data.frame of your own forecasts/actuals across weeks and
#' use [mae()], [coverage()], and [wis()] directly instead.
#'
#' @param root_dir path under which to search recursively for asof_week*
#'   folders (e.g. \code{"output"})
#' @param revision \code{"first"} or \code{"latest"} - see above
#' @return a data.frame, one row per batch whose forecast has since been
#'   realized, with location, time, forecast, actual, and per-point
#'   MAE/coverage (95\% and 50\%) columns.
#' @export
score_batches <- function(root_dir, revision = c("first", "latest")){
  revision <- match.arg(revision)
  batch_dirs <- list.dirs(root_dir, recursive = TRUE)
  batch_dirs <- batch_dirs[grepl("asof_week[0-9]+$", basename(batch_dirs))]
  if (length(batch_dirs) == 0) stop("no asof_week* batches found under ", root_dir)

  location_of <- function(dataset_folder_name){
    sub("_[0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]{4}-[0-9]{2}-[0-9]{2}-trimmed$", "", dataset_folder_name)
  }
  asof_of <- function(batch_folder_name){
    as.numeric(sub("^asof_week", "", batch_folder_name))
  }

  forecasts <- do.call(rbind, lapply(batch_dirs, function(d){
    f <- file.path(d, "epifilter_forecast_next_week.csv")
    if (!file.exists(f)) return(NULL)
    x <- utils::read.csv(f)
    x$batch <- basename(d)
    x$location <- location_of(basename(dirname(d)))
    x
  }))
  observed_all <- do.call(rbind, lapply(batch_dirs, function(d){
    f <- file.path(d, "epifilter_results.csv")
    if (!file.exists(f)) return(NULL)
    x <- utils::read.csv(f)
    x$location <- location_of(basename(dirname(d)))
    x$asof <- asof_of(basename(d))
    x[, c("location", "time", "cases", "asof")]
  }))

  observed_all <- observed_all[order(observed_all$location, observed_all$time, observed_all$asof), ]
  if (revision == "first"){
    observed_all <- observed_all[!duplicated(observed_all[, c("location", "time")]), ]
  } else {
    observed_all <- observed_all[!duplicated(observed_all[, c("location", "time")], fromLast = TRUE), ]
  }
  observed_all$asof <- NULL

  merged <- merge(forecasts, observed_all, by = c("location", "time"), all.x = TRUE)
  merged <- merged[!is.na(merged$cases), ]
  if (nrow(merged) == 0){
    warning("no forecasts have a realized actual value yet (all forecast weeks are still in the future)")
    return(merged)
  }

  merged$abs_error <- abs(merged$cases - merged$pred_mean)
  merged$in_95 <- merged$cases >= merged$lo95 & merged$cases <= merged$hi95
  merged$in_50 <- merged$cases >= merged$lo50 & merged$cases <= merged$hi50

  merged[order(merged$location, merged$time),
         c("location", "batch", "time", "pred_mean", "cases",
           "abs_error", "lo95", "hi95", "in_95", "lo50", "hi50", "in_50")]
}
