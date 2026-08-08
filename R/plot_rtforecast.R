#' Plot method for rtforecast objects
#'
#' Draws to the current graphics device (matches normal R plotting
#' conventions - use [grDevices::png()]/[grDevices::pdf()] yourself
#' around the call if you want a file, the way you would for any base R
#' plot). Ports the "current situation" view from the original
#' `run_rtforecast.R` script.
#'
#' @param x an `"rtforecast"` object, as returned by [rt_forecast()]
#' @param which one of `"Rt"` (filtered vs smoothed R_t), `"forecast"`
#'   (recent observed/predicted cases + one-step forecast), or
#'   `"observed_vs_predicted"` (full-period observed vs in-sample
#'   predicted + forecast)
#' @param n_recent for `which = "forecast"`, how many recent time
#'   points to show (default 8)
#' @param ... passed on to the underlying [graphics::plot()] call
#' @return The input `"rtforecast"` object `x`, returned invisibly.
#'   Called primarily for its side effect of drawing to the current
#'   graphics device.
#' @export
plot.rtforecast <- function(x, which = c("Rt", "forecast", "observed_vs_predicted"),
                             n_recent = 8, ...){
  which <- match.arg(which)
  raw <- x$raw
  time_r <- raw$time[raw$tvec]

  if (which == "Rt"){
    Rfilt <- raw$Rfilt; Rsmooth <- raw$Rsmooth
    graphics::plot(time_r, Rfilt$Rmean, type = 'l', col = 'blue', lwd = 2,
         ylim = c(0, max(Rsmooth$Rhat[2,], Rfilt$Rhat[2,], na.rm = TRUE)),
         xlab = 'Time', ylab = 'R_t', main = 'Filtered (real-time) vs Smoothed R_t', ...)
    graphics::lines(time_r, Rfilt$Rhat[1,], col = 'blue', lty = 2)
    graphics::lines(time_r, Rfilt$Rhat[2,], col = 'blue', lty = 2)
    graphics::lines(time_r, Rsmooth$Rmean, col = 'red', lwd = 2)
    graphics::lines(time_r, Rsmooth$Rhat[1,], col = 'red', lty = 2)
    graphics::lines(time_r, Rsmooth$Rhat[2,], col = 'red', lty = 2)
    graphics::abline(h = 1, col = 'black', lty = 3)
    graphics::legend('topright', c('Filtered','Smoothed','R=1'),
                      col = c('blue','red','black'), lty = c(1,1,3))
    return(invisible(x))
  }

  time_pred <- raw$time[raw$tvec_pred]
  cases_pred <- raw$cases[raw$tvec_pred]
  fc <- x$forecast

  if (which == "forecast"){
    idx <- utils::tail(seq_along(time_pred), n_recent)
    graphics::plot(time_pred[idx], cases_pred[idx], type = 'n',
         xlim = c(min(time_pred[idx]), fc$time + 0.5),
         ylim = c(0, max(c(cases_pred[idx], x$predictions$pred_hi95[idx], fc$hi95), na.rm = TRUE)*1.1),
         xlab = 'Time', ylab = 'Cases',
         main = paste0('Current situation (last ', n_recent, ' points) + forecast'), ...)
    xx <- time_pred[idx]
    graphics::polygon(c(xx, rev(xx)),
             c(x$predictions$pred_lo95[idx], rev(x$predictions$pred_hi95[idx])),
             col = grDevices::rgb(0, 0, 1, 0.15), border = NA)
    graphics::lines(xx, x$predictions$pred_next[idx], col = 'blue', lwd = 2)
    graphics::points(xx, cases_pred[idx], pch = 19, col = 'black', cex = 1.2)
    graphics::abline(v = raw$time[length(raw$time)], col = 'black', lty = 3)
    graphics::points(fc$time, fc$pred_mean, col = 'red', pch = 19, cex = 1.5)
    graphics::arrows(fc$time, fc$lo95, fc$time, fc$hi95, col = 'red', angle = 90, code = 3, length = 0.05, lwd = 2)
    graphics::arrows(fc$time, fc$lo50, fc$time, fc$hi50, col = 'red', angle = 90, code = 3, length = 0, lwd = 4)
    graphics::legend('topleft', c('Observed','Predicted (in-sample)','95% CI','Forecast'),
                      col = c('black','blue', grDevices::rgb(0,0,1,0.4), 'red'),
                      pch = c(19, NA, 15, 19), lty = c(NA, 1, NA, NA), cex = 0.75)
    return(invisible(x))
  }

  if (which == "observed_vs_predicted"){
    all_times <- raw$time; all_cases <- raw$cases; n <- length(all_cases)
    pred_aligned <- rep(NA, n); pred_lo <- rep(NA, n); pred_hi <- rep(NA, n)
    match_idx <- match(raw$tvec_pred, seq_len(n))
    pred_aligned[match_idx] <- x$predictions$pred_next
    pred_lo[match_idx] <- x$predictions$pred_lo95
    pred_hi[match_idx] <- x$predictions$pred_hi95

    graphics::plot(all_times, all_cases, type = 'n',
         xlim = c(min(all_times), fc$time + 0.5),
         ylim = c(0, max(c(all_cases, pred_hi, fc$hi95), na.rm = TRUE)*1.05),
         xlab = 'Time', ylab = 'Cases', main = 'Observed vs predicted + forecast', ...)
    valid <- which(!is.na(pred_lo))
    graphics::polygon(c(all_times[valid], rev(all_times[valid])),
             c(pred_lo[valid], rev(pred_hi[valid])),
             col = grDevices::rgb(0, 0, 1, 0.15), border = NA)
    graphics::lines(all_times, pred_aligned, col = 'blue', lwd = 2)
    graphics::points(all_times, all_cases, pch = 19, col = 'black')
    graphics::points(all_times[1], all_cases[1], col = 'black', pch = 1, cex = 1.5)
    graphics::abline(v = raw$time[length(raw$time)], col = 'black', lty = 3)
    graphics::points(fc$time, fc$pred_mean, col = 'red', pch = 19, cex = 1.3)
    graphics::arrows(fc$time, fc$lo95, fc$time, fc$hi95, col = 'red', angle = 90, code = 3, length = 0.05, lwd = 2)
    graphics::arrows(fc$time, fc$lo50, fc$time, fc$hi50, col = 'red', angle = 90, code = 3, length = 0, lwd = 4)
    graphics::legend('topright', c('Observed','Predicted (in-sample)','95% CI','Forecast'),
                      col = c('black','blue', grDevices::rgb(0,0,1,0.4), 'red'),
                      pch = c(19, NA, 15, 19), lty = c(NA, 1, NA, NA), cex = 0.75)
    return(invisible(x))
  }
}
