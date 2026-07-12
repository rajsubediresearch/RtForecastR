######################################################################
## RtForecastR - one-step-ahead R_t estimation and forecasting
## Author: Raj Subedi
##
## Core algorithm (core/epiFilter.R, epiSmoother.R, recursPredict.R) is
## from EpiFilter, unmodified:
##   Parag KV. (2021) "Improved estimation of time-varying reproduction
##   numbers at low case incidence and between epidemic waves."
##   PLOS Computational Biology 17(9): e1009347.
##   https://github.com/kpzoo/EpiFilter
##
## Generation-interval discretization (core/computeLambda.R) follows the
## standard convention used in GrowthPredict's get_Rt.m:
##   https://github.com/gchowell/forecasting_growthmodels
##   Chowell G, et al. (2024) Scientific Reports 14, 1630.
##
## Everything else here - the out-of-sample forecast extension, the
## three-plot reporting layout, and the batch-organized output structure
## - is original work built on top of the above.
##
## Disease-agnostic: nothing in the core logic is specific to any one
## pathogen. Set DATA_FILE, mean_GI, and var_GI below for your own data.
## The values shipped here (examples/measles_cdmx_example.txt, and a
## measles-typical generation interval) are a worked example, not a
## fixed assumption of the tool itself.
######################################################################

suppressMessages(library("caTools"))

# Source all EpiFilter main functions (including computeLambda.R)
files.sources = list.files(path = "core")
for (i in 1:length(files.sources)) {
  source(paste0(c("core/", files.sources[i]), collapse = ''))
}

# ---- 1. Load data ---------------------------------------------------
# Two columns, no header: time index, case count
DATA_FILE <- "examples/JALISCO_2025-08-18_2026-06-15-trimmed.txt"
data1 <- read.table(DATA_FILE, header = FALSE)
timevect <- data1[[1]]
I <- data1[[2]]
n <- length(I)

# ---- 1b. Set up organized output folder ------------------------------
# ./output/<dataset_name>/asof_week<last_observed_week>/
# Keying the batch folder off the last observed week (rather than, say,
# today's date) means re-running with a new week of data automatically
# lands in its own folder - so nothing gets overwritten, and each
# batch's outputs are tied to exactly what data went into it, which
# matters for actually tracking forecast performance over time.

dataset_name <- tools::file_path_sans_ext(basename(DATA_FILE))
batch_name <- paste0("asof_week", tail(timevect, 1))
outdir <- file.path("output", dataset_name, batch_name)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Small helper so every output-writing call below just needs a filename
outpath <- function(filename) file.path(outdir, filename)

cat("Output folder for this run:", outdir, "\n\n")

# ---- 2. Generation interval (set for your disease/dataset) -----------
# EXAMPLE VALUES BELOW ARE FOR MEASLES - replace mean_GI/var_GI with the
# appropriate literature value for whatever disease your DATA_FILE
# contains. Nothing else in this script is disease-specific.
#
# NOTE: measles' generation interval is typically cited around ~11-12
# days (varies by source), NOT the 5/7-week value in the existing
# options_Rt.m - that file's comment suggests DT-scaled units, worth
# double-checking against your actual DT convention before trusting
# these R_t estimates for anything beyond exploration. Adjust below to
# match your data's actual time units (weeks here, given weekly data).
mean_GI <- 11/7   # EXAMPLE (measles): GI ~11 days, converted to weeks
var_GI  <- (4/7)^2  # placeholder SD ~4 days, converted to weeks - adjust
                     # to a literature value if you have a better estimate

Lam <- computeLambda(I, mean_GI, var_GI)

# ---- 3. EpiFilter setup ----------------------------------------------
Rmin <- 0.01; Rmax <- 10; eta <- 0.1   # eta = state noise (smoothness)
m <- 200
pR0 <- (1/m) * rep(1, m)
Rgrid <- seq(Rmin, Rmax, length.out = m)
a <- 0.025  # -> 95% credible intervals

tvec <- 2:n  # Lambda[1] is undefined, start from t=2

# ---- 4. Run filtering (real-time/causal) ------------------------------
Rfilt <- epiFilter(Rgrid, m, eta, pR0, length(tvec), Lam[tvec], I[tvec], a)
# Rfilt: [Rmed, Rhatci, Rmean, pR, pRup, pstate]

# ---- 5. One-step-ahead predictions (model adequacy check) -------------
Ifilt <- recursPredict(Rgrid, Rfilt[[4]], Lam[tvec], Rfilt[[3]], a)
# Ifilt: [pred, predci]

# ---- 6. Run smoothing (retrospective, uses all data) -------------------
Rsmooth <- epiSmoother(Rgrid, m, Rfilt[[4]], Rfilt[[5]], length(tvec), Rfilt[[6]], a)
# Rsmooth: [Rmed, Rhatci, Rmean, qR]

Ismooth <- recursPredict(Rgrid, Rsmooth[[4]], Lam[tvec], Rsmooth[[3]], a)

# ---- 6b. Genuine out-of-sample forecast: one week beyond the data -------
# Unlike the in-sample one-step-ahead predictions above (every one of
# which has observed data available, since they're computed using R
# estimates from *within* the fitted series), this extends prediction ONE
# step past the last observed week - using only data through the present,
# same underlying Poisson-renewal logic as recursPredict.R, applied by
# hand for this single future point.

# Total infectiousness for the NEXT (unobserved) week, using only
# observed incidence through the present
I_extended <- c(I, NA)                       # placeholder for the next week
Lam_next <- computeLambda(I_extended, mean_GI, var_GI)
Lam_forecast <- tail(Lam_next, 1)             # Lambda for week n+1

# Posterior over R at the most recent (last observed) time point -
# our best current real-time estimate of transmission
pR_last <- Rfilt[[4]][nrow(Rfilt[[4]]), ]

# Predictive distribution for next week's case count, marginalising the
# Poisson renewal rate over the full current posterior on R (identical
# logic to the inner loop of recursPredict.R, applied for one extra step)
Igrid <- 0:800
rate <- Lam_forecast * Rgrid
pI <- sapply(Igrid, function(k) sum(dpois(k, rate) * pR_last))
Fpred <- cumsum(pI) / sum(pI)

forecast_mean  <- Lam_forecast * Rfilt[[3]][length(Rfilt[[3]])]
forecast_lo95  <- Igrid[which(Fpred >= a)[1]]
forecast_hi95  <- Igrid[which(Fpred >= 1 - a)[1]]
forecast_lo50  <- Igrid[which(Fpred >= 0.25)[1]]
forecast_hi50  <- Igrid[which(Fpred >= 0.75)[1]]
forecast_time  <- timevect[n] + 1             # one week beyond last observed

cat("\nOut-of-sample forecast for week", forecast_time, "(no observed data yet):\n")
cat(sprintf("  Predicted cases: %.1f (95%% CI: %d-%d, 50%% CI: %d-%d)\n",
            forecast_mean, forecast_lo95, forecast_hi95, forecast_lo50, forecast_hi50))
cat("  Based on current R_t =", round(Rfilt[[3]][length(Rfilt[[3]])], 3),
    "(", round(Rfilt[[2]][1, ncol(Rfilt[[2]])], 3), "-",
    round(Rfilt[[2]][4, ncol(Rfilt[[2]])], 3), ")\n")

# ---- 7. Elimination probability: P(R_t < 1) at each time point --------
# Computed directly from the posterior grid pR (filtered) / qR (smoothed)
below1_filt   <- Rgrid < 1
below1_smooth <- Rgrid < 1
Pelim_filt   <- apply(Rfilt[[4]],   1, function(p) sum(p[below1_filt]))
Pelim_smooth <- apply(Rsmooth[[4]], 1, function(p) sum(p[below1_smooth]))

# ---- 8. Summary table --------------------------------------------------
# NOTE: recursPredict returns one fewer point than tvec (the first
# filtered time point has no prior R estimate to generate a one-step-
# ahead prediction from), so predictions align to tvec[-1], not tvec.
tvec_pred <- tvec[-1]

results <- data.frame(
  time          = timevect[tvec],
  cases         = I[tvec],
  R_filtered    = round(Rfilt[[3]], 3),
  R_filt_lo95   = round(Rfilt[[2]][1, ], 3),
  R_filt_hi95   = round(Rfilt[[2]][4, ], 3),
  R_smoothed    = round(Rsmooth[[3]], 3),
  R_smooth_lo95 = round(Rsmooth[[2]][1, ], 3),
  R_smooth_hi95 = round(Rsmooth[[2]][4, ], 3),
  P_elim_filt   = round(Pelim_filt, 3),
  P_elim_smooth = round(Pelim_smooth, 3)
)

pred_table <- data.frame(
  time          = timevect[tvec_pred],
  cases         = I[tvec_pred],
  pred_next     = round(Ifilt[[1]], 1),
  pred_lo95     = round(Ifilt[[2]][1, ], 1),
  pred_hi95     = round(Ifilt[[2]][4, ], 1)
)

write.csv(pred_table, outpath("epifilter_predictions.csv"), row.names = FALSE)

# Small standalone record of the forward forecast - designed to be
# compared against next week's actual case count once it's available,
# to track real forecast performance over time
forecast_record <- data.frame(
  time      = forecast_time,
  pred_mean = round(forecast_mean, 1),
  lo95      = forecast_lo95, hi95 = forecast_hi95,
  lo50      = forecast_lo50, hi50 = forecast_hi50,
  R_used    = round(Rfilt[[3]][length(Rfilt[[3]])], 3),
  date_generated = as.character(Sys.Date())
)
write.csv(forecast_record, outpath("epifilter_forecast_next_week.csv"), row.names = FALSE)

write.csv(results, outpath("epifilter_results.csv"), row.names = FALSE)
cat("R_t estimates:\n"); print(results)
cat("\nOne-step-ahead predictions:\n"); print(pred_table)

# ---- 9. Plot 1: Default R_t plot (filtered vs smoothed) -----------------
png(outpath("epifilter_1_default_Rt.png"), width = 900, height = 500)
plot(timevect[tvec], Rfilt[[3]], type = 'l', col = 'blue', lwd = 2,
     ylim = c(0, max(Rsmooth[[2]][4,], Rfilt[[2]][4,], na.rm=TRUE)),
     xlab = 'Time', ylab = 'R_t', main = 'Filtered (real-time) vs Smoothed R_t')
lines(timevect[tvec], Rfilt[[2]][1,], col = 'blue', lty = 2)
lines(timevect[tvec], Rfilt[[2]][4,], col = 'blue', lty = 2)
lines(timevect[tvec], Rsmooth[[3]], col = 'red', lwd = 2)
lines(timevect[tvec], Rsmooth[[2]][1,], col = 'red', lty = 2)
lines(timevect[tvec], Rsmooth[[2]][4,], col = 'red', lty = 2)
abline(h = 1, col = 'black', lty = 3)
legend('topright', c('Filtered','Smoothed','R=1'), col=c('blue','red','black'), lty=c(1,1,3))
dev.off()

# ---- 10. Plot 2: "Current situation" snapshot ----------------------------
# A reusable plotting function so the zoomed and full-period versions stay
# perfectly consistent with each other - same logic, different time window.

plot_situation <- function(idx, title_suffix, filename, show_refline) {
  png(filename, width = 900, height = 500)
  plot(timevect[tvec_pred][idx], I[tvec_pred][idx], type = 'n',
       xlim = c(min(timevect[tvec_pred][idx]), forecast_time + 0.5),
       ylim = c(0, max(c(I[tvec_pred][idx], Ifilt[[2]][2,][idx], forecast_hi95), na.rm=TRUE)*1.1),
       xlab = 'Time', ylab = 'Cases',
       main = paste0('Current situation (', title_suffix, ') + next-week forecast'))
  # 95% prediction interval as a shaded ribbon, not two separate dashed
  # lines - clearer, and avoids any visual confusion with the dotted
  # "last observed week" boundary marker below
  xx <- timevect[tvec_pred][idx]
  polygon(c(xx, rev(xx)), c(Ifilt[[2]][1,][idx], rev(Ifilt[[2]][2,][idx])),
          col = rgb(0, 0, 1, 0.15), border = NA)
  lines(xx, Ifilt[[1]][idx], col = 'blue', lwd = 2)
  points(xx, I[tvec_pred][idx], pch = 19, col = 'black', cex = 1.2)
  abline(v = timevect[n], col = 'black', lty = 3)
  # Optional reference line: last observed value, extended only into the
  # forecast region - a passive visual anchor (not a model prediction) to
  # judge at a glance whether the forecast sits above or below "no change"
  if (show_refline) {
    segments(timevect[n], I[n], forecast_time + 0.5, I[n],
             col = 'gray50', lty = 3, lwd = 1.5)
  }
  points(forecast_time, forecast_mean, col = 'red', pch = 19, cex = 1.5)
  arrows(forecast_time, forecast_lo95, forecast_time, forecast_hi95,
         col = 'red', angle = 90, code = 3, length = 0.05, lwd = 2)
  arrows(forecast_time, forecast_lo50, forecast_time, forecast_hi50,
         col = 'red', angle = 90, code = 3, length = 0, lwd = 4)
  leg_labels <- c('Observed','Predicted (in-sample)','95% CI (in-sample)','Forecast (out-of-sample)','Last observed week')
  leg_col <- c('black','blue', rgb(0,0,1,0.4), 'red','black')
  leg_pch <- c(19,NA,15,19,NA)
  leg_lty <- c(NA,1,NA,NA,3)
  if (show_refline) {
    leg_labels <- c(leg_labels, 'Last observed value (reference)')
    leg_col <- c(leg_col, 'gray50')
    leg_pch <- c(leg_pch, NA)
    leg_lty <- c(leg_lty, 3)
  }
  legend('topleft', leg_labels, col = leg_col, pch = leg_pch, lty = leg_lty, cex = 0.75)
  dev.off()
}
# Zoomed: last n_recent weeks - both with and without the reference line
n_recent <- 8
recent_idx <- tail(seq_along(tvec_pred), n_recent)
plot_situation(recent_idx, paste0('last ', n_recent, ' weeks'),
                outpath('epifilter_2_current_situation.png'), show_refline = FALSE)
plot_situation(recent_idx, paste0('last ', n_recent, ' weeks'),
                outpath('epifilter_2_current_situation_refline.png'), show_refline = TRUE)

# ---- 11. Plot 3: Observed vs predicted, full period ----------------------
# Each week's ACTUAL case count plotted against the prediction made FOR
# that week (using the prior week's R estimate - the "one-step-ahead"
# detail belongs in the write-up's methods section, not the plot itself).
# The very first observed week has no prediction (a prediction needs a
# prior R estimate, which doesn't exist yet at the start of the series) -
# marked with a hollow circle rather than a text label, so the write-up
# can explain it in the caption.
all_times <- timevect
pred_aligned <- rep(NA, n)
pred_lo_aligned <- rep(NA, n); pred_hi_aligned <- rep(NA, n)
match_idx <- match(tvec_pred, seq_len(n))
pred_aligned[match_idx] <- Ifilt[[1]]
pred_lo_aligned[match_idx] <- Ifilt[[2]][1,]
pred_hi_aligned[match_idx] <- Ifilt[[2]][2,]

plot_obs_vs_pred <- function(filename, show_refline) {
  png(filename, width = 1000, height = 500)
  plot(all_times, I, type = 'n',
       xlim = c(min(all_times), forecast_time + 0.5),
       ylim = c(0, max(c(I, pred_hi_aligned, forecast_hi95), na.rm=TRUE)*1.05),
       xlab = 'Time', ylab = 'Cases',
       main = 'Observed vs predicted + forecast')
  # 95% prediction interval as a shaded ribbon (only over the range with a
  # defined prediction - the first week has none, so it's excluded here to
  # avoid NA breaking the polygon)
  valid <- which(!is.na(pred_lo_aligned))
  polygon(c(all_times[valid], rev(all_times[valid])),
          c(pred_lo_aligned[valid], rev(pred_hi_aligned[valid])),
          col = rgb(0, 0, 1, 0.15), border = NA)
  lines(all_times, pred_aligned, col = 'blue', lwd = 2)
  points(all_times, I, pch = 19, col = 'black')
  # First week: hollow circle marks "no prediction available" without
  # needing an on-plot text label
  points(all_times[1], I[1], col = 'black', pch = 1, cex = 1.5)
  abline(v = timevect[n], col = 'black', lty = 3)
  if (show_refline) {
    segments(timevect[n], I[n], forecast_time + 0.5, I[n],
             col = 'gray50', lty = 3, lwd = 1.5)
  }
  points(forecast_time, forecast_mean, col = 'red', pch = 19, cex = 1.3)
  arrows(forecast_time, forecast_lo95, forecast_time, forecast_hi95,
         col = 'red', angle = 90, code = 3, length = 0.05, lwd = 2)
  arrows(forecast_time, forecast_lo50, forecast_time, forecast_hi50,
         col = 'red', angle = 90, code = 3, length = 0, lwd = 4)
  leg_labels <- c('Observed','Predicted (in-sample)','95% CI (in-sample)','Forecast (out-of-sample)','Last observed week')
  leg_col <- c('black','blue', rgb(0,0,1,0.4), 'red','black')
  leg_pch <- c(19,NA,15,19,NA)
  leg_lty <- c(NA,1,NA,NA,3)
  if (show_refline) {
    leg_labels <- c(leg_labels, 'Last observed value (reference)')
    leg_col <- c(leg_col, 'gray50')
    leg_pch <- c(leg_pch, NA)
    leg_lty <- c(leg_lty, 3)
  }
  legend('topright', leg_labels, col = leg_col, pch = leg_pch, lty = leg_lty, cex = 0.75)
  dev.off()
}

plot_obs_vs_pred(outpath('epifilter_3_observed_vs_predicted.png'), show_refline = FALSE)
plot_obs_vs_pred(outpath('epifilter_3_observed_vs_predicted_refline.png'), show_refline = TRUE)

cat("\nDone. Plots + three CSVs written to:", outdir, "\n")
cat("  epifilter_1_default_Rt.png\n")
cat("  epifilter_2_current_situation.png / _refline.png (zoomed, last 8 weeks, + forecast)\n")
cat("  epifilter_3_observed_vs_predicted.png / _refline.png (full period, + forecast)\n")
cat("  epifilter_results.csv, epifilter_predictions.csv, epifilter_forecast_next_week.csv\n")
