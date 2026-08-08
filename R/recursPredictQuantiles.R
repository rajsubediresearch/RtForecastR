#' Bayesian recursive prediction with arbitrary quantile levels (internal)
#'
#' Generalizes [recurs_predict()] to report an arbitrary set of
#' quantile levels rather than just 50%/95%, needed for the weighted
#' interval score ([wis()]). Same Poisson-renewal/posterior
#' marginalisation logic; this only widens what's read off the
#' already-computed CDF, so it never changes point predictions or the
#' 50%/95% values a call to [recurs_predict()] would give, only how
#' many quantiles are reported alongside them.
#'
#' Bracher J, Ray EL, Gneiting T, Reich NG. (2021) "Evaluating epidemic
#' forecasts in an interval format." PLOS Computational Biology 17(2):
#' e1008618. \doi{10.1371/journal.pcbi.1008618}
#'
#' @param Rgrid grid on reproduction numbers
#' @param pR posterior over R (filtered or smoothed)
#' @param Lday total infectiousness at each time point
#' @param Rmean mean R estimate at each time point
#' @param quantile_levels vector of probabilities in (0, 0.5), each
#'   producing a nested (level, 1-level) interval. Default matches
#'   [recurs_predict()]'s 50%/95% reporting; pass e.g.
#'   `seq(0.05, 0.45, by = 0.05)` to add the levels a WIS calculation
#'   typically wants.
#' @param maxI upper bound of the internal prediction grid
#' @return list(pred, quantiles): mean one-step-ahead prediction, and a
#'   named list where `quantiles[[as.character(level)]]` is a 2 x (n-1)
#'   matrix with rows lo/hi for that (level, 1-level) interval
#' @keywords internal
recurs_predict_quantiles <- function(Rgrid, pR, Lday, Rmean,
                                      quantile_levels = c(0.025, 0.25), maxI = 800){

  nday = nrow(pR); m = ncol(pR)
  if (length(Rgrid) != m | length(Lday) != nday){
    stop("Input vectors of incorrect dimension")
  }
  if (any(quantile_levels <= 0 | quantile_levels >= 0.5)){
    stop("quantile_levels must each be strictly between 0 and 0.5 (each defines a (level, 1-level) interval)")
  }

  pred = Lday*Rmean; pred = pred[1:length(pred)-1]

  Igrid = 0:maxI; lenI = length(Igrid);
  if (any(pred > 0.9*max(Igrid))){
    stop(paste0("Epidemic size too large for maxI=", maxI,
                " - increase maxI"))
  }

  quantile_levels = sort(quantile_levels)
  nlev = length(quantile_levels)
  quantiles = vector("list", nlev)
  names(quantiles) = as.character(quantile_levels)
  for (k in 1:nlev) quantiles[[k]] = matrix(0, 2, nday-1)

  Fpred = matrix(0, nday-1, lenI)

  for(i in 1:(nday-1)){
    rate = Lday[i]*Rgrid
    M = outer(rate, Igrid, function(r, k) stats::dpois(k, r))
    pI = as.vector(pR[i, ] %*% M)
    Fpred[i, ] = cumsum(pI)/sum(pI)

    for (k in 1:nlev){
      lev = quantile_levels[k]
      id_lo = which(Fpred[i, ] >= lev)
      id_hi = which(Fpred[i, ] >= 1 - lev)
      quantiles[[k]][1, i] = Igrid[id_lo[1]]
      quantiles[[k]][2, i] = Igrid[id_hi[1]]
    }
  }

  list(pred = pred, quantiles = quantiles)
}
