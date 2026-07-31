#' Bayesian recursive one-step-ahead prediction via EpiFilter (internal)
#'
#' Modified from the original `recursPredict.R` (Parag KV, EpiFilter):
#' the original hardcodes its internal prediction grid to `0:800` and
#' errors out if any predicted mean incidence exceeds 720. This version
#' takes the grid's upper bound as a parameter (`maxI`) instead, so it
#' doesn't break on larger outbreaks. The Poisson-renewal/quantile
#' logic is otherwise unchanged from the original.
#'
#' From: Parag KV. (2021) "Improved estimation of time-varying
#' reproduction numbers at low case incidence and between epidemic
#' waves." PLOS Computational Biology 17(9): e1009347.
#' \doi{10.1371/journal.pcbi.1009347}
#' \url{https://github.com/kpzoo/EpiFilter}
#'
#' @param Rgrid grid on reproduction numbers
#' @param pR posterior over R (filtered or smoothed), from [epi_filter()]/[epi_smoother()]
#' @param Lday total infectiousness at each time point
#' @param Rmean mean R estimate at each time point
#' @param a tail probability defining the reported credible interval
#' @param maxI upper bound of the internal prediction grid; set
#'   comfortably above the data's peak incidence
#' @return list(pred, predInt): mean one-step-ahead prediction, and a
#'   4 x (n-1) matrix of quantiles (rows: lo95, hi95, lo50, hi50)
#' @keywords internal
recurs_predict <- function(Rgrid, pR, Lday, Rmean, a, maxI = 800){

  # Grid size and length of time series
  nday = nrow(pR); m = ncol(pR)
  if (length(Rgrid) != m | length(Lday) != nday){
    stop("Input vectors of incorrect dimension")
  }

  # Mean prediction: Lday[i] => Iday[i+1]
  pred = Lday*Rmean; pred = pred[1:length(pred)-1]

  # Discrete space of possible predictions
  Igrid = 0:maxI; lenI = length(Igrid);

  if (any(pred > 0.9*max(Igrid))){
    stop(paste0("Epidemic size too large for maxI=", maxI,
                " - increase maxI"))
  }

  # Prediction cdf and quantiles (50% and 95%)
  Fpred = matrix(0, nday-1, lenI)
  predInt = matrix(0, 4, nday-1)

  # At every time construct CDF of predictions
  for(i in 1:(nday-1)){
    rate = Lday[i]*Rgrid
    # Vectorized over Igrid: M[g, k] = dpois(Igrid[k], rate[g]) for each
    # R-grid point g and case-count k, then marginalise over R's
    # posterior pR[i, ] in one matrix-vector product - same computation
    # as calling dpois() once per Igrid value in a loop, just done in
    # one vectorized pass instead of lenI separate calls
    M = outer(rate, Igrid, function(r, k) stats::dpois(k, r))
    pI = as.vector(pR[i, ] %*% M)
    Fpred[i, ] = cumsum(pI)/sum(pI)
    id1 = which(Fpred[i, ] >= a); id2 = which(Fpred[i, ] >= 1-a)
    id3 = which(Fpred[i, ] >= 0.25); id4 = which(Fpred[i, ] >= 0.75)
    predInt[1, i] = Igrid[id1[1]]; predInt[2, i] = Igrid[id2[1]]
    predInt[3, i] = Igrid[id3[1]]; predInt[4, i] = Igrid[id4[1]]
  }
  list(pred = pred, predInt = predInt)
}
