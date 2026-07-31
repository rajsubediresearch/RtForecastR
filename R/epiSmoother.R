#' Bayesian recursive smoothing via EpiFilter (internal)
#'
#' Unmodified port of the EpiFilter retrospective (smoothed) estimation
#' algorithm. Must be run after [epi_filter()], using its outputs.
#'
#' From: Parag KV. (2021) "Improved estimation of time-varying
#' reproduction numbers at low case incidence and between epidemic
#' waves." PLOS Computational Biology 17(9): e1009347.
#' \doi{10.1371/journal.pcbi.1009347}
#' \url{https://github.com/kpzoo/EpiFilter}
#'
#' @param Rgrid grid on reproduction numbers
#' @param m size of Rgrid
#' @param pR filtered (causal) posterior over R, from [epi_filter()]
#' @param pRup pre-update posterior, from [epi_filter()]
#' @param nday number of time points
#' @param pstate state-transition matrix, from [epi_filter()]
#' @param a tail probability defining the reported credible interval
#' @return list(Rmed, Rhat, Rmean, qR): median, 50\%/95\% quantiles,
#'   mean, and the smoothed (backward+forward) posterior over R
#' @keywords internal
epi_smoother <- function(Rgrid, m, pR, pRup, nday, pstate, a){

  # Last smoothed distribution same as filtered
  qR = matrix(0, nday, m); qR[nday, ] = pR[nday, ]

  # Main smoothing equation iteratively computed
  for(i in seq(nday-1, 1)){
    # Remove zeros
    pRup[i+1, pRup[i+1, ] == 0] = 10^-8

    # Integral term in smoother
    integ = qR[i+1, ]/pRup[i+1, ]
    integ = integ%*%pstate

    # Smoothed posterior over Rgrid
    qR[i, ] = pR[i, ]*integ
    # Force a normalisation
    qR[i, ] = qR[i, ]/sum(qR[i, ]);
  }

  # Mean, median estimates of R
  Rmean = rep(0, nday); Rmed = Rmean
  # 50% and 95% (depends on a) confidence on R
  Rhat = matrix(0, 4, nday)

  # Compute at every time point
  for (i in 1:nday) {
    # Posterior mean and CDF
    Rmean[i] = qR[i, ]%*%Rgrid
    Rcdf = cumsum(qR[i, ])

    # Quantiles for estimates
    idm = which(Rcdf >= 0.5); Rmed[i] = Rgrid[idm[1]]
    id1 = which(Rcdf >= a, 1); id2 = which(Rcdf >= 1-a, 1)
    id3 = which(Rcdf >= 0.25, 1); id4 = which(Rcdf >= 0.75, 1)
    Rhat[1, i] = Rgrid[id1[1]]; Rhat[2, i] = Rgrid[id2[1]]
    Rhat[3, i] = Rgrid[id3[1]]; Rhat[4, i] = Rgrid[id4[1]]
  }

  # Main outputs: estimates of R and states
  list(Rmed = Rmed, Rhat = Rhat, Rmean = Rmean, qR = qR)
}
