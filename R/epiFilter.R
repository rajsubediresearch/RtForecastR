#' Bayesian recursive filtering via EpiFilter (internal)
#'
#' Unmodified port of the EpiFilter causal (real-time) filtering
#' algorithm.
#'
#' From: Parag KV. (2021) "Improved estimation of time-varying
#' reproduction numbers at low case incidence and between epidemic
#' waves." PLOS Computational Biology 17(9): e1009347.
#' \doi{10.1371/journal.pcbi.1009347}
#' \url{https://github.com/kpzoo/EpiFilter}
#'
#' Assumptions: observation model is the Poisson renewal equation (as
#' in EpiEstim); the reproduction number state-space model is a simple
#' diffusion.
#'
#' @param Rgrid grid on reproduction numbers
#' @param m size of Rgrid
#' @param eta diffusion (state) noise, controls smoothness of R_t
#' @param pR0 prior distribution over Rgrid at t = 1
#' @param nday number of time points
#' @param Lday total infectiousness at each time point (see [compute_lambda()])
#' @param Iday incidence at each time point
#' @param a tail probability defining the reported credible interval
#'   (e.g. 0.025 for a 95\% interval)
#' @return list(Rmed, Rhat, Rmean, pR, pRup, pstate): median, 50\%/95\%
#'   quantiles, mean, causal posterior over R, pre-update posterior,
#'   and the precomputed state-transition matrix
#' @keywords internal
epi_filter <- function(Rgrid, m, eta, pR0, nday, Lday, Iday, a){

  # Probability vector for R and prior
  pR = matrix(0, nday, m); pRup = pR
  pR[1, ] = pR0; pRup[1, ] = pR0

  # Mean and median estimates
  Rmean = rep(0, nday); Rmed = Rmean
  # 50% and 95% (depends on a) confidence on R
  Rhat = matrix(0, 4, nday)

  # Initialise mean
  Rmean[1] = pR[1, ]%*%Rgrid
  # CDF of prior
  Rcdf0 = cumsum(pR0)
  # Initialise quartiles
  idm = which(Rcdf0 >= 0.5, 1); Rmed[1] = Rgrid[idm[1]]
  id1 = which(Rcdf0 >= a, 1); id2 = which(Rcdf0 >= 1-a, 1)
  id3 = which(Rcdf0 >= 0.25, 1); id4 = which(Rcdf0 >= 0.75, 1)
  Rhat[1, 1] = Rgrid[id1[1]]; Rhat[2, 1] = Rgrid[id2[1]]
  Rhat[3, 1] = Rgrid[id3[1]]; Rhat[4, 1] = Rgrid[id4[1]]

  # Precompute state distributions for R transitions
  pstate = matrix(0, m, m);
  for(j in 1:m){
    pstate[j, ] = stats::dnorm(Rgrid[j], Rgrid, sqrt(Rgrid)*eta)
  }

  # Update prior to posterior sequentially
  for(i in 2:nday){
    # Compute mean from Poisson renewal (observation model)
    rate = Lday[i]*Rgrid
    # Probabilities of observations
    pI = stats::dpois(Iday[i], rate)

    # State predictions for R
    pRup[i, ]  = pR[i-1, ]%*%pstate
    # Update to posterior over R
    pR[i, ] = pRup[i, ]*pI
    pR[i, ] = pR[i, ]/sum(pR[i, ])

    # Posterior mean and CDF
    Rmean[i] = pR[i, ]%*%Rgrid
    Rcdf = cumsum(pR[i, ])

    # Quantiles for estimates
    idm = which(Rcdf >= 0.5, 1); Rmed[i] = Rgrid[idm[1]]
    id1 = which(Rcdf >= a, 1); id2 = which(Rcdf >= 1-a, 1)
    id3 = which(Rcdf >= 0.25, 1); id4 = which(Rcdf >= 0.75, 1)
    Rhat[1, i] = Rgrid[id1[1]]; Rhat[2, i] = Rgrid[id2[1]]
    Rhat[3, i] = Rgrid[id3[1]]; Rhat[4, i] = Rgrid[id4[1]]
  }

  # Main outputs: estimates of R and states
  list(Rmed = Rmed, Rhat = Rhat, Rmean = Rmean, pR = pR, pRup = pRup, pstate = pstate)
}
