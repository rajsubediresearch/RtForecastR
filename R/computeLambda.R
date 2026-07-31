#' Total infectiousness (renewal-equation convolution)
#'
#' Computes the "total infectiousness" Lambda_t: the convolution of
#' past incidence with a discretized Gamma(mean, variance) generation
#' interval distribution, as required by [epiFilter()] and
#' [epiSmoother()] alongside raw incidence. Written to avoid a
#' dependency on the EpiEstim package's `overall_infectivity()` for
#' this single step.
#'
#' The Gamma(mean, variance) discretization approach follows the
#' convention used in `get_Rt.m`, part of the GrowthPredict toolbox:
#' Chowell G, Bleichrodt A, Dahal S, Tariq A, Roosa K, Hyman JM, Luo R.
#' (2024) "GrowthPredict: A toolbox and tutorial-based primer for
#' fitting and forecasting growth trajectories using phenomenological
#' growth models." Scientific Reports 14, 1630.
#' \doi{10.1038/s41598-024-51852-8}
#'
#' @param incidence numeric vector of case counts, length n
#' @param mean_GI mean generation interval, in the same time units as
#'   `incidence` (e.g. weeks if `incidence` is weekly)
#' @param var_GI variance of the generation interval, same time units
#' @param max_si maximum serial interval to consider (defaults to
#'   `length(incidence) - 1`)
#' @return numeric vector of length n; `Lambda[1]` is `NA` (no prior
#'   incidence exists to convolve with)
#' @examples
#' data(measles_cdmx)
#' Lambda <- compute_lambda(measles_cdmx$cases, mean_GI = 11/7, var_GI = (4/7)^2)
#' @export
compute_lambda <- function(incidence, mean_GI, var_GI, max_si = NULL) {
  n <- length(incidence)
  if (is.null(max_si)) max_si <- n - 1

  # Gamma(mean, var) discretization, method-of-moments parameterization:
  # shape a = mean^2/var, scale b = var/mean
  b <- var_GI / mean_GI
  a <- mean_GI / b

  # w_s = F(s) - F(s-1), s = 1, ..., max_si  (F = Gamma CDF)
  s <- 0:max_si
  Fs <- stats::pgamma(s, shape = a, scale = b)
  w <- diff(Fs)          # length max_si, w[1] corresponds to s=1, etc.
  w <- w / sum(w)        # normalize to a proper pmf

  Lambda <- rep(NA_real_, n)
  for (t in 2:n) {
    smax <- min(t - 1, max_si)
    Lambda[t] <- sum(incidence[(t - 1):(t - smax)] * w[1:smax])
  }
  Lambda
}
