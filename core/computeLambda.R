######################################################################
## computeLambda.R
## Author: Raj Subedi
##
## Computes the "total infectiousness" Lambda_t (renewal equation
## convolution of past incidence with the discretized generation
## interval distribution), an input EpiFilter (Parag et al.) requires
## alongside raw incidence. Written to avoid a dependency on the
## EpiEstim CRAN package's overall_infectivity() for this step.
##
## The Gamma(mean, variance) discretization approach used here follows
## the standard convention used in get_Rt.m, part of the
## GrowthPredict toolbox by Gerardo Chowell and colleagues:
##   https://github.com/gchowell/forecasting_growthmodels
##   Chowell G, Bleichrodt A, Dahal S, Tariq A, Roosa K, Hyman JM, Luo R.
##   (2024) "GrowthPredict: A toolbox and tutorial-based primer for
##   fitting and forecasting growth trajectories using phenomenological
##   growth models." Scientific Reports 14, 1630.
##   https://doi.org/10.1038/s41598-024-51852-8
######################################################################
#
# Lambda_t = sum_{s=1}^{t-1} I_{t-s} * w_s,  for t = 2, ..., n
# Lambda_1 is undefined/NA (no prior incidence to convolve with).

computeLambda <- function(incidence, mean_GI, var_GI, max_si = NULL) {
  n <- length(incidence)
  if (is.null(max_si)) max_si <- n - 1

  # Gamma(mean, var) discretization, method-of-moments parameterization:
  # shape a = mean^2/var, scale b = var/mean
  b <- var_GI / mean_GI
  a <- mean_GI / b

  # w_s = F(s) - F(s-1), s = 1, ..., max_si  (F = Gamma CDF)
  s <- 0:max_si
  Fs <- pgamma(s, shape = a, scale = b)
  w <- diff(Fs)          # length max_si, w[1] corresponds to s=1, etc.
  w <- w / sum(w)        # normalize to a proper pmf

  Lambda <- rep(NA_real_, n)
  for (t in 2:n) {
    smax <- min(t - 1, max_si)
    Lambda[t] <- sum(incidence[(t - 1):(t - smax)] * w[1:smax])
  }
  Lambda
}
