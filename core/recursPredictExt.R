######################################################################
## Bayesian recursive prediction via EpiFilter (grid-size-extended)
## Modified by: Raj Subedi, from the original recursPredict.R:
#
# From: Parag, KV, (2020) "Improved real-time estimation of reproduction numbers
# at low case incidence and between epidemic waves" BioRxiv.
######################################################################
#
## What changed from the original
## The original recursPredict.R hardcodes its internal prediction grid
## to Igrid = 0:800, and throws "Epidemic size too large" if any
## predicted mean incidence exceeds 0.9*800 = 720. That's fine for
## smaller outbreaks, but breaks outright on larger ones (e.g. an
## outbreak peaking above ~750 weekly cases). This version takes the
## grid's upper bound as a parameter (maxI) instead of a hardcoded
## constant, defaulting to 800 to match the original's behavior when
## not otherwise specified. Everything else - the Poisson renewal logic,
## the quantile extraction - is unchanged from the original.
#
# Notes and assumptions
# - observation model is Poisson renewal equation (as in EpiEstim)
# - reproduction number state space model is a simple diffusion
# - can apply causal or smoothing posteriors over R to predict incidence

# Inputs - grid on reproduction numbers (Rgrid), posterior on R (pR),
# total infectiousness (Lday), mean R esimate using pR (Rmean), confidence
# level (a and 50%), maxI (upper bound of the internal prediction grid -
# set this comfortably above your data's peak incidence)

# Output - mean prediction (pred) and confidence intervals (predInt)

recursPredictExt <- function(Rgrid, pR, Lday, Rmean, a, maxI = 800){

  # Grid size and length of time series
  nday = nrow(pR); m = ncol(pR)
  # Test lengths of inputs
  if (length(Rgrid) != m | length(Lday) != nday){
    stop("Input vectors of incorrect dimension")
  }

  # Mean prediction: Lday[i] => Iday[i+1]
  pred = Lday*Rmean; pred = pred[1:length(pred)-1]

  # Discrete space of possible predictions
  Igrid = 0:maxI; lenI = length(Igrid);

  # Check if close to upper bound
  if (any(pred > 0.9*max(Igrid))){
    stop(paste0("Epidemic size too large for maxI=", maxI,
                " - increase maxI (see recursPredictExt call in run_rtforecast.R)"))
  }

  # Prediction cdf and quantiles (50% and 95%)
  Fpred = matrix(0, nday-1, lenI)
  predInt = matrix(0, 4, nday-1)

  # At every time construct CDF of predictions
  for(i in 1:(nday-1)){
    # Compute rate from Poisson renewal
    rate = Lday[i]*Rgrid
    # Prob of any I marginalised over Rgrid
    pI = rep(0, lenI)

    # Probabilities of observations 1 day ahead
    for(j in 1:lenI){
      # Raw probabilities of Igrid
      pIset = dpois(Igrid[j], rate)
      # Normalised by probs of R
      pI[j] = sum(pIset*pR[i, ])
    }

    # Quantile predictions and CDF at i+1
    Fpred[i, ] = cumsum(pI)/sum(pI)
    id1 = which(Fpred[i, ] >= a); id2 = which(Fpred[i, ] >= 1-a)
    id3 = which(Fpred[i, ] >= 0.25); id4 = which(Fpred[i, ] >= 0.75)

    # Assign prediction results
    predInt[1, i] = Igrid[id1[1]]; predInt[2, i] = Igrid[id2[1]]
    predInt[3, i] = Igrid[id3[1]]; predInt[4, i] = Igrid[id4[1]]
  }
  # Main outputs: mean and 95% predictions
  recursPredictExt = list(pred, predInt)
}
