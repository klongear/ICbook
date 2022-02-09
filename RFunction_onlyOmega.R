### ==================================================
### Function onlyOmega computes the Turnbull estimator
### for an interval-censored variable.
### Note: Closed intervals [Zl, Zr] are assumed.
### Arguments:
###   zl: Vector containing the lower limits
###   zr: Vector containing the upper limits
###   toler: Convergence tolerance
### ==================================================

onlyOmega <- function(zl, zr, toler = 1e-11) {
  ll <- length(zl)
  lr <- length(zr)
  if (min(ll, lr) != max(ll, lr)) {
    stop("Vectors must be of equal length!")
  }
  if (sum(zl > zr) > 0) {
    stop("Zl must be smaller than Zr!")
  }
  n <- length(zl)
  # Possible values of Z
  sj <- round(seq(min(zl), max(zr), 0.01), 2)
  # Length of support of Z
  m <- length(sj)
  # Auxiliary matrices of the data
  zlma <- matrix(rep(zl, m), nrow = n)
  zrma <- matrix(rep(zr, m), nrow = n)
  # Matrix of indicator variables
  smat <- matrix(rep(sj, n), byrow = TRUE, nrow = n)
  alfas <- (zlma <= smat) * (smat <= zrma)

  ## Initial values
  # Omega
  omeHat <- rep(1 / m, m)
  omat <- matrix(rep(omeHat, n), byrow = TRUE, nrow = n)

  # TB algorithm
  repeat {
    omeOld <- omeHat
    omat <- matrix(rep(omeHat, n), byrow = TRUE, nrow = n)
    numerator <- alfas * omat
    denominator <- matrix(rep(rowSums(numerator), m), nrow = n)
    nuumat <- numerator / denominator
    omeHat <- colSums(nuumat) / n

    if (sum((omeHat - omeOld)^2) / sum(omeOld^2) < toler)
      break
  }
  # Returning the estimated probabilities
  cbi <- cbind(sj, omegaHat = round(omeHat, 5))
#  ok <- cbi[, 2] > 0
#  return(cbi[ok, ])
  return(cbi)
}
