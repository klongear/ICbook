### 1. Function GELalgo fits a linear regression model
### with an interval-censored covariate Z (and an additional covariate X)
### using the GEL algorithm (Gomez, Espinal and Lagakos (2003))
### Arguments:
###   y: Vector containing the values of the response variable
###   x: Vector containing the values of the additional covariate
###   zl: Vector containing the lower limits of the interval-censored covariate Z
###   zr: Vector containing the upper limits of the interval-censored covariate Z
###   toler: Convergence tolerance
#################################################################################

GELalgo <- function(y, x, zl, zr, toler = 1e-11){
  ly <- length(y); lx <- length(x); ll <- length(zl); lr <- length(zr)

  if (min(ly, lx, ll, lr) != max(ly, lx, ll, lr))
    stop('Vectors must be of equal length!')
  if (any(zl>=zr))
    stop('Zl must be smaller than Zr!')

  # Number of observations
  n <- length(y)
  # Vector of Z's possible values
  sj <- min(zl):max(zr)
  # Length of Z's support
  m <- length(sj)
  # Auxiliary matrices of the data
  ymat <- matrix(rep(y, m), nrow = n)
  xmat <- matrix(rep(x, m), nrow = n)
  zlma <- matrix(rep(zl, m), nrow = n)
  zrma <- matrix(rep(zr, m), nrow = n)
  # Matrix of indicator variables
  smat <- matrix(rep(sj, n), byrow = TRUE, nrow = n)
  alfas <- (zlma<=smat)*(smat<=zrma)

  ## Initial values
  # Omega
  omeHat <- rep(1/m, m)
  omat <- matrix(rep(omeHat, n), byrow = TRUE, nrow = n)
  # Theta
  zet <- 0.5*(zl+zr)
  lm0 <- lm(y~x+zet)
  alfHat <- lm0$coef[1]
  betHat <- lm0$coef[2]
  gamHat <- lm0$coef[3]
  sigHat <- summary(lm0)$sigma

  # GEL algorithm
  repeat{
    # Updated estimation of omega
    repeat{
      omeOld <- omeHat
      omat <- matrix(rep(omeHat, n), byrow = TRUE, nrow = n)
      numerator <- alfas*dnorm(ymat, alfHat+betHat*xmat+gamHat*smat, sigHat)*omat
      denominator <- matrix(rep(rowSums(numerator), m), nrow = n)
      nuumat <- numerator/denominator
      omeHat <- colSums(nuumat)/n

      if (sum((omeHat-omeOld)^2)/sum(omeOld^2)<toler)
        break
    }

    # Updated estimation of Theta
    alfOld <- alfHat
    betOld <- betHat
    gamOld <- gamHat
    sigOld <- sigHat
    thetOld <- c(alfOld, betOld, gamOld, sigOld)
    lnfunc <- function(alf, bet, gam, sig){
      -sum(log(rowSums(alfas*dnorm(ymat, alf+bet*xmat+gam*smat, sig)*omat) +.00000001))
    }
    # Updated estimation of parameters
    m0 <- mle2(lnfunc, start = list(alf = alfHat, bet = betHat, gam = gamHat, sig = sigHat), lower = c(alf = -Inf, bet = -Inf, gam = -Inf, sig = 0), method = "L-BFGS-B")
    sm0 <- summary(m0)@coef
    alfHat <- sm0[1, 1]
    betHat <- sm0[2, 1]
    gamHat <- sm0[3, 1]
    sigHat <- sm0[4, 1]
    thetHat <- c(alfHat, betHat, gamHat, sigHat)

    if (sum((omeHat-omeOld)^2)/sum(omeOld^2) + sum((thetHat-thetOld)^2)/sum(thetOld^2) < toler)
      break
  }
  # Distribution function of Z
  cbi <- cbind(sj, omegaHat = round(omeHat, 4), FZHat = round(cumsum(omeHat), 4))
  ok <- cbi[, 2]>0
  return(list(Theta = round(thetHat, 4), Omega = cbi[ok, ]))
}



###·===========================================================================
### 2. Function resiMidp computes the Midpoint residuals
### Arguments:
###   y: Response variable
###   x1: First uncensored covariate
###   x2: Second uncensored covariate
###   zl: Lower limits of the interval-censored covariate Z
###   zr: Upper limits of the interval-censored covariate Z
###   alfhat: Estimated model constant
###   b1hat: Estimated parameter of first covariate
###   b2hat: Estimated parameter of second covariate
###   gamhat: Estimated parameter of the interval censored covariate
###   dec: Number of decimal digits
###·===========================================================================

resiMidp <- function(y, x1, x2, zl, zr, alfhat, b1hat, b2hat, gamhat, dec = 3) {
#  auxdf <- data.frame(y, x1, x2, x3, zl, zr)
  if (any(zl > zr)) {
    stop("Zl must not be larger than Zr!")
  }
  zet <- 0.5 * (zl + zr)
  resiMid <- round(y - (alfhat + b1hat * x1 + b2hat * x2 + gamhat * zet), dec)
  return(resiMid)
}

#################################################################################
### 3. Function resiGEL computes the GEL residuals
### Arguments:
###   y: Response variable
###   x1: First uncensored covariate
###   x2: Second uncensored covariate
###   zl: Lower limits of the interval-censored covariate Z
###   zr: Upper limits of the interval-censored covariate Z
###   alfhat: Estimated model constant
###   b1hat: Estimated parameter of first covariate
###   b2hat: Estimated parameter of second covariate
###   gamhat: Estimated parameter of the interval censored covariate
###   omehat: Matrix containing the Turnbull estimate of F_Z
###   dec: Number of decimal digits
#################################################################################

resiGEL <- function(y, x1, x2, zl, zr, alfhat, b1hat, b2hat, gamhat, omehat,
                    dec = 3){
  if (any(zl > zr)) {
    stop("Zl must not be larger than Zr!")
  }
  n <- length(y)
  # Possible values of Z
  sj <- omehat[, 1]
  # Length of support of Z
  m <- length(sj)
  # Auxiliary matrices
  zlma <- matrix(rep(zl, m), nrow = n)
  zrma <- matrix(rep(zr, m), nrow = n)
  omat <- matrix(rep(omehat[, 2], n), byrow = TRUE, nrow = n)
  smat <- matrix(rep(sj, n), byrow = TRUE, nrow = n)
  # Matrix of indicator variables
  alfas <- (zlma <= smat) * (smat <= zrma)
  # Computation of E(Z|Zl, Zr; TB)
  Zexp <- round(rowSums(alfas * smat * omat) / rowSums(alfas * omat), dec)

  resiGel <- round(y - (alfhat + b1hat * x1 + b2hat * x2 + gamhat * Zexp), dec)
  return(resiGel)
}

#################################################################################
### 4. Function resi.ToGo computes the Topp-Gomez residuals
### Arguments:
###   y: Response variable
###   x1: First uncensored covariate
###   x2: Second uncensored covariate
###   zl: Lower limits of the interval-censored covariate Z
###   zr: Upper limits of the interval-censored covariate Z
###   alfhat: Estimated model constant
###   b1hat: Estimated parameter of first covariate
###   b2hat: Estimated parameter of second covariate
###   gamhat: Estimated parameter of the interval censored covariate
###   sighat: Estimated residual standard deviation
###   dec: Number of decimal digits
#################################################################################

resiToGo <- function(y, x1, x2, zl, zr, alfhat, b1hat, b2hat, gamhat, sighat,
                     dec = 3){
  if (any(zl > zr)) {
    stop("Zl must not be larger than Zr!")
  }

  # Computation of A_i and B_i
  ll <- y - (alfhat + b1hat * x1 + b2hat * x2 + gamhat * zl)
  rr <- y - (alfhat + b1hat * x1 + b2hat * x2 + gamhat * zr)

  Ai <- ifelse(ll <= rr, ll, rr)
  Bi <- ifelse(ll <= rr, rr, ll)
  resiTogo <- round((dnorm(Ai / sighat) - dnorm(Bi / sighat)) /
                    (pnorm(Bi / sighat) - pnorm(Ai / sighat)) * sighat, dec)
  return(resiTogo)
}

#########################################################
### 5. Function onlyOmega computes the Turnbull estimator
### for an interval-censored variable.
### Note: Closed intervals [Zl, Zr] are assumed.
### Arguments:
###   zl: Vector containing the lower limits
###   zr: Vector containing the upper limits
###   toler: Convergence tolerance
#########################################################

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
#  ok <- cbi[, 2]>0
#  return(cbi[ok, ])
  return(cbi)
}
