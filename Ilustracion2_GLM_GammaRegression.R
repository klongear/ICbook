### Gómez, Marhuenda & Langohr: "Regression with interval-censored covariates"
### 3.5 Illustration 2: Gamma regression (GLM)
### ===========================================================================
rm(list = ls(all = TRUE))
load("Illustration.RData")

## Auxiliar function for nice output
## =================================
underline <- function(txt, symb = "=") {
  cat("\n", txt, "\n", sep = "")
  cat(rep(symb, nchar(txt)), "\n", sep = "")
}


### GEL algorithm to estimate the parameters of a linear regression model
### with an interval-censored covariate.
### =====================================================================
## Package needed for parameter estimation
## ---------------------------------------
library(bbmle)

## Convergence tolerance (the smaller, the longer takes the algorithm)
## -------------------------------------------------------------------
accu <- 1e-8

## Variables included in the model with data of female volunteers:
##   * Response variable: Glucose
##   * Uncensored covariate 1: Age
##   * Uncensored covariate 2: Energy intake
##   * Interval-censored covariate: Sum of carotenoides
## ---------------------------------------------------------------
library(dplyr)
library(tidyr)

dataF <- datos %>%
  filter(sexo == "Female") %>%
  mutate(zl = round(totstart, 2), zr = round(totend, 2),
         zimp = (totstart + totend) / 2, energiat = energiat / 1000) %>%
  drop_na(glucosa) %>%
  select(patient = paciente, glucose = glucosa, age, energy = energiat,
         zl, zr, zimp)

## Number of observations
## ----------------------
n <- nrow(dataF)

## Support of interval-censored covariate
## --------------------------------------
sj <- round(seq(min(dataF$zl), max(dataF$zr), 0.01), 2)

## Length of the support
## ---------------------
m <- length(sj)

## Auxiliary matrices
## ------------------
ymat <- matrix(rep(dataF$glucose, m), nrow = n)
cv1ma <- matrix(rep(dataF$age, m), nrow = n)
cv2ma <- matrix(rep(dataF$energy, m), nrow = n)
zlma <- matrix(rep(dataF$zl, m), nrow = n)
zrma <- matrix(rep(dataF$zr, m), nrow = n)

## Matrix of the indicator variables
## ---------------------------------
smat <- matrix(rep(sj, n), byrow = TRUE, nrow = n)
alfas <- (zlma <=  smat) * (smat <=  zrma)

## Initiation of the algorithm
## ---------------------------
## 1. Omega [Distribution function of interval-censored covariate]
omeHat <- rep(1 / m, m)
omat <- matrix(rep(omeHat, n), byrow = TRUE, nrow = n)

## 2. Model parameters
glm0 <- glm(glucose ~ age + energy + zimp, dataF, family = Gamma("log"))
alfaHat <- glm0$coef[1]
bet1Hat <- glm0$coef[2]
bet2Hat <- glm0$coef[3]
gammHat <- glm0$coef[4]
tauHat <- summary(glm0)$dispersion
summary(glm0)

## GEL (Gómez, Espinal, Lagakos) algorithm
## =======================================
repeat {
  # Estimation (update) of omega given theta
  # ----------------------------------------
  # Parameters of the gamma distribution in the first step
  scalePm <- exp(alfaHat + bet1Hat * cv1ma + bet2Hat * cv2ma
                         + gammHat * smat) / tauHat
  shapePm <- tauHat
  repeat {
    omeOld <- omeHat
    omat <- matrix(rep(omeHat, n), byrow = TRUE, nrow = n)
    numerator <- alfas * dgamma(ymat, shape = shapePm, scale = scalePm) * omat
    denominator <- matrix(rep(rowSums(numerator), m), nrow = n)
    nuumat <- numerator / denominator
    omeHat <- colSums(nuumat, na.rm = TRUE) / n
    if (sum((omeHat - omeOld)^2) / sum(omeOld^2) < accu)
      break
    }
  omat <- matrix(rep(omeHat, n), byrow = TRUE, nrow = n)

  # Estimation (update) of model parameters given omega
  # ---------------------------------------------------
  alfaOld <- alfaHat
  bet1Old <- bet1Hat
  bet2Old <- bet2Hat
  gammOld <- gammHat
  tauOld <- tauHat
  thetOld <- c(alfaOld, bet1Old, bet2Old, gammOld, tauOld)
  lnfunc <- function(alf, bet1, bet2, gama, tau) {
    -sum(log(rowSums(alfas *
                     dgamma(ymat, shape = tau,
                            scale = exp(alf + bet1 * cv1ma + bet2 * cv2ma +
                                        gama * smat) / tau) * omat) + 0.00000001),
         na.rm = TRUE)
  }
  m0 <- mle2(lnfunc, start = list(alf = alfaHat, bet1 = bet1Hat, bet2 = bet2Hat,
                                  gama = gammHat, tau = tauHat),
                                  lower = c(alf = -Inf, bet1 = -Inf, bet2 = -Inf,
                                            gama = -Inf, tau = 0),
                                  method = "L-BFGS-B")
  alfaHat <- m0@coef[1]
  bet1Hat <- m0@coef[2]
  bet2Hat <- m0@coef[3]
  gammHat <- m0@coef[4]
  tauHat <- m0@coef[5]
  thetHat <- c(alfaHat, bet1Hat, bet2Hat, gammHat, tauHat)

  if (sum((omeHat - omeOld)^2) / sum(omeOld^2) + sum((thetHat - thetOld)^2) /
                                                 sum(thetOld^2) < accu)
    break
}

## Results of the algorithm
## ========================
underline("Gamma regression (using log link) for response variable Glucose:")

resus <- round(summary(m0)@coef, 4)[1:4, ]
rownames(resus) <- c("Intercept", "Age", "Energy intake", "Sum carotenoides")
cil <- round(resus[2:4, 1] + qnorm(c(0.025)) * resus[2:4, 2], 4)
ciu <- round(resus[2:4, 1] + qnorm(c(0.975)) * resus[2:4, 2], 4)
cis <- cbind(resus[2:4, 1], cil, ciu, resus[2:4, 4])
colnames(cis) <- c("Estimate", "Lower", "Upper", "p-value")
{
  underline("Interval-censored covariate", "-")
  cat("\nParameter estimates:\n")
  print(resus)
  cat("\n95% Confidence intervals:\n")
  print(cis)
}

## Cleaning up
## -----------
rm(accu, m0, omeOld, numerator, denominator,
   nuumat, zlma, zrma, ymat, thetOld, lnfunc, cil,
   ciu, glm0, m, n, sj, alfaOld, bet1Old, bet2Old, cv1ma,
   cv2ma, gammOld, tauOld)


## Computation of Pearson residuals
## ================================
source("RFunction_onlyOmega.R")
omhat <- onlyOmega(dataF$zl, dataF$zr)
omat <- matrix(rep(omhat[, 2], n), byrow = TRUE, nrow = n)

dataF <- mutate(dataF,
                Zexp = rowSums(alfas * smat * omat) / rowSums(alfas * omat),
                linPred = alfaHat + bet1Hat * age + bet2Hat * energy +
                          gammHat * Zexp,
                muhat = exp(linPred),
                varhat = exp(linPred)^2 / tauHat,
                peares = (glucose - muhat) / sqrt(varhat)) %>%
  select(-linPred)

## Deviance residuals (see Czado, p. 18)
devi <- with(dataF, -2 * sum(log(glucose / muhat) - (glucose - muhat) / muhat))

## p-value
1 - pchisq(devi * tauHat, n - 3)
