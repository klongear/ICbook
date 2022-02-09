### Gómez, Marhuenda & Langohr: "Regression with interval-censored covariates"
### 3.5 Illustration 3: Logistic regression (GLM)
### ===========================================================================
rm(list = ls(all = TRUE))
load("Datos.RData")

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
##   * Response variable: Obesity (i.e. BMI >= 30 kg/m^2)
##   * Uncensored covariate 1: Age
##   * Uncensored covariate 2: Energy intake
##   * Interval-censored covariate: Sum of carotenoides
## ---------------------------------------------------------------
library(dplyr)
library(tidyr)

dataF <- datos %>%
  filter(sexo == "Female") %>%
  drop_na(glucosa) %>%
  mutate(obesity = factor(imc >= 30, labels = c("No", "Yes")),
         obesnum = as.numeric(obesity) - 1,
         zl = round(totstart, 2), zr = round(totend, 2),
         zimp = (totstart + totend) / 2, energiat = energiat / 1000) %>%
  select(patient = paciente, obesity, obesnum, age, energy = energiat,
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
ymat <- matrix(rep(dataF$obesnum, m), nrow = n)
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
glm0 <- glm(obesnum ~ age + energy + zimp, dataF, family = binomial)
alfaHat <- glm0$coef[1]
bet1Hat <- glm0$coef[2]
bet2Hat <- glm0$coef[3]
gammHat <- glm0$coef[4]


## GEL (Gómez, Espinal, Lagakos) algorithm
## =======================================
repeat {
  # Estimation (update) of omega given theta
  # ----------------------------------------
  # Parameters of the gamma distribution in the first step
  linpred <- alfaHat + bet1Hat * cv1ma + bet2Hat * cv2ma + gammHat * smat
  repeat {
    omeOld <- omeHat
    omat <- matrix(rep(omeHat, n), byrow = TRUE, nrow = n)
    numerator <- alfas * dbinom(ymat, 1, exp(linpred) / (1 + exp(linpred))) * omat
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
  thetOld <- c(alfaOld, bet1Old, bet2Old, gammOld)
  lnfunc <- function(alf, bet1, bet2, gama) {
    -sum(log(rowSums(alfas *
                     exp(alf + bet1 * cv1ma + bet2 * cv2ma + gama * smat)^ymat /
                     (1 + exp(alf + bet1 * cv1ma + bet2 * cv2ma + gama * smat)) *
                     omat) + 0.00000001), na.rm = TRUE)
  }
  m0 <- mle2(lnfunc, start = list(alf = alfaHat, bet1 = bet1Hat, bet2 = bet2Hat,
                                  gama = gammHat),
                                  lower = c(alf = -Inf, bet1 = -Inf, bet2 = -Inf,
                                            gama = -Inf))#,
#                                   method = "L-BFGS-B")
  alfaHat <- m0@coef[1]
  bet1Hat <- m0@coef[2]
  bet2Hat <- m0@coef[3]
  gammHat <- m0@coef[4]
  thetHat <- c(alfaHat, bet1Hat, bet2Hat, gammHat)

  if (sum((omeHat - omeOld)^2) / sum(omeOld^2) + sum((thetHat - thetOld)^2) /
                                                 sum(thetOld^2) < accu)
    break
}

## Results of the algorithm
## ========================
underline(paste("Response variable:", nameY))

resus <- round(summary(m0)@coef, 4)
rownames(resus) <- c("Intercept", namesCovs, nameZ, "Error term")
cil <- round(resus[2:4, 1] + qnorm(c(0.025)) * resus[2:4, 2], 3)
ciu <- round(resus[2:4, 1] + qnorm(c(0.975)) * resus[2:4, 2], 3)
cis <- cbind(resus[2:4, 1], cil, ciu, resus[2:4, 4])
colnames(cis) <- c("Estimate", "Lower", "Upper", "p-value")
{
  underline("Interval-censored covariate", "-")
  cat("\nParameter estimates:\n")
  print(resus)
  cat("\n95% Confidence intervals:\n")
  print(cis)
}


## Results of model fit with midpoint imputation
## =============================================
glm0 <- glm(obesnum ~ age + energy + zimp, dataF, family = binomial)
resusI <- round(summary(glm0)$coef, 4)
rownames(resusI) <- c("Intercept", namesCovs, nameZ)
cisI <- round(cbind(resusI[2:4, 1], confint(glmodI)[2:4, ], resusI[2:4, 4]), 4)
colnames(cisI) <- c("Estimate", "Lower", "Upper", "p-value")
{
  underline("Midpoint estimation", "-")
  cat("\nParameter estimates:\n")
  print(resusI)
  cat("\n95% Confidence intervals:\n")
  print(cisI)
}
