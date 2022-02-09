### Gómez, Marhuenda & Langohr: "Regression with interval-censored covariates"
### 3.5 Illustration 1: Linear regression with log-transformed response
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
  mutate(obesity = factor(imc >= 30, labels = c("No", "Yes")),
         logGlucose = log(glucosa), zl = round(totstart, 2),
         zr = round(totend, 2), zimp = (totstart + totend) / 2,
         energy = energiat/1000) %>%
  drop_na(glucosa) %>%
  select(patient = paciente, glucose = glucosa, logGlucose, age,
         energy, obesity, zl, zr, zimp)

## Descriptive analysis and GofCens figures
## ========================================
# library(compareGroups)
# export2latex(descrTable(~ glucose + logGlucose + age + energy + zl + zr + zimp,
#                         dataF, extra.labels = c("", "", "", ""), method = 1),
#              label = "tab:desi",
#              file = "DesiMean.tex")
#
# export2latex(descrTable(~ glucose + logGlucose + age + energy + obesity + zl + zr + zimp,
#                         dataF, extra.labels = c("", "", "", ""), method = 2),
#              file = "DesiMedian.tex")
#
# library(GofCens)
# cumhazPlot(dataF$glucose, ggplo = TRUE, prnt = FALSE)

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
ymat <- matrix(rep(dataF$logGlucose, m), nrow = n)
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
lm0 <- lm(logGlucose ~ age + energy + zimp, dataF)
alfaHat <- lm0$coef[1]
bet1Hat <- lm0$coef[2]
bet2Hat <- lm0$coef[3]
gammHat <- lm0$coef[4]
sigmHat <- summary(lm0)$sigma


## GEL (Gómez, Espinal, Lagakos) algorithm
## =======================================
repeat {
  # Updated estimation of omega
  # ---------------------------
  repeat {
    omeOld <- omeHat
    omat <- matrix(rep(omeHat, n), byrow = TRUE, nrow = n)
    numerator <- alfas * dnorm(ymat, alfaHat + bet1Hat * cv1ma +
                                     bet2Hat * cv2ma + gammHat * smat,
                                     sigmHat) * omat
    denominator <- matrix(rep(rowSums(numerator), m), nrow = n)
    nuumat <- numerator / denominator
    omeHat <- colSums(nuumat, na.rm = TRUE) / n
    if (sum((omeHat - omeOld)^2) / sum(omeOld^2) < accu)
      break
    }

  # Updated estimation of model parameters
  # --------------------------------------
  alfaOld <- alfaHat
  bet1Old <- bet1Hat
  bet2Old <- bet2Hat
  gammOld <- gammHat
  sigmOld <- sigmHat
  thetOld <- c(alfaOld, bet1Old, bet2Old, gammOld, sigmOld)
  lnfunc <- function(alf, bet1, bet2, gama, sig) {
    -sum(log(rowSums(alfas * dnorm(ymat, alf + bet1 * cv1ma + bet2 * cv2ma +
                                         gama * smat, sig) *
                                   omat) + 0.00000001), na.rm = TRUE)
  }

  m0 <- mle2(lnfunc, start = list(alf = alfaHat, bet1 = bet1Hat, bet2 = bet2Hat,
                                  gama = gammHat, sig = sigmHat),
                                  lower = c(alf = -Inf, bet1 = -Inf, bet2 = -Inf,
                                            gama = -Inf, sig = 0),
                                  method = "L-BFGS-B")
  alfaHat <- m0@coef[1]
  bet1Hat <- m0@coef[2]
  bet2Hat <- m0@coef[3]
  gammHat <- m0@coef[4]
  sigmHat <- m0@coef[5]
  thetHat <- c(alfaHat, bet1Hat, bet2Hat, gammHat, sigmHat)

  if (sum((omeHat - omeOld)^2) / sum(omeOld^2) + sum((thetHat - thetOld)^2) /
                                                 sum(thetOld^2) < accu)
    break
}

## Results of the algorithm
## ========================
underline("Response variable: log(Glucose)")

resus <- round(summary(m0)@coef, 5)
rownames(resus) <- c("Intercept", "Age", "Energy intake", "Sum carotenoides",
                     "Error term")
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
modI <- lm(logGlucose ~ age + energy + zimp, dataF)
resusI <- round(summary(modI)$coef, 5)
rownames(resusI) <- c("Intercept", "Age", "Energy intake", "Sum carotenoides")
cisI <- round(cbind(resusI[2:4, 1], confint(modI)[2:4, ], resusI[2:4, 4]), 6)
colnames(cisI) <- c("Estimate", "Lower", "Upper", "p-value")
{
  underline("Midpoint estimation", "-")
  cat("\nParameter estimates:\n")
  print(resusI)
  cat("\n95% Confidence intervals:\n")
  print(cisI)
}


## Cleaning up
## -----------
rm(accu, m0, omeOld, numerator, denominator,
   nuumat, zlma, zrma, ymat, thetOld, lnfunc, omeHat, cil,
   ciu, lm0, m, n, omat, sj, alfaOld, bet1Old, bet2Old, cv1ma,
   cv2ma, gammOld, sigmOld)


## Goodness of fit with GEL and Topp-Gómez residuals
## =================================================
source("RFunctions_Residuals_2Covs.R")

## GEL residuals
## -------------
omhat <- onlyOmega(dataF$zl, dataF$zr)
gelresi <- with(dataF, resiGEL(logGlucose, age, energy, zl, zr, alfaHat, bet1Hat,
                               bet2Hat, gammHat, omhat))

## Normality plots
## ===============
library(nortest)

pdf("LMlogY/NormalityLogY.pdf", width = 9, height = 9)
par(font = 2, font.lab = 2, font.axis = 2, las = 1,
    cex.main = 2.1, cex.lab = 2, cex.axis = 1.9, mar  =  c(5, 4.5, 4, 2))
qqnorm(gelresi, pch = 19, main = "")
text(-2, 0.6, paste("KS test: p =", round(lillie.test(gelresi)$p.val, 3)),
     pos = 4, cex = 1.75)
qqline(unclass(gelresi), lwd = 2)
dev.off()

## Plots to check homoscesdaticity
## ===============================
library(lmtest)
# zethat <- with(dataF, (zl + zr) / 2)
fittedGEL <- dataF$logGlucose - gelresi
Zexp <- round(rowSums(alfas * smat * omat) / rowSums(alfas * omat), 3)

pdf("LMlogY/HomoscedasticityLogY.pdf", width = 9, height = 9)
par(font = 2, font.lab = 2, font.axis = 2, las = 1,
    cex.main = 2.1, cex.lab = 2, cex.axis = 1.9, mar  =  c(5, 4.5, 4, 2))
plot(gelresi ~ fittedGEL, pch = 19, main = "", xlab = "Fitted values",
     ylab = "Residuals")
text(4.55, 0.6, paste("BP test: p = ",
                      round(bptest(lm(gelresi ~ age + energy + Zexp,
                                      dataF))$p.val, 3)),
     pos = 4, cex = 1.75)
abline(h = 0, lwd = 2)
dev.off()

## Both plots in one pdf file
## --------------------------
pdf("LMlogY/ResidualPlotsLogY.pdf", width = 12, height = 6)
par(mfrow = 1:2, font = 2, font.lab = 2, font.axis = 2, las = 1,
    cex.lab = 1.2, cex.axis = 1.2, mar  =  c(5, 4.5, 2, 2))
qqnorm(gelresi, pch = 19, main = "")
qqline(unclass(gelresi), lwd = 2)
legend("topleft", paste("KS test: p =", round(lillie.test(gelresi)$p.val, 3)),
       bty = "n")

plot(gelresi ~ fittedGEL, pch = 19, main = "", xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0, lwd = 2)
legend("topleft", paste("BP test: p = ",
                        round(bptest(lm(gelresi ~ age + energy + Zexp,
                                        dataF))$p.val, 3)),
       bty = "n")
dev.off()


## Scatterplot
## ===========
pdf("LMlogY/Scatter.pdf", width = 9, height = 7)
par(font = 2, font.lab = 2, font.axis = 2, las = 1,
    cex.lab = 1.2, cex.axis = 1.2, mar  =  c(5, 4.5, 2, 2))
plot(glucose ~ zimp, dataF, pch = 16,
     xlab = expression(bold(paste("Total carotenoid concentration [", mu, "mol/L]",
                                  sep = ""))),
     ylab = "Glucose [mg/dl]")
dev.off()
