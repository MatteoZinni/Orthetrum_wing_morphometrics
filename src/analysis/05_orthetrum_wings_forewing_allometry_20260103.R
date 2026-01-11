## ========================================================================== ##
## Script:      05_orthetrum_wings_forewing_allometry.R
## Author:      Matteo Zinni
## Date:        2025-01-03
## Description: Allometric analysis of Orthetrum forewing shape using
##              geometric morphometrics. The script investigates centroid
##              size variation and size–shape relationships across species
##              and sexes, using full datasets and outlier-filtered subsets.
## ========================================================================== ##

# 05.01.01.01 FOREWING ANALYSIS ------------------------------------------------

# Packages loading ----
message("Loading packages: ", paste(packages, collapse = ", "))
invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

## 05.01.01.01 CENTROID SIZE ---------------------------------------------------

## 05.01.02.01 Full dataset ----------------------------------------------------

# Distribution of centroid size values
hist(forewing_2d$C_size,
     main = "Forewing centroid size distribution (forewing)",
     xlab = "Centroid size")

#### Does centroid size differ among species? ----

# Testing for difference (chi-squared = 34.641, df = 9, p-value = 6.893e-05)
kruskal.test(C_size ~ Species, data = forewing_2d)

# Post-hoc test for pairwise comparison
kwAllPairsDunnTest(C_size ~ Species, 
                   data = forewing_2d,
                   p.adjust.method = "bonferroni"
                   )

# Plot results
boxplot(C_size ~ Species, data = forewing_2d,
        ylab = "Centroid size",
        main = "Forewing centroid size among species \n (full dataset)")

#### Does centroid size differ between sexes? ----

# Testing for difference (chi-squared = 0.0269, df = 1, p-value = 0.8696)
kruskal.test(C_size ~ Sex, data = forewing_2d)

# Post-hoc test for pairwise comparison
kwAllPairsDunnTest(C_size ~ Sex, 
                   data = forewing_2d,
                   p.adjust.method = "bonferroni"
                   )

# Plot results
boxplot(C_size ~ Sex, data = forewing_2d,
        ylab = "Centroid size",
        main = "Forewing centroid size between sexes \n (full dataset)")

# Distribution of centroid size values between sexes 
hist(forewing_2d$C_size[forewing_2d$Sex == 'Male'], 
     col='deepskyblue',
     main='Forewing centroid size between sexes \n Full dataset', 
     xlab='Centroid size')
hist(forewing_2d$C_size[forewing_2d$Sex == 'Female'], 
     col='darkgoldenrod',
     add=TRUE)

## 05.01.03.01 Dataset without species outliers --------------------------------

# Distribution of centroid size values
hist(fwSpOut_2d$C_size)

#### Does centroid size differ among species? ----

# Testing for difference (chi-squared = 33.381, df = 9, p-value = 0.0001146)
kruskal.test(C_size ~ Species, data = fwSpOut_2d)

# Post-hoc test for pairwise comparison
kwAllPairsDunnTest(C_size ~ Species, data = fwSpOut_2d,
                   p.adjust.method = "bonferroni")

# Plot results
boxplot(C_size ~ Species, data = fwSpOut_2d,
        ylab = "Centroid size",
        main = "Centroid size among species")

## 05.01.04.01 Dataset without sex outliers ------------------------------------

# Distribution of centroid size values
hist(fwSxOut_2d$C_size)

#### Does centroid size differ between sexes? ----

# Testing for difference (chi-squared = 0.10185, df = 1, p-value = 0.7496)
kruskal.test(C_size ~ Sex, data = fwSxOut_2d)

# Plot results
boxplot(C_size ~ Sex, data = fwSxOut_2d,
        ylab = "Centroid size",
        main = "Centroid size between sexes")

# Distribution of centroid size values between sexes 
hist(fwSxOut_2d$C_size[fwSxOut_2d$Sex == 'Male'], 
     col='deepskyblue',
     main='Centroid size between sexes', xlab='Centroid size')
hist(fwSxOut_2d$C_size[fwSxOut_2d$Sex == 'Female'], 
     col='darkgoldenrod',
     add=TRUE)

## 05.02.01.01 ALLOMETRY ANALYSIS ----------------------------------------------

## 05.02.02.01 Full dataset ----------------------------------------------------

#### Effect of size on shape ----

# Null model  
fit01_00 <-  procD.lm(coords ~ 1, 
                      RRPP=T, 
                      data = forewing_gdf, 
                      iter = 9999,
                      seed = 123)

# model summary 
summary(fit01_00)

# Simple size on shape covariation model
fit01_01 <- procD.lm(coords~log(Csize),
                     data = forewing_gdf, 
                     RRPP=T,
                     iter = 9999, 
                     seed = 123)

# model summary 
summary(fit01_01)

#### Species-specific allometric trajectories (unique slopes) ----

# Species allometry | full dataset | unique slopes
fit01_02 <- procD.lm(coords~log(Csize)*Species,
                     data = forewing_gdf,
                     SS.type = "I",
                     effect.type = "F",
                     RRPP = TRUE,
                     iter = 9999, 
                     seed = 123)

# model summary 
summary(fit01_02)

#### Common allometric slope with species-specific shape intercepts ----

# Species allometry | full dataset | common slopes
fit01_03 <- procD.lm(coords~log(Csize) + Species,
                     data = forewing_gdf,
                     SS.type = "I",
                     effect.type = "F",
                     RRPP = TRUE,
                     iter = 9999, 
                     seed = 123)

# model summary 
summary(fit01_03)

#### Model comparison and selection ----

# Comparing models with anova: fit01_02 (unique) is the best candidate
anova(fit01_00, fit01_01, fit01_02, fit01_03)

## 05.02.02.02 Pairwise comparison of allometric trajectories ----

# Compute pairwise comparisons of unique allometric slopes among species
fwFd_pw <- pairwise(fit01_02, 
                       groups = forewing_gdf$Species, 
                       covariate = log(forewing_gdf$Csize))

# Differences in slope vector length
summary(fwFd_pw, confidence = 0.95, test.type = "dist", stat.table = FALSE)

# Differences in slope vector correlation (VC = vector correlation)
summary(fwFd_pw, confidence = 0.95, test.type = "VC", stat.table = FALSE)

# Differences in directional vector (DL = direction length)
summary(fwFd_pw, confidence = 0.95, test.type = "DL", stat.table = FALSE)

## 05.02.03.01 Dataset without species outliers --------------------------------

#### Effect of size on shape ----

# Null model  
fit02_00 <-  procD.lm(coords ~ 1, 
                      RRPP=T, 
                      data = fwSpOut_gdf, 
                      iter = 9999,
                      seed = 123)

# model summary 
summary(fit02_00)

# Simple size on shape covariation model
fit02_01 <- procD.lm(coords~log(Csize),
                     data = fwSpOut_gdf, 
                     RRPP=T,
                     iter = 9999, 
                     seed = 123)

# model summary 
summary(fit02_01)

#### Species-specific allometric trajectories (unique slopes) ----

# Species allometry | no species outliers | unique slopes
fit02_02 <- procD.lm(coords~log(Csize)*Species,
                     data = fwSpOut_gdf,
                     SS.type = "I",
                     effect.type = "F",
                     RRPP = TRUE,
                     iter = 9999, 
                     seed = 123)

# model summary 
summary(fit02_02)

#### Common allometric slope with species-specific shape intercepts ----

# Species allometry | no species outliers | common slopes
fit02_03 <- procD.lm(coords~log(Csize)+ Species,
                     data = fwSpOut_gdf,
                     SS.type = "I",
                     effect.type = "F",
                     RRPP = TRUE,
                     iter = 9999, 
                     seed = 123)

# model summary 
summary(fit02_03)

#### Model comparison and selection ----

# Comparing models with anova: fit02_02 (unique) is the best candidate
anova(fit02_00, fit02_01, fit02_02, fit02_03)

# Pairwise comparison
fwSpOut_pw = pairwise(fit02_02, 
                      groups = fwSpOut_gdf$Species, 
                      covariate = log(fwSpOut_gdf$Csize))

## 05.02.03.02 Pairwise comparison of allometric trajectories ----

# Compute pairwise comparisons of unique allometric slopes among species
fwSpOut_pw <- pairwise(fit02_02, 
                       groups = fwSpOut_gdf$Species, 
                       covariate = log(fwSpOut_gdf$Csize))

# Differences in slope vector length
summary(fwSpOut_pw, confidence = 0.95, test.type = "dist", stat.table = FALSE)

# Differences in slope vector correlation (VC = vector correlation)
summary(fwSpOut_pw, confidence = 0.95, test.type = "VC", stat.table = FALSE)

# Differences in directional vector (DL = direction length)
summary(fwSpOut_pw, confidence = 0.95, test.type = "DL", stat.table = FALSE)

## 05.02.04.01 Dataset without sex outliers ------------------------------------

#### Effect of size on shape ----

# Null model  
fit03_00 <-  procD.lm(coords ~ 1, 
                      RRPP = TRUE, 
                      data = fwSxOut_gdf, 
                      iter = 9999,
                      seed = 123)

# Model summary 
summary(fit03_00)

#### Species-specific allometric trajectories (unique slopes) ----

# Species allometry | no sex outliers | unique slopes
fit03_01 <- procD.lm(coords ~ log(Csize),
                     data = fwSxOut_gdf, 
                     RRPP = TRUE,
                     iter = 9999, 
                     seed = 123)

# Model summary 
summary(fit03_01)

#### Common allometric slope with species-specific shape intercepts ----

# Species allometry | no sex outliers | common slopes
fit03_02 <- procD.lm(coords ~ log(Csize)*Sex,
                     data = fwSxOut_gdf,
                     SS.type = "I",
                     effect.type = "F",
                     RRPP = TRUE,
                     iter = 9999, 
                     seed = 123)

# Model summary 
summary(fit03_02)

# Common slope model 
fit03_03 <- procD.lm(coords ~ log(Csize) + Sex,
                     data = fwSxOut_gdf,
                     SS.type = "I",
                     effect.type = "F",
                     RRPP = TRUE,
                     iter = 9999, 
                     seed = 123)

# Model summary 
summary(fit03_03)

#### Model comparison and selection ----

# Comparing models with anova: fit03_02 (unique) is still the best candidate
anova(fit03_01, fit03_02, fit03_03)

## 05.02.04.02 Pairwise comparison of allometric trajectories ----

# Compute pairwise comparisons of unique allometric slopes among species
fwSpOut_pw <- pairwise(fit03_02, 
                       groups = fwSpOut_gdf$Species, 
                       covariate = log(fwSpOut_gdf$Csize))

# Differences in slope vector length
summary(fwSpOut_pw, confidence = 0.95, test.type = "dist", stat.table = FALSE)

# Differences in slope vector correlation (VC = vector correlation)
summary(fwSpOut_pw, confidence = 0.95, test.type = "VC", stat.table = FALSE)

# Differences in directional vector (DL = direction length)
summary(fwSpOut_pw, confidence = 0.95, test.type = "DL", stat.table = FALSE)

## 04.05.01.01 PLOT AND DATA VISUALIZATION -------------------------------------

### 04.05.01.02 Predicted regression lines -------------------------------------

# Displays the predicted regression lines from the model for each group
# (Species). This plot highlights trajectory of shape-size covariation 
# along the Csize variable. 

### Full dataset -----
plot(fit01_02, 
     reg.type = "PredLine", 
     predictor = forewing_gdf$Csize,
     col = fwSpOut_gdf$Species,
     type = "regression")

### Dataset without species outliers -----
plot(fit02_02, 
     reg.type = "PredLine", 
     predictor = fwSpOut_gdf$Csize,
     col = fwSpOut_gdf$Species,
     type = "regression")

### Dataset without sx outliers -----
plot(fit03_02, 
     reg.type = "PredLine", 
     predictor = fwSxOut_gdf$Csize,
     col = fwSxOut_gdf$Species,
     type = "regression")

### 04.05.01.03 Regression scores of individual specimens ----------------------

# Shows the regression scores of individual specimens along the trajectory
# predicted by the model. Similar to a scatter plot, useful for visualizing
# within-group variation relative to the regression line.

### Full dataset -----
plot(fit01_02, 
     reg.type = "RegScore", 
     predictor = log(forewing_gdf$Csize),
     col = forewing_gdf$Species,
     xlab = "Log of Centroid size",
     type = "regression")

### Dataset without species outliers -----
plot(fit02_02, 
     reg.type = "RegScore", 
     predictor = log(fwSpOut_gdf$Csize),
     col = fwSpOut_gdf$Species,
     xlab = "Log of Centroid size",
     type = "regression")

### Dataset without sex outliers -----
plot(fit03_02, 
     reg.type = "RegScore", 
     predictor = log(fwSxOut_gdf$Csize),
     col = fwSxOut_gdf$Species,
     xlab = "Log of Centroid size",
     type = "regression")

### 04.05.01.04 Checks model fit with diagnostic plots -------------------------

# Checks model fit with diagnostic plots similar to lm diagnostics:
# highlights residuals, leverage, outliers, and potential issues in
# the multivariate fit. Useful for assessing model reliability.

### Full dataset -----
plot(fit01_02, 
     reg.type = "diagnostics", 
     predictor = forewing_gdf$Csize,
     col = forewing_gdf$Species)

### Dataset without species outliers -----
plot(fit02_02, 
     reg.type = "diagnostics", 
     predictor = fwSpOut_gdf$Csize,
     col = fwSpOut_gdf$Species)

### Dataset without species outliers -----
plot(fit03_02, 
     reg.type = "diagnostics", 
     predictor = fwSxOut_gdf$Csize,
     col = fwSxOut_gdf$Species)