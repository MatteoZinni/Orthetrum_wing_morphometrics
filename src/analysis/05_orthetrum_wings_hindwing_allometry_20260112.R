## ========================================================================== ##
## Script:      05_orthetrum_wings_hindwing_allometry_20260112.R
## Author:      Matteo Zinni
## Date:        2026-01-12
## Description: Allometric analysis of Orthetrum hindwing shape using
##              geometric morphometrics. The script investigates centroid
##              size variation and size–shape relationships across species
##              and sexes, using full datasets and outlier-filtered subsets.
## ========================================================================== ##

# 05.01.01.01 HINDWING ANALYSIS ------------------------------------------------

# Packages loading ----
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}

# Sourcing the function
source(here::here("src", "functions", "packages_setup_20260111.R"))

# Lunch the function
install_packages()
## 05.01.01.01 CENTROID SIZE ---------------------------------------------------

## 05.01.02.01 Full dataset ----------------------------------------------------

# Distribution of centroid size values
hist(hindwing_2d$C_size,
     main = "Forewing centroid size distribution (hindwing)",
     xlab = "Centroid size")

#### Does centroid size differ among species? ----

# Testing for difference (chi-squared = 34.641, df = 9, p-value = 6.893e-05)
kruskal.test(C_size ~ Species, data = hindwing_2d)

# Post-hoc test for pairwise comparison
kwAllPairsDunnTest(C_size ~ Species, 
                   data = hindwing_2d,
                   p.adjust.method = "bonferroni"
                   )

# Plot results
boxplot(C_size ~ Species, data = hindwing_2d,
        ylab = "Centroid size",
        main = "Forewing centroid size among species \n (full dataset)")

#### Does centroid size differ between sexes? ----

# Testing for difference (chi-squared = 0.0269, df = 1, p-value = 0.8696)
kruskal.test(C_size ~ Sex, data = hindwing_2d)

# Post-hoc test for pairwise comparison
kwAllPairsDunnTest(C_size ~ Sex, 
                   data = hindwing_2d,
                   p.adjust.method = "bonferroni"
                   )

# Plot results
boxplot(C_size ~ Sex, data = hindwing_2d,
        ylab = "Centroid size",
        main = "Forewing centroid size between sexes \n (full dataset)")

# Distribution of centroid size values between sexes 
hist(hindwing_2d$C_size[hindwing_2d$Sex == 'Male'], 
     col='deepskyblue',
     main='Forewing centroid size between sexes \n Full dataset', 
     xlab='Centroid size')
hist(hindwing_2d$C_size[hindwing_2d$Sex == 'Female'], 
     col='darkgoldenrod',
     add=TRUE)

## 05.01.03.01 Dataset without species outliers --------------------------------

# Distribution of centroid size values
hist(hwSpOut_2d$C_size)

#### Does centroid size differ among species? ----

# Testing for difference (chi-squared = 33.381, df = 9, p-value = 0.0001146)
kruskal.test(C_size ~ Species, data = hwSpOut_2d)

# Post-hoc test for pairwise comparison
kwAllPairsDunnTest(C_size ~ Species, data = hwSpOut_2d,
                   p.adjust.method = "bonferroni")

# Plot results
boxplot(C_size ~ Species, data = hwSpOut_2d,
        ylab = "Centroid size",
        main = "Centroid size among species")

## 05.01.04.01 Dataset without sex outliers ------------------------------------

# Distribution of centroid size values
hist(hwSxOut_2d$C_size)

#### Does centroid size differ between sexes? ----

# Testing for difference (chi-squared = 0.10185, df = 1, p-value = 0.7496)
kruskal.test(C_size ~ Sex, data = hwSxOut_2d)

# Plot results
boxplot(C_size ~ Sex, data = hwSxOut_2d,
        ylab = "Centroid size",
        main = "Centroid size between sexes")

# Distribution of centroid size values between sexes 
hist(hwSxOut_2d$C_size[hwSxOut_2d$Sex == 'Male'], 
     col='deepskyblue',
     main='Centroid size between sexes', xlab='Centroid size')
hist(hwSxOut_2d$C_size[hwSxOut_2d$Sex == 'Female'], 
     col='darkgoldenrod',
     add=TRUE)

## 05.02.01.01 ALLOMETRY ANALYSIS ----------------------------------------------

## 05.02.02.01 Full dataset ----------------------------------------------------

#### Effect of size on shape ----

# Null model  
fit04_00 <-  procD.lm(coords ~ 1, 
                      RRPP=T, 
                      data = hindwing_gdf, 
                      iter = 9999,
                      seed = 123)

# model summary 
summary(fit04_00)

# Simple size on shape covariation model
fit04_01 <- procD.lm(coords~log(Csize),
                     data = hindwing_gdf, 
                     RRPP=T,
                     iter = 9999, 
                     seed = 123)

# model summary 
summary(fit04_01)

#### Species-specific allometric trajectories (unique slopes) ----

# Species allometry | full dataset | unique slopes
fit04_02 <- procD.lm(coords~log(Csize)*Species,
                     data = hindwing_gdf,
                     SS.type = "I",
                     effect.type = "F",
                     RRPP = TRUE,
                     iter = 9999, 
                     seed = 123)

# model summary 
summary(fit04_02)

#### Common allometric slope with species-specific shape intercepts ----

# Species allometry | full dataset | common slopes
fit04_03 <- procD.lm(coords~log(Csize) + Species,
                     data = hindwing_gdf,
                     SS.type = "I",
                     effect.type = "F",
                     RRPP = TRUE,
                     iter = 9999, 
                     seed = 123)

# model summary 
summary(fit04_03)

#### Model comparison and selection ----

# Comparing models with anova: fit04_02 (unique) is the best candidate
anova(fit04_00, fit04_01, fit04_02, fit04_03)

## 05.02.02.02 Pairwise comparison of allometric trajectories ----

# Compute pairwise comparisons of unique allometric slopes among species
hwFd_pw <- pairwise(fit04_02, 
                       groups = hindwing_gdf$Species, 
                       covariate = log(hindwing_gdf$Csize))

# Differences in slope vector length
summary(hwFd_pw, confidence = 0.95, test.type = "dist", stat.table = FALSE)

# Differences in slope vector correlation (VC = vector correlation)
summary(hwFd_pw, confidence = 0.95, test.type = "VC", stat.table = FALSE)

# Differences in directional vector (DL = direction length)
summary(hwFd_pw, confidence = 0.95, test.type = "DL", stat.table = FALSE)

## 05.02.03.01 Dataset without species outliers --------------------------------

#### Effect of size on shape ----

# Null model  
fit05_00 <-  procD.lm(coords ~ 1, 
                      RRPP=T, 
                      data = hwSpOut_gdf, 
                      iter = 9999,
                      seed = 123)

# model summary 
summary(fit05_00)

# Simple size on shape covariation model
fit05_01 <- procD.lm(coords~log(Csize),
                     data = hwSpOut_gdf, 
                     RRPP=T,
                     iter = 9999, 
                     seed = 123)

# model summary 
summary(fit05_01)

#### Species-specific allometric trajectories (unique slopes) ----

# Species allometry | no species outliers | unique slopes
fit05_02 <- procD.lm(coords~log(Csize)*Species,
                     data = hwSpOut_gdf,
                     SS.type = "I",
                     effect.type = "F",
                     RRPP = TRUE,
                     iter = 9999, 
                     seed = 123)

# model summary 
summary(fit05_02)

#### Common allometric slope with species-specific shape intercepts ----

# Species allometry | no species outliers | common slopes
fit05_03 <- procD.lm(coords~log(Csize)+ Species,
                     data = hwSpOut_gdf,
                     SS.type = "I",
                     effect.type = "F",
                     RRPP = TRUE,
                     iter = 9999, 
                     seed = 123)

# model summary 
summary(fit05_03)

#### Model comparison and selection ----

# Comparing models with anova: fit05_02 (unique) is the best candidate
anova(fit05_00, fit05_01, fit05_02, fit05_03)

# Pairwise comparison
hwSpOut_pw = pairwise(fit05_02, 
                      groups = hwSpOut_gdf$Species, 
                      covariate = log(hwSpOut_gdf$Csize))

## 05.02.03.02 Pairwise comparison of allometric trajectories ----

# Compute pairwise comparisons of unique allometric slopes among species
hwSpOut_pw <- pairwise(fit05_02, 
                       groups = hwSpOut_gdf$Species, 
                       covariate = log(hwSpOut_gdf$Csize))

# Differences in slope vector length
summary(hwSpOut_pw, confidence = 0.95, test.type = "dist", stat.table = FALSE)

# Differences in slope vector correlation (VC = vector correlation)
summary(hwSpOut_pw, confidence = 0.95, test.type = "VC", stat.table = FALSE)

# Differences in directional vector (DL = direction length)
summary(hwSpOut_pw, confidence = 0.95, test.type = "DL", stat.table = FALSE)

## 05.02.04.01 Dataset without sex outliers ------------------------------------

#### Effect of size on shape ----

# Null model  
fit06_00 <-  procD.lm(coords ~ 1, 
                      RRPP = TRUE, 
                      data = hwSxOut_gdf, 
                      iter = 9999,
                      seed = 123)

# Model summary 
summary(fit06_00)

#### Species-specific allometric trajectories (unique slopes) ----

# Species allometry | no sex outliers | unique slopes
fit06_01 <- procD.lm(coords ~ log(Csize),
                     data = hwSxOut_gdf, 
                     RRPP = TRUE,
                     iter = 9999, 
                     seed = 123)

# Model summary 
summary(fit06_01)

#### Common allometric slope with species-specific shape intercepts ----

# Species allometry | no sex outliers | common slopes
fit06_02 <- procD.lm(coords ~ log(Csize)*Sex,
                     data = hwSxOut_gdf,
                     SS.type = "I",
                     effect.type = "F",
                     RRPP = TRUE,
                     iter = 9999, 
                     seed = 123)

# Model summary 
summary(fit06_02)

# Common slope model 
fit06_03 <- procD.lm(coords ~ log(Csize) + Sex,
                     data = hwSxOut_gdf,
                     SS.type = "I",
                     effect.type = "F",
                     RRPP = TRUE,
                     iter = 9999, 
                     seed = 123)

# Model summary 
summary(fit06_03)

#### Model comparison and selection ----

# Comparing models with anova: fit06_02 (unique) is still the best candidate
anova(fit06_01, fit06_02, fit06_03)

## 05.02.04.02 Pairwise comparison of allometric trajectories ----

# Compute pairwise comparisons of unique allometric slopes among species
hwSxOut_pw <- pairwise(fit06_02, 
                       groups = hwSxOut_gdf$Species, 
                       covariate = log(hwSxOut_gdf$Csize))

# Differences in slope vector length
summary(hwSxOut_pw, confidence = 0.95, test.type = "dist", stat.table = FALSE)

# Differences in slope vector correlation (VC = vector correlation)
summary(hwSxOut_pw, confidence = 0.95, test.type = "VC", stat.table = FALSE)

# Differences in directional vector (DL = direction length)
summary(hwSxOut_pw, confidence = 0.95, test.type = "DL", stat.table = FALSE)

## 04.05.01.01 PLOT AND DATA VISUALIZATION -------------------------------------

### 04.05.01.02 Predicted regression lines -------------------------------------

# Displays the predicted regression lines from the model for each group
# (Species). This plot highlights trajectory of shape-size covariation 
# along the Csize variable. 

### Full dataset -----
plot(fit04_02, 
     reg.type = "PredLine", 
     predictor = hindwing_gdf$Csize,
     col = hwSpOut_gdf$Species,
     type = "regression")

### Dataset without species outliers -----
plot(fit05_02, 
     reg.type = "PredLine", 
     predictor = hwSpOut_gdf$Csize,
     col = hwSpOut_gdf$Species,
     type = "regression")

### Dataset without sx outliers -----
plot(fit06_02, 
     reg.type = "PredLine", 
     predictor = hwSxOut_gdf$Csize,
     col = hwSxOut_gdf$Species,
     type = "regression")

### 04.05.01.03 Regression scores of individual specimens ----------------------

# Shows the regression scores of individual specimens along the trajectory
# predicted by the model. Similar to a scatter plot, useful for visualizing
# within-group variation relative to the regression line.

### Full dataset -----
plot(fit04_02, 
     reg.type = "RegScore", 
     predictor = log(hindwing_gdf$Csize),
     col = hindwing_gdf$Species,
     xlab = "Log of Centroid size",
     type = "regression")

### Dataset without species outliers -----
plot(fit05_02, 
     reg.type = "RegScore", 
     predictor = log(hwSpOut_gdf$Csize),
     col = hwSpOut_gdf$Species,
     xlab = "Log of Centroid size",
     type = "regression")

### Dataset without sex outliers -----
plot(fit06_02, 
     reg.type = "RegScore", 
     predictor = log(hwSxOut_gdf$Csize),
     col = hwSxOut_gdf$Species,
     xlab = "Log of Centroid size",
     type = "regression")