## ========================================================================== ##
## Script:       07_orthetrum_wings_hindwing_shape_CVA_20260113.R             
## Author:       Matteo Zinni                                                 
## Date:         2026-01-13                                                   
## Description: Canonical Variate Analysis (CVA) and Linear Discriminant
##              Analysis (LDA) of Orthetrum hindwing shape using geometric
##              morphometrics. The script investigates species- and sex-level
##              shape discrimination using full datasets and outlier-filtered
##              subsets.
## ========================================================================== ##

# 07.01.01.01 SHAPE ANALYSIS ---------------------------------------------------

# Packages loading ----
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}

# Sourcing the function
source(here::here("src", "functions", "packages_setup_20260111.R"))

# Lunch the function
install_packages()

## 07.02.01.01 CVA ON SHAPE DATA -----------------------------------------------

## Canonical Variate Analysis (CVA) is a supervised method that maximizes
## shape differences among predefined groups. Therefore, results should
## be interpreted as group discrimination rather than unsupervised
## patterns of shape variation.

## 07.02.01.02 Full dataset ----------------------------------------------------

# Run the analysis
hwFdCVA <- CVA(
  hindwing_gpa$coords,
  hindwing_gdf$Species,
  weighting = TRUE,
  tolinv = 1e-10,
  plot = FALSE,
  rounds = 9999,
  cv = TRUE,
  p.adjust.methods = "bonferroni"
)

# Display results
hwFdCVA

#### Extract features from CVA ----

# Percentage of variance explained by canonical axes
hwFdCVA$Var

# Canonical vectors
hwFdCVA$CV

# Canonical scores
hwFdCVA$CVscores

# Scores as data frame
hwFdCVA_df <- as.data.frame(hwFdCVA$CVscores)
hwFdCVA_df$Species <- rownames(hwFdCVA$CVscores)
rownames(hwFdCVA_df) <- NULL

# Plot result with ggplot
ggplot(hwFdCVA_df, aes(x = `CV 1`, y = `CV 2`, color = Species)) +
  geom_point(size = 4, alpha = 0.6) +
  labs(
    x = paste0("CV1 (", round(hwFdCVA$Var[1, 2], 2), "%)"),
    y = paste0("CV2 (", round(hwFdCVA$Var[2, 2], 2), "%)")
  )

#### Group distances ----------------------------------------------------------

## Mahalanobis distances quantify multivariate separation among group means
## accounting for within-group covariance, while Procrustes distances
## represent actual shape differences in morphospace.

# Mahalanobis distances among species means

# Get distances between groups
hwFdCVA$Dist$GroupdistMaha

# Get probabilities after permutation test
hwFdCVA$Dist$probsMaha

# Procrustes distances among species means

# Get distances between groups
hwFdCVA$Dist$GroupdistEuclid

# Get probabilites after permutation test
hwFdCVA$Dist$probsEuclid

## 07.02.02.01 Without species outliers ----------------------------------------

# Run the analysis
hwSpCVA <- CVA(
  hwSpOut_gpa$coords,
  hwSpOut_gpa$Species,
  weighting = TRUE,
  tolinv = 1e-10,
  plot = FALSE,
  rounds = 9999,
  cv = TRUE,
  p.adjust.methods = "bonferroni"
)

# Display results
hwSpCVA

#### Extract features from CVA ----

# Percentage of variance explained by canonical axes
hwSpCVA$Var

# Canonical vectors
hwSpCVA$CV

# Canonical scores
hwSpCVA$CVscores

# Scores as data frame
hwSpCVA_df <- as.data.frame(hwSpCVA$CVscores)
hwSpCVA_df$Species <- rownames(hwSpCVA$CVscores)
rownames(hwSpCVA_df) <- NULL

# Plot result with ggplot
ggplot(hwSpCVA_df, aes(x = `CV 1`, y = `CV 2`, color = Species)) +
  geom_point(size = 4, alpha = 0.6) +
  labs(
    x = paste0("CV1 (", round(hwSpCVA$Var[1, 2], 2), "%)"),
    y = paste0("CV2 (", round(hwSpCVA$Var[2, 2], 2), "%)")
  )

#### Group distances ----------------------------------------------------------

## Mahalanobis distances quantify multivariate separation among group means
## accounting for within-group covariance, while Procrustes distances
## represent actual shape differences in morphospace.

# Mahalanobis distances among species means

# Get distances between groups
hwSpCVA$Dist$GroupdistMaha

# Get probabilities after permutation test
hwSpCVA$Dist$probsMaha

# Procrustes distances among species means

# Get distances between groups
hwSpCVA$Dist$GroupdistEuclid

# Get probabilites after permutation test
hwSpCVA$Dist$probsEuclid

## 07.02.03.01 Without sex outliers --------------------------------------------

# Run the analysis
hwSxCVA <- CVA(
  hwSxOut_gpa$coords,
  hwSxOut_gpa$Species,
  weighting = TRUE,
  tolinv = 1e-10,
  plot = FALSE,
  rounds = 9999,
  cv = TRUE,
  p.adjust.methods = "bonferroni"
)

# Display results
hwSxCVA

#### Extract features from CVA ----

# Percentage of variance explained by canonical axes
hwSxCVA$Var

# Canonical vectors
hwSxCVA$CV

# Canonical scores
hwSxCVA$CVscores

# Scores as data frame
hwSxCVA_df <- as.data.frame(hwSxCVA$CVscores)
hwSxCVA_df$Sxecies <- rownames(hwSxCVA$CVscores)
rownames(hwSxCVA_df) <- NULL

# Plot result with ggplot
ggplot(hwSxCVA_df, aes(x = `CV 1`, y = `CV 2`, color = Species)) +
  geom_point(size = 4, alpha = 0.6) +
  labs(
    x = paste0("CV1 (", round(hwSxCVA$Var[1, 2], 2), "%)"),
    y = paste0("CV2 (", round(hwSxCVA$Var[2, 2], 2), "%)")
  )

#### Group distances ----------------------------------------------------------

## Mahalanobis distances quantify multivariate separation among group means
## accounting for within-group covariance, while Procrustes distances
## represent actual shape differences in morphospace.

# Mahalanobis distances among species means

# Get distances between groups
hwSxCVA$Dist$GroupdistMaha

# Get probabilities after permutation test
hwSxCVA$Dist$probsMaha

# Procrustes distances among species means

# Get distances between groups
hwSxCVA$Dist$GroupdistEuclid

# Get probabilites after permutation test
hwSxCVA$Dist$probsEuclid

## 07.03.01.01 LINEAR DISCRIMINANT ANALYSIS (LDA) ON SHAPE DATA ----------------

## When applied to two groups (e.g. Sex), CVA is mathematically equivalent
## to Linear Discriminant Analysis (LDA).

## 07.03.02.01 Full dataset ----------------------------------------------------

# Run the analysis
hwFdLDA <- CVA(
  hindwing_gpa$coords,
  hindwing_gdf$Sex,
  weighting = TRUE,
  tolinv = 1e-10,
  plot = FALSE,
  rounds = 9999,
  cv = TRUE,
  p.adjust.methods = "bonferroni"
)

# Display results
hwFdLDA

#### Extract features from LDA ----

# Percentage of variance explained by canonical axes
hwFdLDA$Var

# Canonical vectors
hwFdLDA$CV

# Canonical scores
hwFdLDA$CVscores

# Scores as data frame
hwFdLDA_df <- as.data.frame(hwFdLDA$CVscores)
hwFdLDA_df$Sex <- rownames(hwFdLDA$CVscores)
rownames(hwFdLDA_df) <- NULL

# Plot result with ggplot
ggplot(hwFdLDA_df, aes(x=`CV 1`, fill= Sex, color = Sex)) +
  geom_histogram(position="identity", alpha=0.65, bins=20)+
  ylab("Frequency") +
  xlab("1st Canonical Axis")+
  xlim(-6, 6)+
  scale_fill_manual(values=c("darkgoldenrod", "deepskyblue"))+
  scale_color_manual(values=c("black", "black"))

#### Group distances ----------------------------------------------------------

# Mahalanobis Dist. between group means

# Get distances between groups
hwFdLDA$Dist$GroupdistMaha

# Get probabilities after permutation test
hwFdLDA$Dist$probsMaha

# Procustes Dist. between group means 

# Get distances between groups
hwFdLDA$Dist$GroupdistEuclid

# Get probabilites after permutation test
hwFdLDA$Dist$probsEuclid

## 07.03.03.01 Without species outliers ----------------------------------------

# Run the analysis
hwSpLDA <- CVA(
  hwSpOut_gpa$coords,
  hwSpOut_gpa$Sex,
  weighting = TRUE,
  tolinv = 1e-10,
  plot = FALSE,
  rounds = 9999,
  cv = TRUE,
  p.adjust.methods = "bonferroni"
)

# Display results
hwSpLDA

#### Extract features from LDA ----

# Percentage of variance explained by canonical axes
hwSpLDA$Var

# Canonical vectors
hwSpLDA$CV

# Canonical scores
hwSpLDA$CVscores

# Scores as data frame
hwSpLDA_df <- as.data.frame(hwSpLDA$CVscores)
hwSpLDA_df$Sex <- rownames(hwSpLDA$CVscores)
rownames(hwSpLDA_df) <- NULL

# Plot result with ggplot
ggplot(hwSpLDA_df, aes(x=`CV 1`, fill= Sex, color = Sex)) +
  geom_histogram(position="identity", alpha=0.65, bins=20)+
  ylab("Frequency") +
  xlab("1st Canonical Axis")+
  xlim(-6, 6)+
  scale_fill_manual(values=c("darkgoldenrod", "deepskyblue"))+
  scale_color_manual(values=c("black", "black"))

#### Group distances ----------------------------------------------------------

# Mahalanobis Dist. between group means

# Get distances between groups
hwSpLDA$Dist$GroupdistMaha

# Get probabilities after permutation test
hwSpLDA$Dist$probsMaha

# Procustes Dist. between group means 

# Get distances between groups
hwSpLDA$Dist$GroupdistEuclid

# Get probabilites after permutation test
hwSpLDA$Dist$probsEuclid

## 07.03.03.01 Without sex outliers --------------------------------------------

# Run the analysis
hwSxLDA <- CVA(
  hwSxOut_gpa$coords,
  hwSxOut_gpa$Sex,
  weighting = TRUE,
  tolinv = 1e-10,
  plot = FALSE,
  rounds = 9999,
  cv = TRUE,
  p.adjust.methods = "bonferroni"
)

# Display results
hwSxLDA

#### Extract features from LDA ----

# Percentage of variance explained by canonical axes
hwSxLDA$Var

# Canonical vectors
hwSxLDA$CV

# Canonical scores
hwSxLDA$CVscores

# Scores as data frame
hwSxLDA_df <- as.data.frame(hwSxLDA$CVscores)
hwSxLDA_df$Sex <- rownames(hwSxLDA$CVscores)
rownames(hwSxLDA_df) <- NULL

# Plot result with ggplot
ggplot(hwSxLDA_df, aes(x=`CV 1`, fill= Sex, color = Sex)) +
  geom_histogram(position="identity", alpha=0.65, bins=20)+
  ylab("Frequency") +
  xlab("1st Canonical Axis")+
  xlim(-6, 6)+
  scale_fill_manual(values=c("darkgoldenrod", "deepskyblue"))+
  scale_color_manual(values=c("black", "black"))

#### Group distances ----------------------------------------------------------

# Mahalanobis Dist. between group means

# Get distances between groups
hwSxLDA$Dist$GroupdistMaha

# Get probabilities after permutation test
hwSxLDA$Dist$probsMaha

# Procustes Dist. between group means 

# Get distances between groups
hwSxLDA$Dist$GroupdistEuclid

# Get probabilites after permutation test
hwSxLDA$Dist$probsEuclid