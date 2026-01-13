## ========================================================================== ##
## Script:       07_orthetrum_wings_forewing_shape_CVA_20260113.R             
## Author:       Matteo Zinni                                                 
## Date:         2026-01-13                                                   
## Description: Canonical Variate Analysis (CVA) and Linear Discriminant
##              Analysis (LDA) of Orthetrum forewing shape using geometric
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
fwFdCVA <- CVA(
  forewing_gpa$coords,
  forewing_gdf$Species,
  weighting = TRUE,
  tolinv = 1e-10,
  plot = FALSE,
  rounds = 9999,
  cv = TRUE,
  p.adjust.methods = "bonferroni"
)

# Display results
fwFdCVA

#### Extract features from CVA ----

# Percentage of variance explained by canonical axes
fwFdCVA$Var

# Canonical vectors
fwFdCVA$CV

# Canonical scores
fwFdCVA$CVscores

# Scores as data frame
fwFdCVA_df <- as.data.frame(fwFdCVA$CVscores)
fwFdCVA_df$Species <- rownames(fwFdCVA$CVscores)
rownames(fwFdCVA_df) <- NULL

# Plot result with ggplot
ggplot(fwFdCVA_df, aes(x = `CV 1`, y = `CV 2`, color = Species)) +
  geom_point(size = 4, alpha = 0.6) +
  labs(
    x = paste0("CV1 (", round(fwFdCVA$Var[1, 2], 2), "%)"),
    y = paste0("CV2 (", round(fwFdCVA$Var[2, 2], 2), "%)")
  )

#### Group distances ----------------------------------------------------------

## Mahalanobis distances quantify multivariate separation among group means
## accounting for within-group covariance, while Procrustes distances
## represent actual shape differences in morphospace.

# Mahalanobis distances among species means

# Get distances between groups
fwFdCVA$Dist$GroupdistMaha

# Get probabilities after permutation test
fwFdCVA$Dist$probsMaha

# Procrustes distances among species means

# Get distances between groups
fwFdCVA$Dist$GroupdistEuclid

# Get probabilites after permutation test
fwFdCVA$Dist$probsEuclid

## 07.02.02.01 Without species outliers ----------------------------------------

# Run the analysis
fwSpCVA <- CVA(
  fwSpOut_gpa$coords,
  fwSpOut_gpa$Species,
  weighting = TRUE,
  tolinv = 1e-10,
  plot = FALSE,
  rounds = 9999,
  cv = TRUE,
  p.adjust.methods = "bonferroni"
)

# Display results
fwSpCVA

#### Extract features from CVA ----

# Percentage of variance explained by canonical axes
fwSpCVA$Var

# Canonical vectors
fwSpCVA$CV

# Canonical scores
fwSpCVA$CVscores

# Scores as data frame
fwSpCVA_df <- as.data.frame(fwSpCVA$CVscores)
fwSpCVA_df$Species <- rownames(fwSpCVA$CVscores)
rownames(fwSpCVA_df) <- NULL

# Plot result with ggplot
ggplot(fwSpCVA_df, aes(x = `CV 1`, y = `CV 2`, color = Species)) +
  geom_point(size = 4, alpha = 0.6) +
  labs(
    x = paste0("CV1 (", round(fwSpCVA$Var[1, 2], 2), "%)"),
    y = paste0("CV2 (", round(fwSpCVA$Var[2, 2], 2), "%)")
  )

#### Group distances ----------------------------------------------------------

## Mahalanobis distances quantify multivariate separation among group means
## accounting for within-group covariance, while Procrustes distances
## represent actual shape differences in morphospace.

# Mahalanobis distances among species means

# Get distances between groups
fwSpCVA$Dist$GroupdistMaha

# Get probabilities after permutation test
fwSpCVA$Dist$probsMaha

# Procrustes distances among species means

# Get distances between groups
fwSpCVA$Dist$GroupdistEuclid

# Get probabilites after permutation test
fwSpCVA$Dist$probsEuclid

## 07.02.03.01 Without sex outliers --------------------------------------------

# Run the analysis
fwSxCVA <- CVA(
  fwSxOut_gpa$coords,
  fwSxOut_gpa$Species,
  weighting = TRUE,
  tolinv = 1e-10,
  plot = FALSE,
  rounds = 9999,
  cv = TRUE,
  p.adjust.methods = "bonferroni"
)

# Display results
fwSxCVA

#### Extract features from CVA ----

# Percentage of variance explained by canonical axes
fwSxCVA$Var

# Canonical vectors
fwSxCVA$CV

# Canonical scores
fwSxCVA$CVscores

# Scores as data frame
fwSxCVA_df <- as.data.frame(fwSxCVA$CVscores)
fwSxCVA_df$Sxecies <- rownames(fwSxCVA$CVscores)
rownames(fwSxCVA_df) <- NULL

# Plot result with ggplot
ggplot(fwSxCVA_df, aes(x = `CV 1`, y = `CV 2`, color = Species)) +
  geom_point(size = 4, alpha = 0.6) +
  labs(
    x = paste0("CV1 (", round(fwSxCVA$Var[1, 2], 2), "%)"),
    y = paste0("CV2 (", round(fwSxCVA$Var[2, 2], 2), "%)")
  )

#### Group distances ----------------------------------------------------------

## Mahalanobis distances quantify multivariate separation among group means
## accounting for within-group covariance, while Procrustes distances
## represent actual shape differences in morphospace.

# Mahalanobis distances among species means

# Get distances between groups
fwSxCVA$Dist$GroupdistMaha

# Get probabilities after permutation test
fwSxCVA$Dist$probsMaha

# Procrustes distances among species means

# Get distances between groups
fwSxCVA$Dist$GroupdistEuclid

# Get probabilites after permutation test
fwSxCVA$Dist$probsEuclid

## 07.03.01.01 LINEAR DISCRIMINANT ANALYSIS (LDA) ON SHAPE DATA ----------------

## When applied to two groups (e.g. Sex), CVA is mathematically equivalent
## to Linear Discriminant Analysis (LDA).

## 07.03.02.01 Full dataset ----------------------------------------------------

# Run the analysis
fwFdLDA <- CVA(
  forewing_gpa$coords,
  forewing_gdf$Sex,
  weighting = TRUE,
  tolinv = 1e-10,
  plot = FALSE,
  rounds = 9999,
  cv = TRUE,
  p.adjust.methods = "bonferroni"
)

# Display results
fwFdLDA

#### Extract features from LDA ----

# Percentage of variance explained by canonical axes
fwFdLDA$Var

# Canonical vectors
fwFdLDA$CV

# Canonical scores
fwFdLDA$CVscores

# Scores as data frame
fwFdLDA_df <- as.data.frame(fwFdLDA$CVscores)
fwFdLDA_df$Sex <- rownames(fwFdLDA$CVscores)
rownames(fwFdLDA_df) <- NULL

# Plot result with ggplot
ggplot(fwFdLDA_df, aes(x=`CV 1`, fill= Sex, color = Sex)) +
  geom_histogram(position="identity", alpha=0.65, bins=20)+
  ylab("Frequency") +
  xlab("1st Canonical Axis")+
  xlim(-6, 6)+
  scale_fill_manual(values=c("darkgoldenrod", "deepskyblue"))+
  scale_color_manual(values=c("black", "black"))

#### Group distances ----------------------------------------------------------

# Mahalanobis Dist. between group means

# Get distances between groups
fwFdLDA$Dist$GroupdistMaha

# Get probabilities after permutation test
fwFdLDA$Dist$probsMaha

# Procustes Dist. between group means 

# Get distances between groups
fwFdLDA$Dist$GroupdistEuclid

# Get probabilites after permutation test
fwFdLDA$Dist$probsEuclid

## 07.03.03.01 Without species outliers ----------------------------------------

# Run the analysis
fwSpLDA <- CVA(
  fwSpOut_gpa$coords,
  fwSpOut_gpa$Sex,
  weighting = TRUE,
  tolinv = 1e-10,
  plot = FALSE,
  rounds = 9999,
  cv = TRUE,
  p.adjust.methods = "bonferroni"
)

# Display results
fwSpLDA

#### Extract features from LDA ----

# Percentage of variance explained by canonical axes
fwSpLDA$Var

# Canonical vectors
fwSpLDA$CV

# Canonical scores
fwSpLDA$CVscores

# Scores as data frame
fwSpLDA_df <- as.data.frame(fwSpLDA$CVscores)
fwSpLDA_df$Sex <- rownames(fwSpLDA$CVscores)
rownames(fwSpLDA_df) <- NULL

# Plot result with ggplot
ggplot(fwSpLDA_df, aes(x=`CV 1`, fill= Sex, color = Sex)) +
  geom_histogram(position="identity", alpha=0.65, bins=20)+
  ylab("Frequency") +
  xlab("1st Canonical Axis")+
  xlim(-6, 6)+
  scale_fill_manual(values=c("darkgoldenrod", "deepskyblue"))+
  scale_color_manual(values=c("black", "black"))

#### Group distances ----------------------------------------------------------

# Mahalanobis Dist. between group means

# Get distances between groups
fwSpLDA$Dist$GroupdistMaha

# Get probabilities after permutation test
fwSpLDA$Dist$probsMaha

# Procustes Dist. between group means 

# Get distances between groups
fwSpLDA$Dist$GroupdistEuclid

# Get probabilites after permutation test
fwSpLDA$Dist$probsEuclid

## 07.03.03.01 Without sex outliers --------------------------------------------

# Run the analysis
fwSxLDA <- CVA(
  fwSxOut_gpa$coords,
  fwSxOut_gpa$Sex,
  weighting = TRUE,
  tolinv = 1e-10,
  plot = FALSE,
  rounds = 9999,
  cv = TRUE,
  p.adjust.methods = "bonferroni"
)

# Display results
fwSxLDA

#### Extract features from LDA ----

# Percentage of variance explained by canonical axes
fwSxLDA$Var

# Canonical vectors
fwSxLDA$CV

# Canonical scores
fwSxLDA$CVscores

# Scores as data frame
fwSxLDA_df <- as.data.frame(fwSxLDA$CVscores)
fwSxLDA_df$Sex <- rownames(fwSxLDA$CVscores)
rownames(fwSxLDA_df) <- NULL

# Plot result with ggplot
ggplot(fwSxLDA_df, aes(x=`CV 1`, fill= Sex, color = Sex)) +
  geom_histogram(position="identity", alpha=0.65, bins=20)+
  ylab("Frequency") +
  xlab("1st Canonical Axis")+
  xlim(-6, 6)+
  scale_fill_manual(values=c("darkgoldenrod", "deepskyblue"))+
  scale_color_manual(values=c("black", "black"))

#### Group distances ----------------------------------------------------------

# Mahalanobis Dist. between group means

# Get distances between groups
fwSxLDA$Dist$GroupdistMaha

# Get probabilities after permutation test
fwSxLDA$Dist$probsMaha

# Procustes Dist. between group means 

# Get distances between groups
fwSxLDA$Dist$GroupdistEuclid

# Get probabilites after permutation test
fwSxLDA$Dist$probsEuclid