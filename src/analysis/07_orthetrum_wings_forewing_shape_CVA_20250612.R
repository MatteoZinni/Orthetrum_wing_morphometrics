## ========================================================================== ##
## Script:       07_orthetrum_wings_forewing_shape_CVA_20250612.              ##
## Author:       Matteo Zinni                                                 ##
## Date:         2025-06-12                                                   ##
## Description:  A geometric morphometrics approach to investigate            ##
##               phylogenetic and taxonomic relationships among               ##
##               several skimmer species. Forewing Shape Analysis - CVA       ##
## ========================================================================== ##

# 07.01.01.01 SHAPE ANALYSIS ---------------------------------------------------

# Packages loading ----
invisible(lapply(packages, function(pkg) {
  library(pkg, character.only = TRUE)
}))

## 07.02.01.01 CVA ON SHAPE DATA ------------------------------------------------

### Full Dataset ----

# Run the analysis
fwFdCVA <- CVA(forewing_gpa$coords, specimen$Species, 
  weighting = TRUE, tolinv = 1e-10, 
  plot = TRUE, rounds = 9999, cv = TRUE, 
  p.adjust.methods = "bonferroni")

# Display results
fwFdCVA

#### Extract features from CVA ----

# Variance explained by the canonical roots:
fwFdCVA$Var

# Extract CVs
fwFdCVA$CV

# Extract scores
fwFdCVA$CVscores

# Scores as data frame
fwFdCVA_df = as.data.frame(fwFdCVA$CVscores)
fwFdCVA_df$Species = rownames(fwFdCVA$CVscores)
rownames(fwFdCVA_df) = NULL

# Plot result with ggplot
ggplot(fwFdCVA_df, aes(x= `CV 1`, y= `CV 2`)) +
  geom_point(aes(color=Species),  size = 4, alpha=0.5) +
  labs(y = paste("2nd Canonical Axis", 
  paste0("(", paste(round(fwFdCVA$Var[2,2],2),"%", sep = ""), ")")),
  x = paste("1st Canonical Axis", 
  paste0("(", paste(round(fwFdCVA$Var[1,2],2),"%", sep = ""), ")")))

# Mahalanobis Dist. between group means

# Get distances between groups
fwFdCVA$Dist$GroupdistMaha

# Get probabilities after permutation test
fwFdCVA$Dist$probsMaha

# Procustes Dist. between group means

# Get distances between groups
fwFdCVA$Dist$GroupdistEuclid

# Get probabilites after permutation test
fwFdCVA$Dist$probsEuclid

### Without species outliers ----

# Run the analysis
fwSpCVA <- CVA(fwSpOut_gpa$coords, fwSpOut_gdf$Species, 
  weighting = TRUE, tolinv = 1e-10, 
  plot = TRUE, rounds = 9999, cv = TRUE, 
  p.adjust.methods = "bonferroni")

# Display results
fwSpCVA

#### Extract features from CVA ----

# Variance explained by the canonical roots:
fwSpCVA$Var

# Extract CVs
fwSpCVA$CV

# Extract scores
fwSpCVA$CVscores

# Scores as data frame
fwSpCVA_df = as.data.frame(fwSpCVA$CVscores)
fwSpCVA_df$Species = rownames(fwSpCVA$CVscores)
rownames(fwSpCVA_df) = NULL

# Plot result with ggplot
ggplot(fwSpCVA_df, aes(x= `CV 1`, y= `CV 2`)) +
  geom_point(aes(color=Species),  size = 4, alpha=0.5) +
  labs(y = paste("2nd Canonical Axis", 
  paste0("(", paste(round(fwSpCVA$Var[2,2],2),"%", sep = ""), ")")),
  x = paste("1st Canonical Axis", paste0("(", 
  paste(round(fwSpCVA$Var[1,2],2),"%", sep = ""), ")")))

# Mahalanobis Dist. between group means

# Get distances between groups
fwSpCVA$Dist$GroupdistMaha

# Get probabilities after permutation test
fwSpCVA$Dist$probsMaha

# Procustes Dist. between group means

# Get distances between groups
fwSpCVA$Dist$GroupdistEuclid

# Get probabilites after permutation test
fwSpCVA$Dist$probsEuclid

### Without species outliers ----

# Run the analysis
fwSxCVA <- CVA(fwSxOut_gpa$coords, fwSxOut_gdf$Species, 
  weighting = TRUE, tolinv = 1e-10, 
  plot = TRUE, rounds = 9999, cv = TRUE, 
  p.adjust.methods = "bonferroni")

# Display results
fwSxCVA

#### Extract features from CVA ----

# Variance explained by the canonical roots:
fwSxCVA$Var

# Extract CVs
fwSxCVA$CV

# Extract scores
fwSxCVA$CVscores

# Scores as data frame
fwSxCVA_df = as.data.frame(fwSxCVA$CVscores)
fwSxCVA_df$Species = rownames(fwSxCVA$CVscores)
rownames(fwSxCVA_df) = NULL

# Plot result with ggplot
ggplot(fwSxCVA_df, aes(x= `CV 1`, y= `CV 2`)) +
  geom_point(aes(color=Species),  size = 4, alpha=0.5) +
  labs(y = paste("2nd Canonical Axis", 
  paste0("(", paste(round(fwSxCVA$Var[2,2],2),"%", sep = ""), ")")),
  x = paste("1st Canonical Axis", paste0("(", 
  paste(round(fwSxCVA$Var[1,2],2),"%", sep = ""), ")")))

# Mahalanobis Dist. between group means

# Get distances between groups
fwSxCVA$Dist$GroupdistMaha

# Get probabilities after permutation test
fwSxCVA$Dist$probsMaha

# Procustes Dist. between group means

# Get distances between groups
fwSxCVA$Dist$GroupdistEuclid

# Get probabilites after permutation test
fwSxCVA$Dist$probsEuclid

## 07.02.02.01 LDA ON SHAPE DATA ------------------------------------------------

### Full Dataset ----

# Run the analysis
fwFdLDA <- CVA(forewing_gpa$coords, specimen$Sex, 
  weighting = TRUE, tolinv = 1e-10, 
  plot = TRUE, rounds = 9999, cv = TRUE, 
  p.adjust.methods = "bonferroni")

# Display results
fwFdLDA

#### Extract features from CVA ----

# Variance explained by the canonical roots:
fwFdLDA$Var

# Extract CVs
fwFdLDA$CV

# Extract scores
fwSpLDA$CVscores

# Scores as data frame
fwFdLDA_df = as.data.frame(fwFdLDA$CVscores)
fwFdLDA_df$Sex = rownames(fwFdLDA$CVscores)
rownames(fwFdLDA_df) = NULL

# Plot result with ggplot
ggplot(fwFdLDA_df, aes(x=`CV 1`, fill= Sex, color = Sex)) +
  geom_histogram(position="identity", alpha=0.65, bins=20)+
  ylab("Frequency") +
  xlab("1st Canonical Axis")+
  xlim(-6, 6)+
  scale_fill_manual(values=c("darkgoldenrod", "deepskyblue"))+
  scale_color_manual(values=c("black", "black"))

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

### Without species outliers ----

# Run the analysis
fwSpLDA <- CVA(fwSpOut_gpa$coords, fwSpOut_gdf$Sex, 
  weighting = TRUE, tolinv = 1e-10, 
  plot = TRUE, rounds = 9999, cv = TRUE, 
  p.adjust.methods = "bonferroni")

# Display results
fwSpLDA

#### Extract features from CVA ----

# Variance explained by the canonical roots:
fwSpLDA$Var

# Extract CVs
fwSpLDA$CV

# Extract scores
fwSpLDA$CVscores

# Scores as data frame
fwSpLDA_df = as.data.frame(fwSpLDA$CVscores)
fwSpLDA_df$Sex = rownames(fwSpLDA$CVscores)
rownames(fwSpLDA_df) = NULL

# Plot result with ggplot
ggplot(fwSpLDA_df, aes(x=`CV 1`, fill= Sex, color = Sex)) +
  geom_histogram(position="identity", alpha=0.65, bins=20)+
  ylab("Frequency") +
  xlab("1st Canonical Axis")+
  xlim(-6, 6)+
  scale_fill_manual(values=c("darkgoldenrod", "deepskyblue"))+
  scale_color_manual(values=c("black", "black"))

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

### Without sex outliers ----

# Run the analysis
fwSxLDA <- CVA(fwSxOut_gpa$coords, fwSxOut_gdf$Sex, 
               weighting = TRUE, tolinv = 1e-10, 
               plot = TRUE, rounds = 9999, cv = TRUE, 
               p.adjust.methods = "bonferroni")

# Display results
fwSxLDA

#### Extract features from CVA ----

# Variance explained by the canonical roots:
fwSxLDA$Var

# Extract CVs
fwSxLDA$CV

# Extract scores
fwSxLDA$CVscores

# Scores as data frame
fwSxLDA_df = as.data.frame(fwSxLDA$CVscores)
fwSxLDA_df$Sex = rownames(fwSxLDA$CVscores)
rownames(fwSxLDA_df) = NULL

# Plot result with ggplot
ggplot(fwSxLDA_df, aes(x=`CV 1`, fill= Sex, color = Sex)) +
  geom_histogram(position="identity", alpha=0.65, bins=20)+
  ylab("Frequency") +
  xlab("1st Canonical Axis")+
  xlim(-6, 6)+
  scale_fill_manual(values=c("darkgoldenrod", "deepskyblue"))+
  scale_color_manual(values=c("black", "black"))

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