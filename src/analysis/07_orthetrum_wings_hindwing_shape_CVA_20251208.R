## ========================================================================== ##
## Script:       07_orthetrum_wings_hindwing_shape_CVA_20250612.              ##
## Author:       Matteo Zinni                                                 ##
## Date:         2025-12-08                                                   ##
## Description:  A geometric morphometrics approach to investigate            ##
##               phylogenetic and taxonomic relationships among               ##
##               several skimmer species. Hindwing Shape Analysis - CVA       ##
## ========================================================================== ##

# 07.01.01.01 SHAPE ANALYSIS ---------------------------------------------------

# Packages loading ----
invisible(lapply(packages, function(pkg) {
  library(pkg, character.only = TRUE)
}))

## 07.02.01.01 CVA ON SHAPE DATA ------------------------------------------------

### Full Dataset ----

# Run the analysis
hwFdCVA <- CVA(hindwing_gpa$coords, specimen$Species, 
  weighting = TRUE, tolinv = 1e-10, 
  plot = TRUE, rounds = 9999, cv = TRUE, 
  p.adjust.methods = "bonferroni")

# Display results
hwFdCVA

#### Extract features from CVA ----

# Variance explained by the canonical roots:
hwFdCVA$Var

# Extract CVs
hwFdCVA$CV

# Extract scores
hwFdCVA$CVscores

# Scores as data frame
hwFdCVA_df = as.data.frame(hwFdCVA$CVscores)
hwFdCVA_df$Species = rownames(hwFdCVA$CVscores)
rownames(hwFdCVA_df) = NULL

# Plot result with ggplot
ggplot(hwFdCVA_df, aes(x= `CV 1`, y= `CV 2`)) +
  geom_point(aes(color=Species),  size = 4, alpha=0.5) +
  labs(y = paste("2nd Canonical Axis", 
  paste0("(", paste(round(hwFdCVA$Var[2,2],2),"%", sep = ""), ")")),
  x = paste("1st Canonical Axis", 
  paste0("(", paste(round(hwFdCVA$Var[1,2],2),"%", sep = ""), ")")))

# Mahalanobis Dist. between group means

# Get distances between groups
hwFdCVA$Dist$GroupdistMaha

# Get probabilities after permutation test
hwFdCVA$Dist$probsMaha

# Procustes Dist. between group means

# Get distances between groups
hwFdCVA$Dist$GroupdistEuclid

# Get probabilites after permutation test
hwFdCVA$Dist$probsEuclid

sink("session_info.txt")  # tutti i messaggi vengono scritti nel file
sessionInfo()
sink()    
### Without species outliers ----

# Run the analysis
hwSpCVA <- CVA(hwSpOut_gpa$coords, hwSpOut_gdf$Species, 
  weighting = TRUE, tolinv = 1e-10, 
  plot = TRUE, rounds = 9999, cv = TRUE, 
  p.adjust.methods = "bonferroni")

# Display results
hwSpCVA

#### Extract features from CVA ----

# Variance explained by the canonical roots:
hwSpCVA$Var

# Extract CVs
hwSpCVA$CV

# Extract scores
hwSpCVA$CVscores

# Scores as data frame
hwSpCVA_df = as.data.frame(hwSpCVA$CVscores)
hwSpCVA_df$Species = rownames(hwSpCVA$CVscores)
rownames(hwSpCVA_df) = NULL

# Plot result with ggplot
ggplot(hwSpCVA_df, aes(x= `CV 1`, y= `CV 2`)) +
  geom_point(aes(color=Species),  size = 4, alpha=0.5) +
  labs(y = paste("2nd Canonical Axis", 
  paste0("(", paste(round(hwSpCVA$Var[2,2],2),"%", sep = ""), ")")),
  x = paste("1st Canonical Axis", paste0("(", 
  paste(round(hwSpCVA$Var[1,2],2),"%", sep = ""), ")")))

# Mahalanobis Dist. between group means

# Get distances between groups
hwSpCVA$Dist$GroupdistMaha

# Get probabilities after permutation test
hwSpCVA$Dist$probsMaha

# Procustes Dist. between group means

# Get distances between groups
hwSpCVA$Dist$GroupdistEuclid

# Get probabilites after permutation test
hwSpCVA$Dist$probsEuclid

### Without species outliers ----

# Run the analysis
hwSxCVA <- CVA(hwSxOut_gpa$coords, hwSxOut_gdf$Species, 
  weighting = TRUE, tolinv = 1e-10, 
  plot = TRUE, rounds = 9999, cv = TRUE, 
  p.adjust.methods = "bonferroni")

# Display results
hwSxCVA

#### Extract features from CVA ----

# Variance explained by the canonical roots:
hwSxCVA$Var

# Extract CVs
hwSxCVA$CV

# Extract scores
hwSxCVA$CVscores

# Scores as data frame
hwSxCVA_df = as.data.frame(hwSxCVA$CVscores)
hwSxCVA_df$Species = rownames(hwSxCVA$CVscores)
rownames(hwSxCVA_df) = NULL

# Plot result with ggplot
ggplot(hwSxCVA_df, aes(x= `CV 1`, y= `CV 2`)) +
  geom_point(aes(color=Species),  size = 4, alpha=0.5) +
  labs(y = paste("2nd Canonical Axis", 
  paste0("(", paste(round(hwSxCVA$Var[2,2],2),"%", sep = ""), ")")),
  x = paste("1st Canonical Axis", paste0("(", 
  paste(round(hwSxCVA$Var[1,2],2),"%", sep = ""), ")")))

# Mahalanobis Dist. between group means

# Get distances between groups
hwSxCVA$Dist$GroupdistMaha

# Get probabilities after permutation test
hwSxCVA$Dist$probsMaha

# Procustes Dist. between group means

# Get distances between groups
hwSxCVA$Dist$GroupdistEuclid

# Get probabilites after permutation test
hwSxCVA$Dist$probsEuclid

## 07.02.02.01 LDA ON SHAPE DATA -----------------------------------------------

### Full Dataset ----

# Run the analysis
hwFdLDA <- CVA(hindwing_gpa$coords, specimen$Sex, 
  weighting = TRUE, tolinv = 1e-10, 
  plot = TRUE, rounds = 9999, cv = TRUE, 
  p.adjust.methods = "bonferroni")

# Display results
hwFdLDA

#### Extract features from CVA ----

# Variance explained by the canonical roots:
hwFdLDA$Var

# Extract CVs
hwFdLDA$CV

# Extract scores
hwSpLDA$CVscores

# Scores as data frame
hwFdLDA_df = as.data.frame(hwFdLDA$CVscores)
hwFdLDA_df$Sex = rownames(hwFdLDA$CVscores)
rownames(hwFdLDA_df) = NULL

# Plot result with ggplot
ggplot(hwFdLDA_df, aes(x=`CV 1`, fill= Sex, color = Sex)) +
  geom_histogram(position="identity", alpha=0.65, bins=20)+
  ylab("Frequency") +
  xlab("1st Canonical Axis")+
  xlim(-6, 6)+
  scale_fill_manual(values=c("darkgoldenrod", "deepskyblue"))+
  scale_color_manual(values=c("black", "black"))

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

### Without species outliers ----

# Run the analysis
hwSpLDA <- CVA(hwSpOut_gpa$coords, hwSpOut_gdf$Sex, 
  weighting = TRUE, tolinv = 1e-10, 
  plot = TRUE, rounds = 9999, cv = TRUE, 
  p.adjust.methods = "bonferroni")

# Display results
hwSpLDA

#### Extract features from CVA ----

# Variance explained by the canonical roots:
hwSpLDA$Var

# Extract CVs
hwSpLDA$CV

# Extract scores
hwSpLDA$CVscores

# Scores as data frame
hwSpLDA_df = as.data.frame(hwSpLDA$CVscores)
hwSpLDA_df$Sex = rownames(hwSpLDA$CVscores)
rownames(hwSpLDA_df) = NULL

# Plot result with ggplot
ggplot(hwSpLDA_df, aes(x=`CV 1`, fill= Sex, color = Sex)) +
  geom_histogram(position="identity", alpha=0.65, bins=20)+
  ylab("Frequency") +
  xlab("1st Canonical Axis")+
  xlim(-6, 6)+
  scale_fill_manual(values=c("darkgoldenrod", "deepskyblue"))+
  scale_color_manual(values=c("black", "black"))

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

### Without sex outliers ----

# Run the analysis
hwSxLDA <- CVA(hwSxOut_gpa$coords, hwSxOut_gdf$Sex, 
  weighting = TRUE, tolinv = 1e-10, 
  plot = TRUE, rounds = 9999, cv = TRUE, 
  p.adjust.methods = "bonferroni")

# Display results
hwSxLDA

#### Extract features from CVA ----

# Variance explained by the canonical roots:
hwSxLDA$Var

# Extract CVs
hwSxLDA$CV

# Extract scores
hwSxLDA$CVscores

# Scores as data frame
hwSxLDA_df = as.data.frame(hwSxLDA$CVscores)
hwSxLDA_df$Sex = rownames(hwSxLDA$CVscores)
rownames(hwSxLDA_df) = NULL

# Plot result with ggplot
ggplot(hwSxLDA_df, aes(x=`CV 1`, fill= Sex, color = Sex)) +
  geom_histogram(position="identity", alpha=0.65, bins=20)+
  ylab("Frequency") +
  xlab("1st Canonical Axis")+
  xlim(-6, 6)+
  scale_fill_manual(values=c("darkgoldenrod", "deepskyblue"))+
  scale_color_manual(values=c("black", "black"))

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

files <- list.files(all.files = TRUE, recursive = TRUE, full.names = TRUE)
sizes <- file.info(files)[, "size"] / 1024^2  # dimensioni in MB
cbind(files, sizes)


# Mostra solo i file >50 MB (per sicurezza)
big_files <- data.frame(file = files, size_MB = sizes)
big_files[big_files$size_MB > 50, ]

# Rimuove i file temporanei grandi
system("git rm --cached .RData")
system("git rm --cached .RDataTmp")