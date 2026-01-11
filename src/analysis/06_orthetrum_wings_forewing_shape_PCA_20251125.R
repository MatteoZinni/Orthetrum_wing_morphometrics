## ========================================================================== ##
## Script:       06_orthetrum_wings_forewing_shape_PCA_20251126.              ##
## Author:       Matteo Zinni                                                 ##
## Date:         2025-04-24                                                   ##
## Description:  A geometric morphometrics approach to investigate            ##
##               phylogenetic and taxonomic relationships among               ##
##               several skimmer species. Forewing Shape Analysis - PCA.      ##
## ========================================================================== ##

# 06.01.01.01 SHAPE ANALYSIS ---------------------------------------------------

# Packages loading ----
invisible(lapply(packages, function(pkg) {
  library(pkg, character.only = TRUE)
}))

## 06.02.01.01 PCA ON SHAPE DATA -----------------------------------------------

### PCA on the full dataset ----
fwFdPCA = gm.prcomp(forewing_gpa$coords, transform = T)

# summary
summary(fwFdPCA)

#### Extract features from PCA ----

# Scores
fwFdPCA_scores = as.data.frame(fwFdPCA$x)
fwFdPCA_scores$Sex = forewing_2d$Sex
fwFdPCA_scores$Species = forewing_2d$Species

# Test for sex separation along PC1 
kruskal.test(fwFdPCA_scores$Comp1 ~ fwFdPCA_scores$Sex)

# Test for sex separation along PC2
kruskal.test(fwFdPCA_scores$Comp2 ~ fwFdPCA_scores$Sex)

# Test for sex separation along PC3
kruskal.test(fwFdPCA_scores$Comp3 ~ fwFdPCA_scores$Sex)

# Test for species separation along PC1 
kruskal.test(fwFdPCA_scores$Comp1 ~ fwFdPCA_scores$Species)

# Test for species separation along PC2
kruskal.test(fwFdPCA_scores$Comp2 ~ fwFdPCA_scores$Species)

# Test for species separation along PC3
kruskal.test(fwFdPCA_scores$Comp3 ~ fwFdPCA_scores$Species)

# Post-hoc test for pairwise comparison
kwAllPairsDunnTest(Comp1 ~ Species, data = fwFdPCA_scores,
                   p.adjust.method = "bonferroni")

# Post-hoc test for pairwise comparison
kwAllPairsDunnTest(Comp2 ~ Species, data = fwFdPCA_scores,
                   p.adjust.method = "bonferroni")

# Post-hoc test for pairwise comparison
kwAllPairsDunnTest(Comp3 ~ Species, data = fwFdPCA_scores,
                   p.adjust.method = "bonferroni")

plotRefToTarget(mshape(forewing_gpa$coords[,,specimen$Sex=="Male"]),
                mshape(forewing_gpa$coords[,,specimen$Sex=="Female"]),
                mag = 1,
                links = fw_links)


# Scoreplot with ggplot2: Species
ggplot2::ggplot(fwFdPCA_scores, aes(x=Comp1, y=Comp2, color = Species)) + 
  geom_point(shape = 16, alpha=0.50, size = 5)+
  scale_x_continuous(limits = c(-0.1, 0.1)) + 
  xlab("PC1 (30.66%)") + 
  ylab("PC2 (20.93%)")

# Scoreplot with ggplot2: Sex
ggplot2::ggplot(fwFdPCA_scores, aes(x=Comp1, y=Comp2, color = Sex)) + 
  geom_point(shape = 16, alpha=0.50, size = 5)+
  scale_x_continuous(limits = c(-0.1, 0.1)) + 
  xlab("PC1 (30.66%)") + 
  ylab("PC2 (20.93%)")

# Loadings

# extract data
fwFdPCA_loadings = fwFdPCA$rotation

# edit rownames
rownames(fwFdPCA_loadings) = paste(rep(c("x","y"),20),
  c(paste(0, rep(1:9, each = 2),sep = ""), rep(10:20, each = 2)), sep = "")

# plot loadings for first PC
barplot(fwFdPCA_loadings[,1], las = 2, ylim = c(-0.6,0.6), 
  main = "PC1 Loadings")

# plot loadings for second PC
barplot(fwFdPCA_loadings[,2], las = 2, ylim = c(-0.6,0.6), 
  main = "PC2 Loadings")

# plot loadings for third PC
barplot(fwFdPCA_loadings[,3], las = 2, ylim = c(-0.6,0.6), 
  main = "PC3 Loadings")

# Eigenvalues

# sd of the principal components (square roots of the eigenvalues)
fwFdPCA_sdev = fwFdPCA$sdev

# eigenvalues
fwFdPCA_eigs = fwFdPCA$sdev^2 

# screeplot to choose how many PCs retain
screeplot(fwFdPCA, type = "l", main = "PCA screeplot")

# Variance explained and fidelity of representation
fwFdPCA_vrnc = fwFdPCA$sdev^2
fwFdPCA_vrnc_prop = fwFdPCA_vrnc/sum(fwFdPCA_vrnc)

fwFdPCA_fidelity = rbind(fwFdPCA_sdev,fwFdPCA_vrnc,
  fwFdPCA_vrnc_prop, cumsum(fwFdPCA_vrnc_prop))
  rownames(fwFdPCA_fidelity)=c("Standard deviation", 
  "Variance", "Proportion of Variance ", "Cumulative Proportion")
  colnames(fwFdPCA_fidelity)=colnames(fwFdPCA_scores[1:36])
  
round(fwFdPCA_fidelity,3)

# Barplot of explained variance
barplot(fwFdPCA_vrnc_prop*100, xlab='PC', ylab='Variance (%)', 
  main = "Explained variance",
  names.arg=1:length(fwFdPCA_vrnc_prop), las=1, 
  ylim=c(0, 40),
  col='gray')
  abline(h=0)

### PCA on the dataset without species outliers ----
fwSpPCA = gm.prcomp(fwSpOut_gpa$coords, transform = T)

# summary
summary(fwSpPCA)

#### Extract features from PCA ----

# Scores
fwSpPCA_scores = as.data.frame(fwSpPCA$x)
fwSpPCA_scores$Sex = fwSpOut_2d$Sex
fwSpPCA_scores$Species = fwSpOut_2d$Species

# Scoreplot with ggplot2: Species
ggplot2::ggplot(fwSpPCA_scores, aes(x=Comp1, y=Comp2, color = Species)) + 
  geom_point(shape = 16, alpha=0.50, size = 5)+
  scale_x_continuous(limits = c(-0.1, 0.1)) + 
  xlab("PC1 (30.74%)") + 
  ylab("PC2 (20.08%)")

# Scoreplot with ggplot2: Sex
ggplot2::ggplot(fwSpPCA_scores, aes(x=Comp1, y=Comp2, color = Sex)) + 
  geom_point(shape = 16, alpha=0.50, size = 5)+
  scale_x_continuous(limits = c(-0.1, 0.1)) + 
  xlab("PC1 (30.74%)") + 
  ylab("PC2 (20.08%)")

# Loadings

# extract data
fwSpPCA_loadings = fwSpPCA$rotation

# edit rownames
rownames(fwSpPCA_loadings) = paste(rep(c("x","y"),20),
  c(paste(0, rep(1:9, each = 2),sep = ""), 
  rep(10:20, each = 2)), sep = "")

# plot loadings for first PC
barplot(fwSpPCA_loadings[,1], las = 2, ylim = c(-0.6,0.6), 
  main = "PC1 Loadings")

# plot loadings for second PC
barplot(fwSpPCA_loadings[,2], las = 2, ylim = c(-0.6,0.6), 
  main = "PC2 Loadings")

# plot loadings for third PC
barplot(fwSpPCA_loadings[,3], las = 2, ylim = c(-0.6,0.6), 
  main = "PC3 Loadings")

# Eigenvalues

# sd of the principal components (square roots of the eigenvalues)
fwSpPCA_sdev = fwSpPCA$sdev

# eigenvalues
fwSpPCA_eigs = fwSpPCA$sdev^2 

# screeplot to choose how many PCs retain
screeplot(fwSpPCA, type = "l", main = "PCA screeplot")

# Variance explained and fidelity of representation
fwSpPCA_vrnc = fwSpPCA$sdev^2
fwSpPCA_vrnc_prop = fwSpPCA_vrnc/sum(fwSpPCA_vrnc)

fwSpPCA_fidelity= rbind(fwSpPCA_sdev,fwSpPCA_vrnc,fwSpPCA_vrnc_prop,
  cumsum(fwSpPCA_vrnc_prop))
  rownames(fwSpPCA_fidelity)=c("Standard deviation", "Variance",
  "Proportion of Variance ", "Cumulative Proportion")
  colnames(fwSpPCA_fidelity)=colnames(fwSpPCA_scores[1:36])

round(fwSpPCA_fidelity,3)

# Barplot of explained variance
barplot(fwSpPCA_vrnc_prop*100, xlab='PC', ylab='Variance (%)', 
  main = "Explained variance",
  names.arg=1:length(fwSpPCA_vrnc_prop), las=1, 
  ylim=c(0, 40),
  col='gray')
  abline(h=0)

### PCA on the dataset without sex outliers ----
fwSxPCA = gm.prcomp(fwSxOut_gpa$coords, transform = T)

# summary
summary(fwSxPCA)

#### Extract features from PCA ----

# Scores
fwSxPCA_scores = as.data.frame(fwSxPCA$x)
fwSxPCA_scores$Sex = fwSxOut_2d$Sex
fwSxPCA_scores$Species = fwSxOut_2d$Species

# Scoreplot with ggplot2: Species
ggplot2::ggplot(fwSxPCA_scores, aes(x=Comp1, y=Comp2, color = Species)) + 
  geom_point(shape = 16, alpha=0.50, size = 5)+
  scale_x_continuous(limits = c(-0.1, 0.1)) + 
  xlab("PC1 (29.16%)") + 
  ylab("PC2 (18.22%)")

# Scoreplot with ggplot2: Sex
ggplot2::ggplot(fwSxPCA_scores, aes(x=Comp1, y=Comp2, color = Sex)) + 
  geom_point(shape = 16, alpha=0.50, size = 5)+
  scale_x_continuous(limits = c(-0.1, 0.1)) + 
  xlab("PC1 (29.16%)") + 
  ylab("PC2 (18.22%)")

# Loadings

# extract data
fwSxPCA_loadings = fwSxPCA$rotation

# edit rownames
rownames(fwSxPCA_loadings) = paste(rep(c("x","y"),20),
  c(paste(0, rep(1:9, each = 2),sep = ""), rep(10:20, each = 2)), sep = "")

# plot loadings for first PC
barplot(fwSxPCA_loadings[,1], las = 2, ylim = c(-0.6,0.6), 
  main = "PC1 Loadings")

# plot loadings for second PC
barplot(fwSxPCA_loadings[,2], las = 2, ylim = c(-0.6,0.6), 
  main = "PC2 Loadings")

# plot loadings for third PC
barplot(fwSxPCA_loadings[,3], las = 2, ylim = c(-0.6,0.6), 
  main = "PC3 Loadings")

# Eigenvalues

# sd of the principal components (square roots of the eigenvalues)
fwSxPCA_sdev = fwSxPCA$sdev

# eigenvalues
fwSxPCA_eigs = fwSxPCA$sdev^2 

# screeplot to choose how many PCs retain
screeplot(fwSxPCA, type = "l", main = "PCA screeplot")

# Variance explained and fidelity of representation
fwSxPCA_vrnc = fwSxPCA$sdev^2
fwSxPCA_vrnc_prop = fwSxPCA_vrnc/sum(fwSxPCA_vrnc)

fwSxPCA_fidelity = rbind(fwSxPCA_sdev,fwSxPCA_vrnc,
  fwSxPCA_vrnc_prop, cumsum(fwSxPCA_vrnc_prop))
  rownames(fwSxPCA_fidelity)=c("Standard deviation", 
  "Variance", "Proportion of Variance ","Cumulative Proportion")
  colnames(fwSxPCA_fidelity)=colnames(fwSxPCA_scores[1:36])
  
round(fwSxPCA_fidelity,3)

# Barplot of explained variance
barplot(fwSxPCA_vrnc_prop*100, xlab='PC', ylab='Variance (%)', 
  main = "Explained variance",
  names.arg=1:length(fwSxPCA_vrnc_prop), las=1, 
  ylim=c(0, 40),
  col='gray')
abline(h=0)

## 06.02.02.01 SHAPE VARIATION -------------------------------------------------

### Shape Variation: max deformation across PCA axes ----

# Setting graphic paramteres with gridpar
GP1 <- gridPar(tar.pt.bg = 1, tar.pt.size = 1.5)
GP2 <- gridPar(pt.bg = 24 , pt.size = 2, tar.pt.bg = 607, tar.pt.size = 1.5)

# Shape deformation along PCs 

# Computing meanshapes for each dataset

# Full dataset
fwFd_ref <- mshape(forewing_gpa$coords)

# Without species outlier
fwSp_ref <- mshape(fwSpOut_gpa$coords)

# Without sex outlier
fwSx_ref <- mshape(fwSxOut_gpa$coords)

#### Maximum deformation as TPS ----  

##### Full dataset ----
pc1_fwFd_max = fwFdPCA$shapes$shapes.comp1$max  
pc2_fwFd_max = fwFdPCA$shapes$shapes.comp2$max
pc3_fwFd_max = fwFdPCA$shapes$shapes.comp3$max

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwFd_pc1_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwFd_ref, pc1_fwFd_max,
  mag=2, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="",sub = "")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwFd_pc2_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwFd_ref, pc2_fwFd_max,
  mag=2, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="",sub = "")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwFd_pc3_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwFd_ref, pc3_fwFd_max,
  mag=3, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="",sub = "")
dev.off()

##### Without species outlier ----
pc1_fwSp_max = fwSpPCA$shapes$shapes.comp1$max  
pc2_fwSp_max = fwSpPCA$shapes$shapes.comp2$max
pc3_fwSp_max = fwSpPCA$shapes$shapes.comp3$max

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwSp_pc1_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSp_ref, pc1_fwSp_max,
  mag=1, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum TPS",sub = "Max PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwSp_pc2_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSp_ref, pc2_fwSp_max,
  mag=1, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum TPS",sub = "Max PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwSp_pc3_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSp_ref, pc3_fwSp_max,
  mag=1, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum TPS",sub = "Max PC3")
dev.off()

##### Without sex outlier ----
pc1_fwSx_max = fwSpPCA$shapes$shapes.comp1$max  
pc2_fwSx_max = fwSpPCA$shapes$shapes.comp2$max
pc3_fwSx_max = fwSpPCA$shapes$shapes.comp3$max

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwSx_pc1_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSx_ref, pc1_fwSx_max,
  mag=1, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum TPS",sub = "Max PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwSx_pc2_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSx_ref, pc2_fwSx_max,
  mag=1, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum TPS",sub = "Max PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwSx_pc3_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSx_ref, pc3_fwSx_max,
  mag=1, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum TPS",sub = "Max PC3")
dev.off()

#### Maximum deformation as Vector  ----

##### Full dataset ----

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwFd_pc1_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwFd_ref, pc1_fwFd_max,
  mag=1.5, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="",sub = "")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwFd_pc2_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwFd_ref, pc2_fwFd_max,
  mag=2, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="",sub = "")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwFd_pc3_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwFd_ref, pc3_fwFd_max,
  mag=2, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="",sub = "")
dev.off()

##### Without species outlier ----

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwSp_pc1_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSp_ref, pc1_fwSp_max,
  mag=1, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum vectors",sub = "Max PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwSp_pc2_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSp_ref, pc2_fwSp_max,
  mag=1, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum vectors",sub = "Max PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwSp_pc3_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSp_ref, pc3_fwSp_max,
  mag=1, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum vectors",sub = "Max PC3")
dev.off()

##### Without sex outlier ----

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwSx_pc1_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSx_ref, pc1_fwSx_max,
  mag=1, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum vectors",sub = "Max PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwSx_pc2_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSx_ref, pc2_fwSx_max,
  mag=1, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum vectors",sub = "Max PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwSx_pc3_max", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSx_ref, pc3_fwSx_max,
  mag=1, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum vectors",sub = "Max PC3")
dev.off()

### Shape Variation: minimum deformation across PCA axes ----

#### Minimum deformation as TPS ----  

##### Full dataset ----
pc1_fwFd_min = fwFdPCA$shapes$shapes.comp1$min  
pc2_fwFd_min = fwFdPCA$shapes$shapes.comp2$min
pc3_fwFd_min = fwFdPCA$shapes$shapes.comp3$min

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwFd_pc1_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwFd_ref, pc1_fwFd_min,
  mag=1.5, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="",sub = "")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwFd_pc2_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwFd_ref, pc2_fwFd_min,
  mag=1, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum TPS",sub = "Min PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwFd_pc3_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwFd_ref, pc3_fwFd_min,
  mag=1, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum TPS",sub = "Min PC3")
dev.off()

##### Without species outlier ----
pc1_fwSp_min = fwSpPCA$shapes$shapes.comp1$min  
pc2_fwSp_min = fwSpPCA$shapes$shapes.comp2$min
pc3_fwSp_min = fwSpPCA$shapes$shapes.comp3$min

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwSp_pc1_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSp_ref, pc1_fwSp_min,
  mag=1, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum TPS",sub = "Min PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwSp_pc2_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSp_ref, pc2_fwSp_min,
  mag=1, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum TPS",sub = "Min PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwSp_pc3_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSp_ref, pc3_fwSp_min,
  mag=1, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum TPS",sub = "Min PC3")
dev.off()

##### Without sex outlier ----
pc1_fwSx_min = fwSpPCA$shapes$shapes.comp1$min  
pc2_fwSx_min = fwSpPCA$shapes$shapes.comp2$min
pc3_fwSx_min = fwSpPCA$shapes$shapes.comp3$min

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwSx_pc1_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSx_ref, pc1_fwSx_min,
  mag=1, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum TPS",sub = "Min PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwSx_pc2_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSx_ref, pc2_fwSx_min,
  mag=1, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum TPS",sub = "Min PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_fwSx_pc3_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSx_ref, pc3_fwSx_min,
  mag=1, 
  method = "TPS", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum TPS",sub = "Min PC3")
dev.off()

#### Minimum deformation as Vector  ----

##### Full dataset ----

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwFd_pc1_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwFd_ref, pc1_fwFd_min,
  mag=2, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="",sub = "")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwFd_pc2_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwFd_ref, pc2_fwFd_min,
  mag=2, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="",sub = "")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwFd_pc3_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwFd_ref, pc3_fwFd_min,
  mag=2, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="",sub = "")
dev.off()

##### Without species outlier ----

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwSp_pc1_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSp_ref, pc1_fwSp_min,
  mag=1, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum vectors",sub = "Min PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwSp_pc2_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSp_ref, pc2_fwSp_min,
  mag=1, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum vectors",sub = "Min PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwSp_pc3_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSp_ref, pc3_fwSp_min,
  mag=1, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum vectors",sub = "Min PC3")
dev.off()

##### Without sex outlier ----

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwSx_pc1_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSx_ref, pc1_fwSx_min,
  mag=1, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum vectors",sub = "Min PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwSx_pc2_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSx_ref, pc2_fwSx_min,
  mag=1, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum vectors",sub = "Min PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_fwSx_pc3_min", ".tif", sep=""), 
  res=1000, 
  compression = "lzw", 
  height=10, 
  width=10, 
  units="in")
plotRefToTarget(fwSx_ref, pc3_fwSx_min,
  mag=1, 
  method = "vector", 
  gridPars = GP1,
  links = fw_links)
title(main="Forewing Orthetrum vectors",sub = "Min PC3")
dev.off()