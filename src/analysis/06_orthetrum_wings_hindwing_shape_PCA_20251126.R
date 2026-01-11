## ========================================================================== ##
## Script:       06_orthetrum_wings_hindwing_shape_PCA_20251126.              ##
## Author:       Matteo Zinni                                                 ##
## Date:         2025-04-24                                                   ##
## Description:  A geometric morphometrics approach to investigate            ##
##               phylogenetic and taxonomic relationships among               ##
##               several skimmer species. Hindwing Shape Analysis - PCA.      ##
## ========================================================================== ##

# 06.01.01.01 SHAPE ANALYSIS ---------------------------------------------------

# Packages loading ----
invisible(lapply(packages, function(pkg) {
  library(pkg, character.only = TRUE)
}))

## 06.02.01.01 PCA ON SHAPE DATA -----------------------------------------------

### PCA on the full dataset ----
hwFdPCA = gm.prcomp(hindwing_gpa$coords, transform = T)

# summary
summary(hwFdPCA)

#### Extract features from PCA ----

# Scores
hwFdPCA_scores = as.data.frame(hwFdPCA$x)
hwFdPCA_scores$Sex = hindwing_2d$Sex
hwFdPCA_scores$Species = hindwing_2d$Species

# Scoreplot with ggplot2: Species
hw_scoreplot_sp <- ggplot(hwFdPCA_scores, aes(x=Comp1, y=Comp2, color = Species)) + 
  geom_point(shape = 16, alpha=0.50, size = 5)+
  scale_x_continuous(limits = c(-0.1, 0.1)) + 
  text = element_text(family = "Times") +
  xlab("PC1 (50.82%)") +
  ylab("PC2 (22.33%)")






install.packages("showtext")
library(extrafont)

library(showtext)

font_add("Times New Roman", "C:/Windows/Fonts/times.ttf")
showtext_auto()



font_import()

loadfonts(device = "win")

hw_scoreplot_sp <- ggplot2::ggplot(
  hwFdPCA_scores, 
  aes(x = Comp1, y = Comp2, color = Species)
) + 
  geom_point(shape = 16, alpha = 0.50, size = 5) +
  scale_x_continuous(limits = c(-0.1, 0.1)) +
  xlab("PC1 (50.82%)") + 
  ylab("PC2 (22.33%)") +
  theme_classic() +
  theme(text = element_text(family = "Times"))

hw_scoreplot_sp <- ggplot(
  hwFdPCA_scores,
  aes(x = Comp1, y = Comp2, color = Species)
) + 
  geom_point(shape = 16, alpha = 0.50, size = 5) +
  scale_x_continuous(limits = c(-0.1, 0.1)) +
  xlab("PC1 (50.82%)") + 
  ylab("PC2 (22.33%)") +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10)
  )



hw_scoreplot_sp <- ggplot(
  hwFdPCA_scores,
  aes(x = Comp1, y = Comp2, color = Species)
) +
  geom_point(shape = 16, alpha = 0.5, size = 10) +
  scale_x_continuous(limits = c(-0.1, 0.1)) +
  xlab("PC1 (50.82%)") +
  ylab("PC2 (22.33%)") +
  theme(
    text = element_text(family = "Times New Roman", size = 20),  # dimensione base più grande
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 22)
  )


ggsave("hw_scoreplot_sp.jpg",
       plot = hw_scoreplot_sp,
       dpi = 600,
       width = 15,
       height = 15,
       units = "cm")






# Carica Times New Roman
font_add("Times New Roman", "C:/Windows/Fonts/times.ttf")
showtext_auto()

# Il tuo plot originale con il tema che ti piace
hw_scoreplot_sp <- ggplot(
  hwFdPCA_scores,
  aes(x = Comp1, y = Comp2, color = Species)
) +
  geom_point(shape = 16, alpha = 0.50, size = 5) +
  scale_x_continuous(limits = c(-0.1, 0.1)) +
  xlab("PC1 (50.82%)") +
  ylab("PC2 (22.33%)") +
  theme(
    # Manteniamo il tema originale
    text = element_text(family = "Times New Roman", size = 36),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 26),
    legend.title = element_text(size = 28),
    legend.text = element_text(size = 26)
  )

# Esporta PDF vettoriale
ggsave(
  "hw_scoreplot_species.jpg",
  plot = hw_scoreplot_sp,
  width = 15, height = 15, units = "cm"
)





# Il tuo plot originale con il tema che ti piace
hw_scoreplot_sx <- ggplot(
  hwFdPCA_scores,
  aes(x = Comp1, y = Comp2, color = Sex)
) +
  geom_point(shape = 16, alpha = 0.50, size = 10) +
  scale_x_continuous(limits = c(-0.1, 0.1)) +
  xlab("PC1 (50.82%)") +
  ylab("PC2 (22.33%)") +
  theme(
    # Manteniamo il tema originale
    text = element_text(family = "Times New Roman", size = 36),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 26),
    legend.title = element_text(size = 28),
    legend.text = element_text(size = 26)
  )

# Esporta PDF vettoriale
ggsave(
  "hw_scoreplot_sex.jpg",
  plot = hw_scoreplot_sx,
  width = 15, height = 15, units = "cm"
)








# Scoreplot with ggplot2: Sex
hw_scoreplot_sx <- ggplot2::ggplot(hwFdPCA_scores, aes(x=Comp1, y=Comp2, color = Sex)) + 
  geom_point(shape = 16, alpha=0.50, size = 5)+
  scale_x_continuous(limits = c(-0.1, 0.1)) + 
  xlab("PC1 (50.82%)") +
  ylab("PC2 (22.33%)") 


ggsave("hw_scoreplot_sx.tiff",
       plot = hw_scoreplot_sx,
       dpi = 600,
       width = 15,
       height = 15,
       units = "cm",
       compression = "lzw")

# Test for sex separation along PC1 
kruskal.test(hwFdPCA_scores$Comp1 ~ hwFdPCA_scores$Sex)

# Test for sex separation along PC2
kruskal.test(hwFdPCA_scores$Comp2 ~ hwFdPCA_scores$Sex)

# Test for species separation along PC1 
kruskal.test(hwFdPCA_scores$Comp1 ~ hwFdPCA_scores$Species)

# Test for species separation along PC2
kruskal.test(hwFdPCA_scores$Comp2 ~ hwFdPCA_scores$Species)

# Post-hoc test for pairwise comparison
kwAllPairsDunnTest(Comp1 ~ Species, data = hwFdPCA_scores,
                   p.adjust.method = "bonferroni")

# Loadings

# extract data
hwFdPCA_loadings = hwFdPCA$rotation

hwFdPCA_loadings[,1]

# edit rownames
rownames(hwFdPCA_loadings) = paste(rep(c("x","y"),20),
                                   c(paste(0, rep(1:9, each = 2),sep = ""), rep(10:20, each = 2)), sep = "")

# plot loadings for first PC
barplot(hwFdPCA_loadings[,1], las = 2, ylim = c(-0.6,0.6), 
        main = "PC1 Loadings")

# plot loadings for second PC
barplot(hwFdPCA_loadings[,2], las = 2, ylim = c(-0.6,0.6), 
        main = "PC2 Loadings")

# plot loadings for third PC
barplot(hwFdPCA_loadings[,3], las = 2, ylim = c(-0.6,0.6), 
        main = "PC3 Loadings")

# Eigenvalues

# sd of the principal components (square roots of the eigenvalues)
hwFdPCA_sdev = hwFdPCA$sdev

# eigenvalues
hwFdPCA_eigs = hwFdPCA$sdev^2 

# screeplot to choose how many PCs retain
screeplot(hwFdPCA, type = "l", main = "PCA screeplot")

# Variance explained and fidelity of representation
hwFdPCA_vrnc = hwFdPCA$sdev^2
hwFdPCA_vrnc_prop = hwFdPCA_vrnc/sum(hwFdPCA_vrnc)

hwFdPCA_fidelity = rbind(hwFdPCA_sdev,hwFdPCA_vrnc,
                         hwFdPCA_vrnc_prop, cumsum(hwFdPCA_vrnc_prop))
rownames(hwFdPCA_fidelity)=c("Standard deviation", 
                             "Variance", "Proportion of Variance ", "Cumulative Proportion")
colnames(hwFdPCA_fidelity)=colnames(hwFdPCA_scores[1:36])

round(hwFdPCA_fidelity,3)

# Barplot of explained variance
barplot(hwFdPCA_vrnc_prop*100, xlab='PC', ylab='Variance (%)', 
        main = "Explained variance",
        names.arg=1:length(hwFdPCA_vrnc_prop), las=1, 
        ylim=c(0, 40),
        col='gray')
abline(h=0)

### PCA on the dataset without species outliers ----
hwSpPCA = gm.prcomp(hwSpOut_gpa$coords, transform = T)

# summary
summary(hwSpPCA)

#### Extract features from PCA ----

# Scores
hwSpPCA_scores = as.data.frame(hwSpPCA$x)
hwSpPCA_scores$Sex = hwSpOut_2d$Sex
hwSpPCA_scores$Species = hwSpOut_2d$Species

# Scoreplot with ggplot2: Species
ggplot2::ggplot(hwSpPCA_scores, aes(x=Comp1, y=Comp2, color = Species)) + 
  geom_point(shape = 16, alpha=0.50, size = 5)+
  scale_x_continuous(limits = c(-0.1, 0.1)) + 
  xlab("PC1 (30.74%)") + 
  ylab("PC2 (20.08%)")

# Scoreplot with ggplot2: Sex
ggplot2::ggplot(hwSpPCA_scores, aes(x=Comp1, y=Comp2, color = Sex)) + 
  geom_point(shape = 16, alpha=0.50, size = 5)+
  scale_x_continuous(limits = c(-0.1, 0.1)) + 
  xlab("PC1 (30.74%)") + 
  ylab("PC2 (20.08%)")

# Loadings

# extract data
hwSpPCA_loadings = hwSpPCA$rotation

# edit rownames
rownames(hwSpPCA_loadings) = paste(rep(c("x","y"),20),
                                   c(paste(0, rep(1:9, each = 2),sep = ""), 
                                     rep(10:20, each = 2)), sep = "")

# plot loadings for first PC
barplot(hwSpPCA_loadings[,1], las = 2, ylim = c(-0.6,0.6), 
        main = "PC1 Loadings")

# plot loadings for second PC
barplot(hwSpPCA_loadings[,2], las = 2, ylim = c(-0.6,0.6), 
        main = "PC2 Loadings")

# plot loadings for third PC
barplot(hwSpPCA_loadings[,3], las = 2, ylim = c(-0.6,0.6), 
        main = "PC3 Loadings")

# Eigenvalues

# sd of the principal components (square roots of the eigenvalues)
hwSpPCA_sdev = hwSpPCA$sdev

# eigenvalues
hwSpPCA_eigs = hwSpPCA$sdev^2 

# screeplot to choose how many PCs retain
screeplot(hwSpPCA, type = "l", main = "PCA screeplot")

# Variance explained and fidelity of representation
hwSpPCA_vrnc = hwSpPCA$sdev^2
hwSpPCA_vrnc_prop = hwSpPCA_vrnc/sum(hwSpPCA_vrnc)

hwSpPCA_fidelity= rbind(hwSpPCA_sdev,hwSpPCA_vrnc,hwSpPCA_vrnc_prop,
                        cumsum(hwSpPCA_vrnc_prop))
rownames(hwSpPCA_fidelity)=c("Standard deviation", "Variance",
                             "Proportion of Variance ", "Cumulative Proportion")
colnames(hwSpPCA_fidelity)=colnames(hwSpPCA_scores[1:36])

round(hwSpPCA_fidelity,3)

# Barplot of explained variance
barplot(hwSpPCA_vrnc_prop*100, xlab='PC', ylab='Variance (%)', 
        main = "Explained variance",
        names.arg=1:length(hwSpPCA_vrnc_prop), las=1, 
        ylim=c(0, 40),
        col='gray')
abline(h=0)

### PCA on the dataset without sex outliers ----
hwSxPCA = gm.prcomp(hwSxOut_gpa$coords, transform = T)

# summary
summary(hwSxPCA)

#### Extract features from PCA ----

# Scores
hwSxPCA_scores = as.data.frame(hwSxPCA$x)
hwSxPCA_scores$Sex = hwSxOut_2d$Sex
hwSxPCA_scores$Species = hwSxOut_2d$Species

# Scoreplot with ggplot2: Species
ggplot2::ggplot(hwSxPCA_scores, aes(x=Comp1, y=Comp2, color = Species)) + 
  geom_point(shape = 16, alpha=0.50, size = 5)+
  scale_x_continuous(limits = c(-0.1, 0.1)) + 
  xlab("PC1 (29.16%)") + 
  ylab("PC2 (18.22%)")

# Scoreplot with ggplot2: Sex
ggplot2::ggplot(hwSxPCA_scores, aes(x=Comp1, y=Comp2, color = Sex)) + 
  geom_point(shape = 16, alpha=0.50, size = 5)+
  scale_x_continuous(limits = c(-0.1, 0.1)) + 
  xlab("PC1 (29.16%)") + 
  ylab("PC2 (18.22%)")

# Loadings

# extract data
hwSxPCA_loadings = hwSxPCA$rotation

# edit rownames
rownames(hwSxPCA_loadings) = paste(rep(c("x","y"),20),
                                   c(paste(0, rep(1:9, each = 2),sep = ""), rep(10:20, each = 2)), sep = "")

# plot loadings for first PC
barplot(hwSxPCA_loadings[,1], las = 2, ylim = c(-0.6,0.6), 
        main = "PC1 Loadings")

# plot loadings for second PC
barplot(hwSxPCA_loadings[,2], las = 2, ylim = c(-0.6,0.6), 
        main = "PC2 Loadings")

# plot loadings for third PC
barplot(hwSxPCA_loadings[,3], las = 2, ylim = c(-0.6,0.6), 
        main = "PC3 Loadings")

# Eigenvalues

# sd of the principal components (square roots of the eigenvalues)
hwSxPCA_sdev = hwSxPCA$sdev

# eigenvalues
hwSxPCA_eigs = hwSxPCA$sdev^2 

# screeplot to choose how many PCs retain
screeplot(hwSxPCA, type = "l", main = "PCA screeplot")

# Variance explained and fidelity of representation
hwSxPCA_vrnc = hwSxPCA$sdev^2
hwSxPCA_vrnc_prop = hwSxPCA_vrnc/sum(hwSxPCA_vrnc)

hwSxPCA_fidelity = rbind(hwSxPCA_sdev,hwSxPCA_vrnc,
                         hwSxPCA_vrnc_prop, cumsum(hwSxPCA_vrnc_prop))
rownames(hwSxPCA_fidelity)=c("Standard deviation", 
                             "Variance", "Proportion of Variance ","Cumulative Proportion")
colnames(hwSxPCA_fidelity)=colnames(hwSxPCA_scores[1:36])

round(hwSxPCA_fidelity,3)

# Barplot of explained variance
barplot(hwSxPCA_vrnc_prop*100, xlab='PC', ylab='Variance (%)', 
        main = "Explained variance",
        names.arg=1:length(hwSxPCA_vrnc_prop), las=1, 
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
hwFd_ref <- mshape(hindwing_gpa$coords)

# Without species outlier
hwSp_ref <- mshape(hwSpOut_gpa$coords)

# Without sex outlier
hwSx_ref <- mshape(hwSxOut_gpa$coords)

#### Maximum deformation as TPS ----  

##### Full dataset ----
pc1_hwFd_max = hwFdPCA$shapes$shapes.comp1$max  
pc2_hwFd_max = hwFdPCA$shapes$shapes.comp2$max
pc3_hwFd_max = hwFdPCA$shapes$shapes.comp3$max

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwFd_pc1_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwFd_ref, pc1_hwFd_max,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Max PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwFd_pc2_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwFd_ref, pc2_hwFd_max,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Max PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwFd_pc3_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwFd_ref, pc3_hwFd_max,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Max PC3")
dev.off()

##### Without species outlier ----
pc1_hwSp_max = hwSpPCA$shapes$shapes.comp1$max  
pc2_hwSp_max = hwSpPCA$shapes$shapes.comp2$max
pc3_hwSp_max = hwSpPCA$shapes$shapes.comp3$max

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwSp_pc1_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSp_ref, pc1_hwSp_max,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Max PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwSp_pc2_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSp_ref, pc2_hwSp_max,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Max PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwSp_pc3_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSp_ref, pc3_hwSp_max,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Max PC3")
dev.off()

##### Without sex outlier ----
pc1_hwSx_max = hwSpPCA$shapes$shapes.comp1$max  
pc2_hwSx_max = hwSpPCA$shapes$shapes.comp2$max
pc3_hwSx_max = hwSpPCA$shapes$shapes.comp3$max

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwSx_pc1_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSx_ref, pc1_hwSx_max,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Max PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwSx_pc2_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSx_ref, pc2_hwSx_max,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Max PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwSx_pc3_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSx_ref, pc3_hwSx_max,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Max PC3")
dev.off()

#### Maximum deformation as Vector  ----

##### Full dataset ----

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwFd_pc1_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwFd_ref, pc1_hwFd_max,
                mag=1, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="",sub = "")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwFd_pc2_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwFd_ref, pc2_hwFd_max,
                mag=1, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="",sub = "")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwFd_pc3_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwFd_ref, pc3_hwFd_max,
                mag=1, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum vectors",sub = "Max PC3")
dev.off()

##### Without species outlier ----

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwSp_pc1_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSp_ref, pc1_hwSp_max,
                mag=1, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum vectors",sub = "Max PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwSp_pc2_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSp_ref, pc2_hwSp_max,
                mag=1, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum vectors",sub = "Max PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwSp_pc3_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSp_ref, pc3_hwSp_max,
                mag=1, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum vectors",sub = "Max PC3")
dev.off()

##### Without sex outlier ----

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwSx_pc1_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSx_ref, pc1_hwSx_max,
                mag=1, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum vectors",sub = "Max PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwSx_pc2_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSx_ref, pc2_hwSx_max,
                mag=1, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum vectors",sub = "Max PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwSx_pc3_max", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSx_ref, pc3_hwSx_max,
                mag=1, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum vectors",sub = "Max PC3")
dev.off()

### Shape Variation: minimum deformation across PCA axes ----

#### Minimum deformation as TPS ----  

##### Full dataset ----
pc1_hwFd_min = hwFdPCA$shapes$shapes.comp1$min  
pc2_hwFd_min = hwFdPCA$shapes$shapes.comp2$min
pc3_hwFd_min = hwFdPCA$shapes$shapes.comp3$min

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwFd_pc1_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwFd_ref, pc1_hwFd_min,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Min PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwFd_pc2_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwFd_ref, pc2_hwFd_min,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Min PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwFd_pc3_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwFd_ref, pc3_hwFd_min,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Min PC3")
dev.off()

##### Without species outlier ----
pc1_hwSp_min = hwSpPCA$shapes$shapes.comp1$min  
pc2_hwSp_min = hwSpPCA$shapes$shapes.comp2$min
pc3_hwSp_min = hwSpPCA$shapes$shapes.comp3$min

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwSp_pc1_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSp_ref, pc1_hwSp_min,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Min PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwSp_pc2_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSp_ref, pc2_hwSp_min,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Min PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwSp_pc3_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSp_ref, pc3_hwSp_min,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Min PC3")
dev.off()

##### Without sex outlier ----
pc1_hwSx_min = hwSpPCA$shapes$shapes.comp1$min  
pc2_hwSx_min = hwSpPCA$shapes$shapes.comp2$min
pc3_hwSx_min = hwSpPCA$shapes$shapes.comp3$min

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwSx_pc1_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSx_ref, pc1_hwSx_min,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Min PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwSx_pc2_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSx_ref, pc2_hwSx_min,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Min PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_tps_hwSx_pc3_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSx_ref, pc3_hwSx_min,
                mag=1, 
                method = "TPS", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum TPS",sub = "Min PC3")
dev.off()

#### Minimum deformation as Vector  ----

##### Full dataset ----

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwFd_pc1_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwFd_ref, pc1_hwFd_min,
                mag=1, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="",sub = "")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwFd_pc2_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwFd_ref, pc2_hwFd_min,
                mag=1.5, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="",sub = "")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwFd_pc3_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwFd_ref, pc3_hwFd_min,
                mag=1.5, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="",sub = "")
dev.off()

##### Without species outlier ----

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwSp_pc1_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSp_ref, pc1_hwSp_min,
                mag=1, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum vectors",sub = "Min PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwSp_pc2_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSp_ref, pc2_hwSp_min,
                mag=1, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum vectors",sub = "Min PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwSp_pc3_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSp_ref, pc3_hwSp_min,
                mag=1, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum vectors",sub = "Min PC3")
dev.off()

##### Without sex outlier ----

# PC1
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwSx_pc1_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSx_ref, pc1_hwSx_min,
                mag=1, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum vectors",sub = "Min PC1")
dev.off()

# PC2
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwSx_pc2_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSx_ref, pc2_hwSx_min,
                mag=1, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum vectors",sub = "Min PC2")
dev.off()

# PC3
plot.new()
tiff(paste0(output_figs,"orthetrum_shape_vec_hwSx_pc3_min", ".tif", sep=""), 
     res=1000, 
     compression = "lzw", 
     height=10, 
     width=10, 
     units="in")
plotRefToTarget(hwSx_ref, pc3_hwSx_min,
                mag=1, 
                method = "vector", 
                gridPars = GP1,
                links = hw_links)
title(main="Hindwing Orthetrum vectors",sub = "Min PC3")
dev.off()