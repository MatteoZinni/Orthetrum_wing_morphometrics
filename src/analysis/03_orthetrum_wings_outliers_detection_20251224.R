## ========================================================================== ##
## Script:      03_orthetrum_wings_outlier_detection_20251224.R
## Author:      Matteo Zinni
## Date:        2025-12-24
## Description: Detect and remove outliers in Orthetrum wing landmark data
##              using geometric morphometric methods.
## ========================================================================== ##

# 03.01.01.01 OUTLIER DETECTION ------------------------------------------------

# Packages loading ----
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}

# Sourcing the function
source(here::here("src", "functions", "packages_setup_20260111.R"))

# Lunch the function
install_packages()

## 03.01.02.01 LOOKING FOR OUTLIER ---------------------------------------------
# NOTE: inspect.outliers = TRUE enables interactive inspection;
# final outlier lists are manually curated below.

### Outliers in the whole dataset ----

# Forewing
plotOutliers(forewing_gpa$coords, inspect.outliers = TRUE)

# Hindwing
plotOutliers(hindwing_gpa$coords, inspect.outliers = TRUE)

### Outliers within groups ----

#### Forewing ----

# Species
plotOutliers(forewing_gpa$coords, groups = specimen$Species)

# Sex
plotOutliers(forewing_gpa$coords, groups = specimen$Sex)

# List of outliers within species for forewing landmarks
fw_out_sp <- c("ORTALB010", "ORTBRU007", "ORTCAN003", 
               "LIBDEP005", "LIBDEP003", "ORTSAB012")

# List of outliers within sex for forewing landmarks
fw_out_sx <- c("ORTJUL001", "ORTCAN003", "ORTGLA002",
               "ORTJUL011", "ORTJUL010", "ORTALB010", 
               "ORTBRU007", "ORTJUL006", "ORTALB007")

#### Hindwing ----

# Species
plotOutliers(hindwing_gpa$coords, groups = specimen$Species)

# Sex
plotOutliers(hindwing_gpa$coords, groups = specimen$Sex)

# List of outliers within species for hindwing landmarks
hw_out_sp <- c("ORTALB003", "LIBDEP003", "ORTGLA005", 
               "ORTGLA002", "ORTJUL010", "ORTJUL011", 
               "ORTSAB010", "ORTSAB011", "ORTTRI001")

# List of outliers within sex for hindwing landmarks
hw_out_sx <- c("ORTJUL005", "ORTTRI010", "ORTALB004", "ORTTRI001", 
               "ORTJUL011", "ORTJUL010", "ORTALB003")

# List t store outliers
outliers <- list(
  forewing = list(
    species = fw_out_sp,
    sex     = fw_out_sx
  ),
  hindwing = list(
    species = hw_out_sp,
    sex     = hw_out_sx
  )
)

## 03.01.03.01 OUTLIER REMOVAL -------------------------------------------------

### Removing outliers within species ----

# Forewing landmarks data
fwSpOut <- forewing[, , !(dimnames(forewing)[[3]] %in% outliers[[1]]$species)]

# Hindwing landmarks data
hwSpOut <- hindwing[, , !(dimnames(hindwing)[[3]] %in% outliers[[2]]$species)]

### Removing outliers within sexes ----

# Forewing landmarks data
fwSxOut <- forewing[, , !(dimnames(forewing)[[3]] %in% outliers[[1]]$sex)]

# Hindwing landmarks data
hwSxOut <- hindwing[, , !(dimnames(hindwing)[[3]] %in% outliers[[2]]$sex)]