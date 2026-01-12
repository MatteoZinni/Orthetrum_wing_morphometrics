## ========================================================================== ##
## Script:      04_orthetrum_wings_data_preparation_20260112.R
## Author:      Matteo Zinni
## Date:        2026-01-12
## Description: Prepare Orthetrum wing landmark data for analysis, including
##              GPA, outlier-filtered datasets, 2D array export, and creation
##              of geomorph data frames.
## ========================================================================== ##

# 04.01.01.01 DATA PREPARATION -------------------------------------------------

# Packages loading ----
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}

# Sourcing the function
source(here::here("src", "functions", "packages_setup_20260111.R"))

# Lunch the function
install_packages()

# Load all custom functions from the "data_preparation.R" script.
# This allows us to use the function prepare_2d_array()
source(here::here("src", "functions", "data_preparation_20251224.R"))

## 04.02.02.01 PREPARE FULL DATASET --------------------------------------------

### 04.02.02.02 Convert Data into 2D Array -----

#### Forewing data ----

# Convert data into dataframe
forewing_2d <- prepare_2d_array(forewing_gpa, specimen, 20)

# Export forewing data
write.csv(
  forewing_2d,
  file = here::here(
    "data", "processed",
    paste0("forewing_2d_", format(Sys.time(), "%Y%m%d"), ".csv")
  ),
  row.names = TRUE
)

#### Hindwing data ----

# Convert data into dataframe
hindwing_2d <- prepare_2d_array(hindwing_gpa, specimen, 22)

# Export hindwing data
write.csv(
  hindwing_2d,
  file = here::here(
    "data", "processed",
    paste0("hindwng_2d_", format(Sys.time(), "%Y%m%d"), ".csv")
  ),
  row.names = TRUE
)

## 04.03.03.01 PREPARE DATA WITH OUTLIER WITHIN SPECIES REMOVED ----------------

### Forewing data ----

# Generalised Procrustes Analysis without outliers
fwSpOut_gpa <- gpagen(fwSpOut, curves = NULL, surfaces = NULL, PrinAxes = TRUE, 
                      max.iter = NULL, ProcD = TRUE, Proj = TRUE, 
                      print.progress = FALSE)

### Hindwing data ----

# Generalised Procrustes Analysis without outliers
hwSpOut_gpa <- gpagen(hwSpOut, curves = NULL, surfaces = NULL, PrinAxes = TRUE, 
                      max.iter = NULL, ProcD = TRUE, Proj = TRUE, 
                      print.progress = FALSE)

#### 04.03.03.02 Convert Data into 2D Array ----

##### Forewing data ----

# Convert data into dataframe
fwSpOut_2d <- prepare_2d_array(fwSpOut_gpa, specimen, 20)

# Export forewing data
write.csv(
  fwSpOut_2d,
  file = here::here(
    "data", "processed",
    paste0("fwSpOut_2d_", format(Sys.time(), "%Y%m%d"), ".csv")
  ),
  row.names = TRUE
)

##### Hindwing data ----

# Convert data into dataframe
hwSpOut_2d <- prepare_2d_array(hwSpOut_gpa, specimen, 22)

# Export
write.csv(
  hwSpOut_2d,
  file = here::here(
    "data", "processed",
    paste0("hwSpOut_2d_", format(Sys.time(), "%Y%m%d"), ".csv")
  ),
  row.names = TRUE
)

## 04.03.04.01 PREPARE DATA WITH OUTLIER WITHIN SEXES REMOVED ------------------

### Forewing data ----
fwSxOut_gpa <- gpagen(fwSxOut, curves = NULL, surfaces = NULL, PrinAxes = TRUE, 
                      max.iter = NULL, ProcD = TRUE, Proj = TRUE, 
                      print.progress = FALSE)

### Hindwing data ----
hwSxOut_gpa <- gpagen(hwSxOut, curves = NULL, surfaces = NULL, PrinAxes = TRUE, 
                      max.iter = NULL, ProcD = TRUE, Proj = TRUE, 
                      print.progress = FALSE)

#### 04.03.04.02 Convert Data into 2D Array ----

##### Forewing data ----

# Convert data into dataframe
fwSxOut_2d <- prepare_2d_array(fwSxOut_gpa, specimen, 20)

# Export
write.csv(
  fwSxOut_2d,
  file = here::here(
    "data", "processed",
    paste0("fwSxOut_2d_", format(Sys.time(), "%Y%m%d"), ".csv")
  ),
  row.names = TRUE
)

##### Hindwing data ----

# Convert data into dataframe
hwSxOut_2d <- prepare_2d_array(hwSxOut_gpa, specimen, 22)

# Export
write.csv(
  hwSxOut_2d,
  file = here::here(
    "data", "processed",
    paste0("hwSxOut_2d_", format(Sys.time(), "%Y%m%d"), ".csv")
  ),
  row.names = TRUE
)

## 04.04.01.01 CREATION OF GEOMORPH DATA FRAME ---------------------------------

### Full dataset ----

# Forewing
forewing_gdf <- geomorph.data.frame(forewing_gpa, 
                                    Species = specimen$Species, 
                                    Sex     = specimen$Sex)

# Hindwing
hindwing_gdf <- geomorph.data.frame(hindwing_gpa, 
                                    Species = specimen$Species, 
                                    Sex     = specimen$Sex)

### Without species outliers ----

# Forewing
fwSpOut_gdf <- geomorph.data.frame(fwSpOut_gpa, 
                                   Species = fwSpOut_2d$Species, 
                                   Sex     = fwSpOut_2d$Sex)

# Hindwing
hwSpOut_gdf <- geomorph.data.frame(hwSpOut_gpa, 
                                   Species = hwSpOut_2d$Species, 
                                   Sex     = hwSpOut_2d$Sex)

### Without sex outliers ----

# Forewing
fwSxOut_gdf <- geomorph.data.frame(fwSxOut_gpa, 
                                   Species = fwSxOut_2d$Species, 
                                   Sex     = fwSxOut_2d$Sex)

# Hindwing
hwSxOut_gdf <- geomorph.data.frame(hwSxOut_gpa, 
                                   Species = hwSxOut_2d$Species, 
                                   Sex     = hwSxOut_2d$Sex)

## 04.05.01.01 VISUALIZATION AND LANDMARK CHECKS ------------------------------

# Exploratory plots and figures for landmark configuration inspection.
# These plots are used for quality control and figure generation for the paper.

## Superimposed landmarks ----

### Forewing data ----

# Full dataset
plotAllSpecimens(
  forewing_gpa$coords,
  mean = TRUE,
  links = NULL,
  label = F,
  plot_param = list()
)

# Without species outlier
plotAllSpecimens(
  fwSpOut_gpa$coords,
  mean = TRUE,
  links = NULL,
  label = F,
  plot_param = list()
)

# Without sex outlier
plotAllSpecimens(
  fwSxOut_gpa$coords,
  mean = TRUE,
  links = NULL,
  label = F,
  plot_param = list()
)

### Hindwing data ----

# Full dataset
plotAllSpecimens(
  hindwing_gpa$coords,
  mean = TRUE,
  links = NULL,
  label = F,
  plot_param = list()
)

# Without species outlier
plotAllSpecimens(
  hwSpOut_gpa$coords,
  mean = TRUE,
  links = NULL,
  label = F,
  plot_param = list()
)

# Without sex outlier
plotAllSpecimens(
  hwSpOut_gpa$coords,
  mean = TRUE,
  links = NULL,
  label = F,
  plot_param = list()
)

## Landmark positions and numebring ----

# Forewing mean shape
plot(mshape(forewing_gpa$coords), links = fw_links)


# Hindwing mean shape
plot(mshape(hindwing_gpa$coords), links = hw_links)


### Figure for manuscript (high resolution TIFF) ----

# Export High Resolution TIFF
tiff(
  filename = here::here("output", "figures", "forewing_plot.tiff"),
  width = 1692,
  height = 875,
  units = "px",
  res = 600,
  compression = "lzw"
)

# Hindwing mean shape
plot(mshape(hindwing_gpa$coords), links = hw_links)
dev.off()

# Export High Resolution TIFF
tiff(
  filename = here::here("output", "figures", "hindwing_plot.tiff"),
  width = 1692,
  height = 875,
  units = "px",
  res = 600,
  compression = "lzw"
)

# Hindwing mean shape
plot(mshape(hindwing_gpa$coords), links = hw_links)
dev.off()