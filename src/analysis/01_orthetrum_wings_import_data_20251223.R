## ========================================================================== ##
## Script:      01_orthetrum_wings_import_data_20260111.R
## Author:      Matteo Zinni
## Date:        2026-01-11
## Description: Import landmark and specimen data for Orthetrum wings
##              as part of a geometric morphometrics workflow.
## ========================================================================== ##

# 01.01.01.01 IMPORT DATA ------------------------------------------------------

# Packages loading ----
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}

# Sourcing the function
source(here::here("src", "functions", "packages_setup_20260111.R"))

# Lunch the function
install_packages()

## Landmarks data ----

# Reading tps file for forewing data
forewing <- readland.tps(
  file = here::here("data", "raw", "orthetrum_forewing.tps"),
  specID = "ID",
  readcurves = FALSE,
  warnmsg = TRUE
)

# Reading tps file for hindwing data
hindwing <- readland.tps(
  file = here::here("data", "raw", "orthetrum_hindwing.tps"),
  specID = "ID",
  readcurves = FALSE,
  warnmsg = TRUE
)

## Specimen data ----
specimen <- read.csv(
  here::here("data", "raw", "orthetrum_specimen.csv"),
  header = TRUE,
  dec = ","
)