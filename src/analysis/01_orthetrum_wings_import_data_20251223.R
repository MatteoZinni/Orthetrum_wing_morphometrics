## ========================================================================== ##
## Script:      01_orthetrum_wings_import_data_2025123.R
## Author:      Matteo Zinni
## Date:        2025-12-23
## Description: Import landmark and specimen data for Orthetrum wings
##              as part of a geometric morphometrics workflow.
## ========================================================================== ##

# 01.01.01.01 IMPORT DATA ------------------------------------------------------

# Packages loading ----
message("Loading packages: ", paste(packages, collapse = ", "))
invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

## Landmarks data ----

# Reading tps file for forewing data
forewing <- readland.tps(file.path(paths$data_base, "orthetrum_forewing.tps"),
                         specID = "ID",
                         readcurves = FALSE,
                         warnmsg = TRUE)

# Reading tps file for hindwing data
hindwing <- readland.tps(file.path(paths$data_base, "orthetrum_hindwing.tps"),
                         specID = "ID",
                         readcurves = FALSE,
                         warnmsg = TRUE)

## Specimen data ----
specimen <- read.csv(file.path(paths$data_base, "orthetrum_specimen.csv"),
                     header = TRUE, 
                     dec = ",")

## Phylogenetic tree ----
orthetrum_phy <- ape::read.nexus(file.path(paths$data_base, 
                                           "orthetrum_phylogenetics.nex"),
                                 tree.names = NULL,
                                 force.multi = FALSE)