## ========================================================================== ##
## Script:       02_orthetrum_wings_data_editing_20251223.R
## Author:       Matteo Zinni
## Date:         2025-12-23
## Description: Edit landmark and specimen data for Orthetrum wings
##              as part of a geometric morphometrics workflow.
## ========================================================================== ##

# 02.01.01.01 DATA EDITING  ----------------------------------------------------

# Packages loading ----
message("Loading packages: ", paste(packages, collapse = ", "))
invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

## 02.01.02.01 Edit data ----

### Specimen data ----

# Cast Species to factor 
specimen$Species <- as.factor(specimen$Species)

# Cast Sex to factor 
specimen$Sex <- as.factor(specimen$Sex)

### Landmarks data ----

#### Link Landmarks for Wireframe ----       

# Link forewing landmarks to build wireframe
fw_link_start <- c(1,2,3,4,5,6,7,8,9,10,12,13,14,16,17,18,19,14)
fw_link_end   <- c(2,3,4,5,6,7,8,9,10,11,13,14,12,17,18,19,16,12)

# Create matrix for forewing links
fw_links <- cbind(fw_link_start, fw_link_end)
fw_links <- as.matrix(fw_links)

# Link hindwing landmarks to build wireframe
hw_link_start <- c(1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,18,19,20,21)
hw_link_end   <- c(2,3,4,5,6,7,8,9,10,11,12,13,15,16,14,19,20,21,18)

# Create matrix for hindwing links
hw_links <- cbind(hw_link_start, hw_link_end)
hw_links <- as.matrix(hw_links)

#### Rotate flipped specimen ----

# Forewing

# Rotate landmarks
fw_rotated_lnmrks <- rotate.coords(forewing[,,c("ORTCOE009", "ORTTRI009",
                                                "ORTTRI010", "ORTJUL003",
                                                "ORTJUL009", "ORTPRU009")], 
                                   "flipY")

# List of specimen                                         
fw_rotated_spcmn <- c("ORTCOE009", "ORTTRI009", "ORTTRI010", 
                      "ORTJUL003", "ORTJUL009", "ORTPRU009")      

# Replace with flipped specimen
forewing[,,fw_rotated_spcmn] <- fw_rotated_lnmrks

# Hindwing

# Rotate landmarks
hw_rotated_lnmrks <- rotate.coords(hindwing[,,c("ORTGLA004", "ORTGLA006",
                                                "ORTPRU007", "ORTJUL009",
                                                "ORTJUL003", "ORTSAB003",
                                                "ORTTRI009")], "flipY")

# List of specimen                                         
hw_rotated_spcmn <- c("ORTGLA004", "ORTGLA006", "ORTPRU007", "ORTJUL009",
                      "ORTJUL003", "ORTSAB003", "ORTTRI009")

# Replace with flipped specimen
hindwing[,,hw_rotated_spcmn] <- hw_rotated_lnmrks

## 02.01.03.01 DATA STANDARDIZATION --------------------------------------------

### Generalised Procrustes Analysis ----

# forewing data
message("Running GPA for forewing...")
forewing_gpa <- gpagen(forewing, curves = NULL, 
                       surfaces = NULL, PrinAxes = TRUE, 
                       max.iter = NULL, ProcD = TRUE, 
                       Proj = TRUE, print.progress = F)
message("Forewing GPA completed.")

# hindwing data 
message("Running GPA for hindwing...")
hindwing_gpa <- gpagen(hindwing, curves = NULL, 
                       surfaces = NULL, PrinAxes = TRUE, 
                       max.iter = NULL, ProcD = TRUE, 
                       Proj = TRUE, print.progress = F)
message("Forewing GPA completed.")

### Check landmarks data structure ----

# number of landmarks
forewing_gpa$p
hindwing_gpa$p 

print(forewing_gpa$p)
print(hindwing_gpa$p)

# number of dimension
forewing_gpa$k
hindwing_gpa$k

print(forewing_gpa$k)
print(hindwing_gpa$k)

# number of landmarks
forewing_gpa$consensus
hindwing_gpa$consensus

print(forewing_gpa$consensus)
print(hindwing_gpa$consensus)

# centroid size
forewing_gpa$Csize
hindwing_gpa$Csize

print(forewing_gpa$Csize)
print(hindwing_gpa$Csize)

### Plot superimposed landmarks ----

# Forewing
plot(forewing_gpa)
title(main = "Forewing \n Superimposed landmarks")

# Hindwing
plot(hindwing_gpa)
title(main = "Hindwing \n Superimposed landmarks")