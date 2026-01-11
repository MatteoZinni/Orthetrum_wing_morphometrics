#==============================================================================#
# Script: data_preparation.R
# Author: Matteo Zinni
# Date: 2025-12-24
# Purpose: Collection of functions for processing geometric morphometric data
#==============================================================================#

#==============================================================================#
#' Convert GPA coordinates into a 2D dataframe
#'
#' This function takes a GPA object from geomorph and converts the coordinates
#' into a 2D dataframe, adding centroid size, species, and sex metadata.
#'
#' @param gpa A GPA object (output of gpagen)
#' @param specimen_data A dataframe containing specimen metadata with columns:
#'   Specimen, Species, Sex
#' @param n_landmarks Integer. Number of landmarks in the configuration
#'
#' @return A dataframe with 2D coordinates, C_size, Species, and Sex columns
#'
#' @examples
#' # forewing_2d <- prepare_2d_array(forewing_gpa, specimen, 20)
#'
#' @export
#==============================================================================#
prepare_2d_array <- function(gpa, specimen_data, n_landmarks) {
  # Convert GPA coordinates to 2D dataframe
  df <- as.data.frame(two.d.array(gpa$coords, sep = "."))
  
  # Set column names
  colnames(df) <- paste(
    rep(c("x","y"), n_landmarks),
    rep(sprintf("%0.2d", 1:n_landmarks), each = 2),
    sep = ""
  )
  
  # Add specimen metadata
  df$C_size  <- as.numeric(gpa$Csize)
  df$Species <- specimen_data[specimen_data$Specimen %in% rownames(df), "Species"]
  df$Sex     <- specimen_data[specimen_data$Specimen %in% rownames(df), "Sex"]
  
  return(df)
}
