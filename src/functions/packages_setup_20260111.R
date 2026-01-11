#==============================================================================#
# Script: install_packages.R
# Author: Matteo Zinni
# Date: 2026-01-11
# Purpose: Install and load required R packages for the project
#==============================================================================#

#==============================================================================#
#' Install and load required R packages
#'
#' This function reads a list of required packages from a text file
#' (default: required_packages.txt), installs any missing packages,
#' and loads all packages into the current R session.
#'
#' The function is designed for reproducible research workflows and
#' should be sourced at the beginning of analysis scripts.
#'
#' @param pkg_file Character. Path to a text file containing one package name
#'   per line.
#'
#' @return Invisibly returns TRUE if all packages are installed and loaded.
#'
#' @examples
#' # source("R/install_packages.R")
#' # install_packages()
#'
#' @export
#==============================================================================#
install_packages <- function(pkg_file = "required_packages.txt") {
  
  if (!file.exists(pkg_file)) {
    stop("Package list file not found: ", pkg_file)
  }
  
  # Read package names
  pkgs <- scan(pkg_file, what = character(), quiet = TRUE)
  pkgs <- pkgs[pkgs != ""]
  
  # Identify missing packages
  installed <- rownames(installed.packages())
  missing <- setdiff(pkgs, installed)
  
  # Install missing packages
  if (length(missing) > 0) {
    install.packages(missing)
  }
  
  # Load all packages
  invisible(lapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }))
  
  invisible(TRUE)
}