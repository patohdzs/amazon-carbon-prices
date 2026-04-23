# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: MASTERFILE SCRIPT TO SOURCE ALL CARBON CAPTURE CURVE ANALYSIS SCRIPTS
# AUTHOR: LEONARDO DE CAMPOS GOMES
#


library(tictoc)

# Start timer
tic(msg = "carbon_capture_curves/_masterfile.R script", log = TRUE)

# Prepare data for carbon capture curves
source("rsrc/analysis/carbon_capture_curves/01_prepare_data.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

# Analyze and plot carbon capture curves
source("rsrc/analysis/carbon_capture_curves/02_analysis.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

# End timer
toc(log = TRUE)
