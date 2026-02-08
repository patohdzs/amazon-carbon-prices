# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: MASTERFILE SCRIPT TO SOURCE ALL CALIBRATION SCRIPTS
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

library(tictoc)

# Start timer
tic(msg = "_masterfile.R script", log = TRUE)

# Create data output directory
if (!dir.exists("data/processed")) {
  dir.create("data/processed", recursive = TRUE)
}

# Construct monthly commodity real prices indices
source("rsrc/processing/prep_cattle_price_index.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

# Extract sample of mapbiomas 30m-pixels and recover full panel (1985-2019)
source("rsrc/processing/prep_pixel_sample.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

# Calculate area of 30m-pixels
source("rsrc/processing/prep_pixel_area.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

# Categorize 30m-pixels
source("rsrc/processing/prep_pixel_categories.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

# Obtain aboveground biomass for 30m-pixels classified as primary forest
source("rsrc/processing/prep_pixel_biomass.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

# Define municipal-level sample
source("rsrc/processing/prep_muni_sample.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

# Obtain aboveground biomass in primary forest at municipal level
source("rsrc/processing/prep_muni_biomass.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

# Merge clean sources into full municipal-level data set
source("rsrc/processing/prep_muni_data.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

# Prepare data to estimate agricultural emissions factor
source("rsrc/processing/prep_state_emissions.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

# Prepare rasters for mapbiomas variables
source("rsrc/processing/prep_rst_land_use.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

# Prepare generate aggregated sample of interest
source("rsrc/processing/prep_rst_biome.R", encoding = "UTF-8", echo = TRUE)

# Clear environment
rm(list = ls())

# End timer
toc(log = TRUE)
