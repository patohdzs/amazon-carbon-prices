# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: MASTERFILE SCRIPT TO SOURCE ALL RAW2CLEAN SCRIPTS
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

library(tictoc)

# START TIMER
tic(msg = "_masterfile_raw2clean.R script", log = TRUE)

# Create folder
if (!dir.exists("data/clean")) {
  dir.create("data/clean", recursive = TRUE)
}
if (!dir.exists("data/calibration")) {
  dir.create("data/calibration", recursive = TRUE)
}
if (!dir.exists("data/calibration/hmc")) {
  dir.create("data/calibration/hmc", recursive = TRUE)
}

# SOURCE
# TREAT RAW DATA - AMAZON BIOME BOUNDARY (IBGE - 2019)
source("rsrc/cleaning/clean_amazon_biome_boundary.R", encoding = "UTF-8", echo = TRUE)

# clear environment
rm(list = ls())

# CLEAN RAW DATA SHAPEFILE OF MUNI DIVISION (2015)
source("rsrc/cleaning/clean_muni_boundaries.R", encoding = "UTF-8", echo = TRUE)

# clear environment
rm(list = ls())

# CLEAN RAW AGRICULTURAL COMMODITY PRICES FROM SEAB-PR
source("rsrc/cleaning/clean_commodity_prices.R", encoding = "UTF-8", echo = TRUE)

# clear environment
rm(list = ls())

# CLEAN RAW DEFLATOR IPA-DI - FGV
source("rsrc/cleaning/clean_deflator.R", encoding = "UTF-8", echo = TRUE)

# clear environment
rm(list = ls())

# CLEAN RAW CATTLE SOLD - AGRICULTURAL CENSUS 2017 (IBGE)
source("rsrc/cleaning/clean_cattle_sold.R", encoding = "UTF-8", echo = TRUE)

# clear environment
rm(list = ls())

# CLEAN RAW AGRICULTURAL USE AREA - AGRICULTURAL CENSUS 2017 (IBGE)
source("rsrc/cleaning/clean_agricultural_use_area.R", encoding = "UTF-8", echo = TRUE)

# clear environment
rm(list = ls())

# clear environment
rm(list = ls())

# clear environment
rm(list = ls())

# CLEAN RAW SEEG AGRICULTURAL EMISSION DATA - LEGAL AMAZON STATES 1990-2019
source("rsrc/cleaning/clean_emissions.R", encoding = "UTF-8", echo = TRUE)

# clear environment
rm(list = ls())

# TREAT RAW DATA EMISSIONS AND GDP BY COUNTRY (WORLD BANK)
source("rsrc/cleaning/clean_emission_kuznets.R", encoding = "UTF-8", echo = TRUE)

# clear environment
rm(list = ls())

# TREAT RAW DATA CARBON PRICE (WORLD BANK)
source("rsrc/cleaning/clean_carbon_price.R", encoding = "UTF-8", echo = TRUE)

# clear environment
rm(list = ls())

# CLEAN RAW BRAZILIAN AMAZON LAND USE DATA TO PROVIDE MUNICIPALITY-BY-YEAR PANEL (MAPBIOMAS COLLECTION 5)
source("rsrc/cleaning/clean_muni_land_use_cover.R", encoding = "UTF-8", echo = TRUE)

# clear environment
rm(list = ls())

# UNIFY MONTHLY TIF FILES - HISTORICAL TEMPERATURE (WORLD CLIM)
source("rsrc/cleaning/clean_temperature.R", encoding = "UTF-8", echo = TRUE)

# clear environment
rm(list = ls())

# UNIFY MONTHLY TIF FILES - HISTORICAL PRECIPITATION (WORLD CLIM)
source("rsrc/cleaning/clean_precipitation.R", encoding = "UTF-8", echo = TRUE)

# clear environment
rm(list = ls())

# RECLASSIFY PIXELS OUTSIDE BIOME BOUNDARY - LAND USE AND COVER (MAPBIOMAS - COL5)
source("rsrc/cleaning/clean_land_use_cover.R", encoding = "UTF-8", echo = TRUE)

# clear environment
rm(list = ls())

# DOWNLOAD ABOVERGROUND BIOMASS/CARBON DATA (ESA BIOMASS - 2010, 2017, 2018)
# DOWNLOAD PROCESS IS OPTIONAL GIVEN THAT THE DATA IS PROVIDED
# source("rsrc/cleaning/abovegroundBiomassESA_download.R", encoding = "UTF-8", echo = TRUE)

# END TIMER
toc(log = TRUE)
