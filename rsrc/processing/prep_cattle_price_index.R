# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CONSTRUCT TRIMESTRAL COMMODITY REAL PRICES INDICES
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -


library(tictoc)
library(tidyverse)
library(conflicted)
library(sjlabelled)

# Start timer
tic(msg = "prep_cattle_price_index.R script", log = TRUE)

# Load agricultural commodity prices
load("data/clean/commodity_prices.Rdata")

# Load deflator
load("data/clean/deflator.Rdata")

# Subset deflator in january 2017
aux_01_2017 <- deflator[deflator$date == "2017-01-01", ]$deflator_ipa

# Transform nominal monthly commodity prices to real standardized prices
# Merge data sets, set deflator base year, and transform to real price
# Finally, construct year, month and trimester variables
cattle_price_index <-
  commodity_prices %>%
  left_join(deflator) %>%
  mutate(deflator_ipa = deflator_ipa / aux_01_2017) %>%
  mutate(price_real_mon_cattle = price_cattle / deflator_ipa) %>%
  mutate(
    month = month(date),
    year = year(date)
  ) %>%
  select(date, year, month, starts_with("price_real"), price_cattle)

# Set labels
set_label(cattle_price_index$date) <- "monthly date"
set_label(cattle_price_index$year) <- "calendar year"
set_label(cattle_price_index$month) <- "month indicator"
set_label(cattle_price_index$price_real_mon_cattle) <- "real average monthly price, cattle (R$ per @, constant 01/2017)"

brl_to_usd <- 3.192
cattle_price_index_subset <- 
  cattle_price_index %>%
  dplyr::filter(year >= 1995 & year <= 2017) %>%
  mutate(price_real_mon_cattle=price_real_mon_cattle/brl_to_usd)
  


write_csv(cattle_price_index_subset, file = "data/calibration/seriesPriceCattle_prepared.csv")

# Save data set
save(cattle_price_index, file = "data/processed/cattle_price_index.Rdata")

# Write to csv file
write_csv(cattle_price_index, file = "data/calibration/cattle_price_index.csv")

# End timer
toc(log = TRUE)
