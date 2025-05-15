# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: PARAMETERS CALIBRATION (78 SITES MODEL)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

library(sf)
library(units)
library(terra)
library(readxl)
library(tidyverse)
library(tictoc)
library(conflicted)

conflicts_prefer(dplyr::filter)

# Start timer
tic(msg = "calibrate_78_sites_model.R script", log = TRUE)

# Load spatial municipal sample
load("data/processed/spatial_muni_sample.Rdata")

# Load municipal data
load("data/processed/muni_data.Rdata")

# Load state emissions data
load("data/processed/state_emissions.Rdata")

# Load cattle price index
load("data/processed/cattle_price_index.Rdata")

# Load raster data (amazon biome share, pixel areas, and land uses)
rsts <- rast(
  list.files(
    "data/processed/",
    pattern = ".tif",
    full.names = TRUE
  )
)

# Transform shares to areas
rsts$area_amazon_biome <-
  rsts$share_amazon_biome * rsts$pixel_area_ha

rsts$area_forest_1995 <-
  rsts$share_forest_1995 * rsts$pixel_area_ha

rsts$area_agricultural_use_1995 <-
  rsts$share_agricultural_use_1995 * rsts$pixel_area_ha

rsts$area_other_1995 <-
  rsts$share_other_1995 * rsts$pixel_area_ha

rsts$area_forest_2017 <-
  rsts$share_forest_2017 * rsts$pixel_area_ha

# rsts$area_pasture_quality_1_2017 <-
#   rsts$share_pasture_quality_1_2017 * rsts$pixel_area_ha

# rsts$area_pasture_quality_2_2017 <-
#   rsts$share_pasture_quality_2_2017 * rsts$pixel_area_ha

# rsts$area_pasture_quality_3_2017 <-
#   rsts$share_pasture_quality_3_2017 * rsts$pixel_area_ha

rsts$area_agricultural_use_2017 <-
  rsts$share_agricultural_use_2017 * rsts$pixel_area_ha

rsts$area_other_2017 <-
  rsts$share_other_2017 * rsts$pixel_area_ha

rsts$area_forest_2008 <-
  rsts$share_forest_2008 * rsts$pixel_area_ha

rsts$area_agricultural_use_2008 <-
  rsts$share_agricultural_use_2008 * rsts$pixel_area_ha

rsts$area_other_2008 <-
  rsts$share_other_2008 * rsts$pixel_area_ha

# Aggregate into larger sites
rsts <- aggregate(
  rsts,
  fact = 4,
  fun = sum,
  na.rm = TRUE
)

# Extract variables as polygons and project data
calib_df <- as.polygons(rsts, dissolve = FALSE) %>%
  st_as_sf() %>%
  st_transform(5880)

# Remove sites with less than 3% overlap with the amazon biome
calib_df <- calib_df %>%
  filter(area_amazon_biome / pixel_area_ha >= 0.03)

# Add id variable
calib_df$id <- seq_len(nrow(calib_df))

# Change names to z and add zbar
calib_df <- calib_df %>%
  rename(site_area_ha = pixel_area_ha) %>%
  mutate(
    z_1995 = area_agricultural_use_1995,
    z_2008 = area_agricultural_use_2008,
    z_2017 = area_agricultural_use_2017,
    zbar_1995 = area_forest_1995 + z_1995,
    zbar_2008 = area_forest_2008 + z_2008,
    zbar_2017 = area_forest_2017 + z_2017,
  )

# # Compute bad pasture quality shares
# calib_df <- calib_df %>%
#   mutate(
#     area_pasture_2017 = area_pasture_quality_1_2017 +
#       area_pasture_quality_2_2017 +
#       area_pasture_quality_3_2017
#   ) %>%
#   mutate(
#     share_low_pq = area_pasture_quality_1_2017 / area_pasture_2017,
#   ) %>%
#   select(id, contains("area"), contains("z"), contains("pq"))

# Convert biomass to CO2 equivalent
muni_data <- muni_data %>%
  mutate(co2e_ha_2017 = (agb_2017 / 2) * (44 / 12))

# Regress log-gamma on geographic covariates
gamma_reg <- lm(
  formula = log(co2e_ha_2017) ~
    log(historical_precip) +
    log(historical_temp) +
    log(lat) +
    log(lon),
  data = muni_data,
  na.action = na.exclude
)

# Compute fitted values and transform back into levels
muni_data <- muni_data %>%
  mutate(co2e_ha_2017_fitted = exp(predict(gamma_reg, .)))

# Match municipalities with sites
site_level_gamma <- calib_df %>%
  st_intersection(muni_data) %>%
  select(
    id,
    muni_code,
    muni_area,
    co2e_ha_2017_fitted
  )

# Set area units
site_level_gamma$muni_site_area <-
  st_area(site_level_gamma) %>%
  set_units(ha) %>%
  unclass()

# Drop spatial feature
site_level_gamma <- site_level_gamma %>%
  st_drop_geometry()

# Average carbon density on primary forest areas by site
site_level_gamma <- site_level_gamma %>%
  group_by(id) %>%
  summarise(
    gamma = weighted.mean(
      co2e_ha_2017_fitted,
      w = muni_site_area,
      na.rm = TRUE
    )
  )

# Add site gamma to calibration variables
calib_df <- calib_df %>% left_join(site_level_gamma)

# Clean environment
rm(site_level_gamma)

# Estimate of alpha same as in the global model
calib_df <- calib_df %>%
  mutate(alpha = 1 - (1 - 0.99)^(1 / 100))

# Set discount rate (delta) of 2%
calib_df <- calib_df %>%
  mutate(delta = 0.02)

# Calculate average net emission factor
#   from agricultural use across years and states
avg_net_emission_factor <- state_emissions %>%
  ungroup() %>%
  summarise(
    nef = sum(net_emissions_co2e) / sum(agriculturalUse_area)
  ) %>%
  pull(nef)

# Estimate kappa
calib_df <- calib_df %>%
  mutate(kappa = avg_net_emission_factor)

# Clean environment
rm(avg_net_emission_factor)

# Zeta is calibrated such that the marginal cost of changing land use
#   (zeta*forest_to_pasture_transition_area) matches the forest to
#   pasture transition cost estimated by Araujo, Costa and Sant'Anna (2022),
#   reported on column 4 of table 4 (on the right).

# We transform to dollars using an FX rate of 4.14 (December 2019)
aux_transition_cost <- 1614.54 / 4.14

# The forest area in 2008 represented 72% of the Legal Amazon area,
#   which covers 501,506,775 ha, so the transition in hectares is
#   0.065*0.72*501,506,775 = 23,470,517,
#   resulting in an annual average of 2,347,052 ha.
aux_transition_area <- (0.065 * 0.72 * 501506775) / (2017 - 2008 + 1)
zeta <- aux_transition_cost / aux_transition_area

# Alternative value based on https://shorturl.at/jgp9i
zeta_alt <- 483 / aux_transition_area

# Estimate of zeta same as in the global model
calib_df <- calib_df %>%
  mutate(zeta = zeta, zeta_alt = zeta_alt)

# x_2017 is estimated as the carbon stock stored in forest areas,
#   assuming that all forests are primary.
calib_df <- calib_df %>%
  mutate(
    x_1995 = gamma * (zbar_1995 - z_1995),
    x_2008 = gamma * (zbar_2008 - z_2008),
    x_2017 = gamma * (zbar_2017 - z_2017)
  )

# Remove outliers from municipal data
muni_data <- muni_data[-c(142, 106, 112), ]

# Filter out observations missing distance to capital
muni_data <- muni_data %>%
  filter(!is.na(distance))

# Theta regression
theta_reg <- muni_data %>%
  filter(slaughter_value_per_ha_2017 > 0) %>%
  lm(
    formula = log(slaughter_value_per_ha_2017) ~
      lat +
      I(lat^2) +
      historical_temp +
      I(historical_temp^2) +
      historical_precip +
      distance +
      log(farm_gate_price_2017),
    na.action = na.exclude,
    weights = pasture_area_2017
  )

# Extract fitted values
muni_data <- muni_data %>%
  mutate(slaughter_value_per_ha_fitted = exp(predict(theta_reg, .)))

# Match municipalities with sites
site_level_theta <- calib_df %>%
  st_intersection(muni_data) %>%
  select(
    id,
    muni_code,
    muni_area,
    pasture_area_2017,
    slaughter_value_per_ha_fitted
  )

# Calculate municipality areas inside each site
site_level_theta$muni_site_area <-
  st_area(site_level_theta) %>%
  set_units(ha) %>%
  unclass()

# Drop spatial feature
site_level_theta <- site_level_theta %>%
  st_drop_geometry()

# Extract average cattle price index for 2017
# BRL to USD (Commercial ER - selling - average - annual - 2017 - ipeadata))
brl_to_usd <- 3.192

mean_pa_2017 <-
  cattle_price_index %>%
  filter(year == 2017) %>%
  group_by(year) %>%
  summarise(mean_price_2017 = mean(price_real_mon_cattle) / brl_to_usd) %>%
  pull(mean_price_2017)

calib_df <- calib_df %>%
  mutate(mean_pa_2017 = mean_pa_2017)

# Calculate theta and pasture area by site
# (for each muni adjust the value by the share of the muni area inside the site)
site_level_theta <- site_level_theta %>%
  group_by(id) %>%
  summarise(
    theta = weighted.mean(
      slaughter_value_per_ha_fitted / mean_pa_2017,
      w = muni_site_area,
      na.rm = TRUE
    ),
    pasture_area_2017 = sum(pasture_area_2017 * (muni_site_area / muni_area), na.rm = TRUE),
  )

# Add thetas and pasture_area to calibration data
calib_df <- calib_df %>% left_join(site_level_theta)

# Clean environment
rm(site_level_theta)

# Filter out missing vlaues
calib_df <- calib_df %>%
  filter(!is.na(theta))

# Add site ID's
calib_df <- calib_df %>%
  mutate(id = seq_len(nrow(calib_df)))

# Save site boundaries
calib_df %>%
  select(id) %>%
  st_write(
    "data/calibration/grid_78_sites.geojson",
    driver = "GeoJSON",
    delete_dsn = TRUE
  )

# Save calibration data
save(calib_df, file = "data/calibration/calibration_78_sites.Rdata")

# Save calibration data CSV
calib_df %>%
  st_drop_geometry() %>%
  write_csv(file = "data/calibration/calibration_78_sites.csv")

# End timer
toc(log = TRUE)
# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: PARAMETERS CALIBRATION (78 SITES MODEL)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

library(sf)
library(units)
library(terra)
library(readxl)
library(tidyverse)
library(tictoc)
library(conflicted)

conflicts_prefer(dplyr::filter)

# Start timer
tic(msg = "calibrate_78_sites_model.R script", log = TRUE)

# Load spatial municipal sample
load("data/processed/spatial_muni_sample.Rdata")

# Load municipal data
load("data/processed/muni_data.Rdata")

# Load state emissions data
load("data/processed/state_emissions.Rdata")

# Load cattle price index
load("data/processed/cattle_price_index.Rdata")

# Load raster data (amazon biome share, pixel areas, and land uses)
rsts <- rast(
  list.files(
    "data/processed/",
    pattern = ".tif",
    full.names = TRUE
  )
)

# Transform shares to areas
rsts$area_amazon_biome <-
  rsts$share_amazon_biome * rsts$pixel_area_ha

rsts$area_forest_1995 <-
  rsts$share_forest_1995 * rsts$pixel_area_ha

rsts$area_agricultural_use_1995 <-
  rsts$share_agricultural_use_1995 * rsts$pixel_area_ha

rsts$area_other_1995 <-
  rsts$share_other_1995 * rsts$pixel_area_ha

rsts$area_forest_2017 <-
  rsts$share_forest_2017 * rsts$pixel_area_ha

# rsts$area_pasture_quality_1_2017 <-
#   rsts$share_pasture_quality_1_2017 * rsts$pixel_area_ha

# rsts$area_pasture_quality_2_2017 <-
#   rsts$share_pasture_quality_2_2017 * rsts$pixel_area_ha

# rsts$area_pasture_quality_3_2017 <-
#   rsts$share_pasture_quality_3_2017 * rsts$pixel_area_ha

rsts$area_agricultural_use_2017 <-
  rsts$share_agricultural_use_2017 * rsts$pixel_area_ha

rsts$area_other_2017 <-
  rsts$share_other_2017 * rsts$pixel_area_ha

rsts$area_forest_2008 <-
  rsts$share_forest_2008 * rsts$pixel_area_ha

rsts$area_agricultural_use_2008 <-
  rsts$share_agricultural_use_2008 * rsts$pixel_area_ha

rsts$area_other_2008 <-
  rsts$share_other_2008 * rsts$pixel_area_ha

# Aggregate into larger sites
rsts <- aggregate(
  rsts,
  fact = 4,
  fun = sum,
  na.rm = TRUE
)

# Extract variables as polygons and project data
calib_df <- as.polygons(rsts, dissolve = FALSE) %>%
  st_as_sf() %>%
  st_transform(5880)

# Remove sites with less than 3% overlap with the amazon biome
calib_df <- calib_df %>%
  filter(area_amazon_biome / pixel_area_ha >= 0.03)

# Add id variable
calib_df$id <- seq_len(nrow(calib_df))

# Change names to z and add zbar
calib_df <- calib_df %>%
  rename(site_area_ha = pixel_area_ha) %>%
  mutate(
    z_1995 = area_agricultural_use_1995,
    z_2008 = area_agricultural_use_2008,
    z_2017 = area_agricultural_use_2017,
    zbar_1995 = area_forest_1995 + z_1995,
    zbar_2008 = area_forest_2008 + z_2008,
    zbar_2017 = area_forest_2017 + z_2017,
  )

# # Compute bad pasture quality shares
# calib_df <- calib_df %>%
#   mutate(
#     area_pasture_2017 = area_pasture_quality_1_2017 +
#       area_pasture_quality_2_2017 +
#       area_pasture_quality_3_2017
#   ) %>%
#   mutate(
#     share_low_pq = area_pasture_quality_1_2017 / area_pasture_2017,
#   ) %>%
#   select(id, contains("area"), contains("z"), contains("pq"))

# Convert biomass to CO2 equivalent
muni_data <- muni_data %>%
  mutate(co2e_ha_2017 = (agb_2017 / 2) * (44 / 12))

# Regress log-gamma on geographic covariates
gamma_reg <- lm(
  formula = log(co2e_ha_2017) ~
    log(historical_precip) +
    log(historical_temp) +
    log(lat) +
    log(lon),
  data = muni_data,
  na.action = na.exclude
)

# Compute fitted values and transform back into levels
muni_data <- muni_data %>%
  mutate(co2e_ha_2017_fitted = exp(predict(gamma_reg, .)))

# Match municipalities with sites
site_level_gamma <- calib_df %>%
  st_intersection(muni_data) %>%
  select(
    id,
    muni_code,
    muni_area,
    co2e_ha_2017_fitted
  )

# Set area units
site_level_gamma$muni_site_area <-
  st_area(site_level_gamma) %>%
  set_units(ha) %>%
  unclass()

# Drop spatial feature
site_level_gamma <- site_level_gamma %>%
  st_drop_geometry()

# Average carbon density on primary forest areas by site
site_level_gamma <- site_level_gamma %>%
  group_by(id) %>%
  summarise(
    gamma = weighted.mean(
      co2e_ha_2017_fitted,
      w = muni_site_area,
      na.rm = TRUE
    )
  )

# Add site gamma to calibration variables
calib_df <- calib_df %>% left_join(site_level_gamma)

# Clean environment
rm(site_level_gamma)

# Estimate of alpha same as in the global model
calib_df <- calib_df %>%
  mutate(alpha = 1 - (1 - 0.99)^(1 / 100))

# Set discount rate (delta) of 2%
calib_df <- calib_df %>%
  mutate(delta = 0.02)

# Calculate average net emission factor
#   from agricultural use across years and states
avg_net_emission_factor <- state_emissions %>%
  ungroup() %>%
  summarise(
    nef = sum(net_emissions_co2e) / sum(agriculturalUse_area)
  ) %>%
  pull(nef)

# Estimate kappa
calib_df <- calib_df %>%
  mutate(kappa = avg_net_emission_factor)

# Clean environment
rm(avg_net_emission_factor)

# Zeta is calibrated such that the marginal cost of changing land use
#   (zeta*forest_to_pasture_transition_area) matches the forest to
#   pasture transition cost estimated by Araujo, Costa and Sant'Anna (2022),
#   reported on column 4 of table 4 (on the right).

# We transform to dollars using an FX rate of 4.14 (December 2019)
aux_transition_cost <- 1614.54 / 4.14

# The forest area in 2008 represented 72% of the Legal Amazon area,
#   which covers 501,506,775 ha, so the transition in hectares is
#   0.065*0.72*501,506,775 = 23,470,517,
#   resulting in an annual average of 2,347,052 ha.
aux_transition_area <- (0.065 * 0.72 * 501506775) / (2017 - 2008 + 1)
zeta <- aux_transition_cost / aux_transition_area

# Alternative value based on https://shorturl.at/jgp9i
zeta_alt <- 483 / aux_transition_area

# Estimate of zeta same as in the global model
calib_df <- calib_df %>%
  mutate(zeta = zeta, zeta_alt = zeta_alt)

# x_2017 is estimated as the carbon stock stored in forest areas,
#   assuming that all forests are primary.
calib_df <- calib_df %>%
  mutate(
    x_1995 = gamma * (zbar_1995 - z_1995),
    x_2008 = gamma * (zbar_2008 - z_2008),
    x_2017 = gamma * (zbar_2017 - z_2017)
  )

# Remove outliers from municipal data
muni_data <- muni_data[-c(142, 106, 112), ]

# Filter out observations missing distance to capital
muni_data <- muni_data %>%
  filter(!is.na(distance))

# Theta regression
theta_reg <- muni_data %>%
  filter(slaughter_value_per_ha_2017 > 0) %>%
  lm(
    formula = log(slaughter_value_per_ha_2017) ~
      lat +
      I(lat^2) +
      historical_temp +
      I(historical_temp^2) +
      historical_precip +
      distance +
      log(farm_gate_price_2017),
    na.action = na.exclude,
    weights = pasture_area_2017
  )

# Extract fitted values
muni_data <- muni_data %>%
  mutate(slaughter_value_per_ha_fitted = exp(predict(theta_reg, .)))

# Match municipalities with sites
site_level_theta <- calib_df %>%
  st_intersection(muni_data) %>%
  select(
    id,
    muni_code,
    muni_area,
    pasture_area_2017,
    slaughter_value_per_ha_fitted
  )

# Calculate municipality areas inside each site
site_level_theta$muni_site_area <-
  st_area(site_level_theta) %>%
  set_units(ha) %>%
  unclass()

# Drop spatial feature
site_level_theta <- site_level_theta %>%
  st_drop_geometry()

# Extract average cattle price index for 2017
# BRL to USD (Commercial ER - selling - average - annual - 2017 - ipeadata))
brl_to_usd <- 3.192

mean_pa_2017 <-
  cattle_price_index %>%
  filter(year == 2017) %>%
  group_by(year) %>%
  summarise(mean_price_2017 = mean(price_real_mon_cattle) / brl_to_usd) %>%
  pull(mean_price_2017)

calib_df <- calib_df %>%
  mutate(mean_pa_2017 = mean_pa_2017)

# Calculate theta and pasture area by site
# (for each muni adjust the value by the share of the muni area inside the site)
site_level_theta <- site_level_theta %>%
  group_by(id) %>%
  summarise(
    theta = weighted.mean(
      slaughter_value_per_ha_fitted / mean_pa_2017,
      w = muni_site_area,
      na.rm = TRUE
    ),
    pasture_area_2017 = sum(pasture_area_2017 * (muni_site_area / muni_area), na.rm = TRUE),
  )

# Add thetas and pasture_area to calibration data
calib_df <- calib_df %>% left_join(site_level_theta)

# Clean environment
rm(site_level_theta)

# Filter out missing vlaues
calib_df <- calib_df %>%
  filter(!is.na(theta))

# Add site ID's
calib_df <- calib_df %>%
  mutate(id = seq_len(nrow(calib_df)))

# Save site boundaries
calib_df %>%
  select(id) %>%
  st_write(
    "data/calibration/grid_78_sites.geojson",
    driver = "GeoJSON",
    delete_dsn = TRUE
  )

# Save calibration data
save(calib_df, file = "data/calibration/calibration_78_sites.Rdata")

# Save calibration data CSV
calib_df %>%
  st_drop_geometry() %>%
  write_csv(file = "data/calibration/calibration_78_sites.csv")

# End timer
toc(log = TRUE)
