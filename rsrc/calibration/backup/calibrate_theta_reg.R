library(sf)
library(units)
library(tictoc)
library(readxl)
library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::filter)

tic(msg = "calibrate_theta_reg.R script", log = TRUE)

# Load variables at the muni level to calibrate theta
load("data/processed/muni_data.Rdata")

# Drop outliers and separate geometry from df
geo <- st_geometry(muni_data)[-c(142, 106, 112)]
muni_data <- as.data.frame(muni_data)[-c(142, 106, 112), ]

# Re-attatch geometry
muni_data <- muni_data %>%
  st_sf(geometry = geo)

# Filter out NaN distances
muni_data <- muni_data %>%
  filter(!is.na(distance))

# Add vector of ones, add polynomial terms, and scale
df <- muni_data %>%
  mutate(
    X1 = 1,
    historical_temp_sq = historical_temp^2,
    lat_sq = lat^2,
    log_slaughter_value_per_ha = log(slaughter_value_per_ha_2017),
    log_farm_gate_price = log(farm_gate_price_2017),
    weights = pasture_area_2017
  ) %>%
  mutate(across(
    c(
      historical_precip,
      historical_temp,
      historical_temp_sq,
      lat,
      lat_sq,
      log_farm_gate_price,
      distance
    ),
    scale
  )) %>%
  select(
    X1,
    historical_precip,
    historical_temp,
    historical_temp_sq,
    lat,
    lat_sq,
    log_farm_gate_price,
    distance,
    log_slaughter_value_per_ha,
    weights
  )

# Output municipality-level regression data
out_file <- "data/calibration/hmc/theta_reg_muni_data.geojson"
st_write(df, out_file, driver = "GeoJSON", delete_dsn = TRUE)

# Output site-level regression data
for (n in list(78, 1043)) {
  id_df <- st_read(sprintf("data/calibration/grid_%d_sites.geojson", n))

  # Project data to site level
  site_level_df <- df %>%
    st_intersection(id_df)

  # Set area units to hectares
  site_level_df$muni_site_area <-
    st_area(site_level_df) %>%
    set_units(ha) %>%
    unclass()

  # Write to file
  st_write(
    site_level_df,
    sprintf("data/calibration/hmc/theta_reg_%d_sites_data.geojson", n),
    driver = "GeoJSON",
    delete_dsn = TRUE
  )
}

# End timer
toc(log = TRUE)
