library(sf)
library(tictoc)
library(tidyverse)

tictoc::tic(msg = "calibrate_gamma_reg.R script", log = TRUE)

# Load variables at the muni level to calibrate gamma
load("data/processed/muni_data.Rdata")

# Convert biomass into CO2e, add column of ones, take logs, and scale
df <- muni_data %>%
  mutate(co2e_ha_2017 = (agb_2017 / 2) * (44 / 12)) %>%
  mutate(
    X1 = 1,
    log_historical_precip = log(historical_precip),
    log_historical_temp = log(historical_temp),
    log_lat = log(lat),
    log_lon = log(lon),
    log_co2e_ha_2017 = log(co2e_ha_2017)
  ) %>%
  mutate(
    log_historical_precip = scale(log_historical_precip),
    log_historical_temp = scale(log_historical_temp),
    log_lat = scale(log_lat),
    log_lon = scale(log_lon)
  ) %>%
  select(
    X1,
    log_historical_precip,
    log_historical_temp,
    log_lat,
    log_lon,
    log_co2e_ha_2017
  )


# Output municipality-level regression data
st_write(df,
  "data/calibration/hmc/gamma_reg_muni_data.geojson",
  driver = "GeoJSON",
  delete_dsn = TRUE
)


# Output site-level regression data
for (n in list(78, 1043)) {
  # Get site boundaries
  id_df <- st_read(sprintf("data/calibration/grid_%d_sites.geojson", n))

  # Project data to site level
  site_level_df <- df %>%
    st_intersection(id_df)

  # Set area units to hectares
  site_level_df$muni_site_area <-
    st_area(site_level_df) %>%
    units::set_units(ha) %>%
    unclass()

  # Write to file
  st_write(
    site_level_df,
    sprintf("data/calibration/hmc/gamma_reg_%d_sites_data.geojson", n),
    driver = "GeoJSON",
    delete_dsn = TRUE
  )
}

# END TIMER
tictoc::toc(log = TRUE)
