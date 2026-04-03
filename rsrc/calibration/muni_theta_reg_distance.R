library(stargazer)

library(terra)

library(sf)
library(conflicted)
library(units)
library(dplyr)
conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)


# Load municipal data
load("data/processed/muni_data.Rdata")

load("data/processed/cattle_price_index.Rdata")

# Load calibration data set
load("data/calibration/calibration_78_sites.Rdata")

# load basin geo data
map_basin <- st_read(
  dsn = "data/raw/mapbiomas/basin",
  layer = "BASIN_LEVEL_2_PNRH"
)

calib_df <- calib_df %>%
  st_transform(st_crs(map_basin))

muni_data <- muni_data %>%
  st_transform(st_crs(map_basin))

centroids <- muni_data %>%
  st_centroid() %>%
  st_coordinates()

muni_data$lon <- centroids[, "X"]
muni_data$lat <- centroids[, "Y"]






map_basin_tf <- map_basin %>%
  select(FEATURE_ID)

muni_with_basin_full <- st_join(muni_data, map_basin_tf, join = st_nearest_feature)

muni_to_reassign <- muni_with_basin_full %>%
  filter(FEATURE_ID %in% c(11, 100, 112, 154, 256)) %>%
  select(-FEATURE_ID)


map_basin_tf_reassign <- map_basin %>%
  select(FEATURE_ID) %>%
  filter(!FEATURE_ID %in% c(11, 100, 107, 112, 154, 256))

muni_to_reassign <- st_join(muni_to_reassign, map_basin_tf_reassign, join = st_nearest_feature)

muni_with_basin <- bind_rows(
  muni_with_basin_full %>% filter(!FEATURE_ID %in% c(11, 100, 112, 154, 256)),
  muni_to_reassign
)




muni_data <- muni_with_basin


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
summary(theta_reg)






# Load necessary library
library(lme4)
theta_reg_re <- muni_data %>%
  filter(slaughter_value_per_ha_2017 > 0) %>%
  lmer(
    formula = log(slaughter_value_per_ha_2017) ~
      lat +
      I(lat^2) +
      historical_temp +
      I(historical_temp^2) +
      historical_precip +
      distance +
      log(farm_gate_price_2017) +
      (1 | FEATURE_ID), # Random intercept for FEATURE_ID
    na.action = na.exclude,
    weights = pasture_area_2017 / mean(pasture_area_2017)
  )

# Show the summary of the random effects model
summary(theta_reg_re)

random_effects <- ranef(theta_reg_re)$FEATURE_ID %>%
  mutate(obs_id = row_number())



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
    theta_iid = weighted.mean(
      slaughter_value_per_ha_fitted / mean_pa_2017,
      w = muni_site_area,
      na.rm = TRUE
    ),
    pasture_area_2017_i = sum(pasture_area_2017 * (muni_site_area / muni_area), na.rm = TRUE),
  )

# Add thetas and pasture_area to calibration data
calib_df <- calib_df %>%
  inner_join(site_level_theta, by = "id")




save(calib_df, file = "data/calibration/theta_calibration_78_sites.Rdata")


muni_data <- muni_data %>%
  mutate(group_id = as.numeric(factor(FEATURE_ID)), muni_id = muni_code)


df_fit <- muni_data %>%
  select(lat, historical_temp, historical_precip, distance, farm_gate_price_2017, group_id, muni_id) %>%
  mutate(X1 = 1, sq_lat = lat^2, sq_temp = historical_temp^2, log_gate_price = log(farm_gate_price_2017)) %>%
  select(
    X1,
    lat,
    sq_lat,
    historical_temp,
    sq_temp,
    historical_precip,
    distance,
    log_gate_price,
    group_id,
    muni_id
  )

df_reg <- muni_data %>%
  filter(slaughter_value_per_ha_2017 > 0) %>%
  filter(!is.na(slaughter_value_per_ha_2017)) %>%
  select(pasture_area_2017, slaughter_value_per_ha_2017, lat, historical_temp, historical_precip, distance, farm_gate_price_2017, group_id) %>%
  mutate(X1 = 1, sq_lat = lat^2, sq_temp = historical_temp^2, log_gate_price = log(farm_gate_price_2017), log_slaughter = log(slaughter_value_per_ha_2017), weights = pasture_area_2017) %>%
  select(
    X1,
    lat,
    sq_lat,
    historical_temp,
    sq_temp,
    historical_precip,
    distance,
    log_gate_price,
    log_slaughter,
    weights,
    group_id
  )


scaling_params_fit <- df_reg %>%
  st_drop_geometry() %>%
  summarise(
    mean_lat = mean(lat, na.rm = TRUE),
    sd_lat = sd(lat, na.rm = TRUE),
    mean_sq_lat = mean(sq_lat, na.rm = TRUE),
    sd_sq_lat = sd(sq_lat, na.rm = TRUE),
    mean_historical_temp = mean(historical_temp, na.rm = TRUE),
    sd_historical_temp = sd(historical_temp, na.rm = TRUE),
    mean_sq_temp = mean(sq_temp, na.rm = TRUE),
    sd_sq_temp = sd(sq_temp, na.rm = TRUE),
    mean_historical_precip = mean(historical_precip, na.rm = TRUE),
    sd_historical_precip = sd(historical_precip, na.rm = TRUE),
    mean_distance = mean(distance, na.rm = TRUE),
    sd_distance = sd(distance, na.rm = TRUE),
    mean_log_gate_price = mean(log_gate_price, na.rm = TRUE),
    sd_log_gate_price = sd(log_gate_price, na.rm = TRUE)
  )


df_reg_scaled <- df_reg %>%
  mutate(
    lat = (lat - scaling_params_fit$mean_lat) / scaling_params_fit$sd_lat,
    sq_lat = (sq_lat - scaling_params_fit$mean_sq_lat) / scaling_params_fit$sd_sq_lat,
    historical_temp = (historical_temp - scaling_params_fit$mean_historical_temp) / scaling_params_fit$sd_historical_temp,
    sq_temp = (sq_temp - scaling_params_fit$mean_sq_temp) / scaling_params_fit$sd_sq_temp,
    historical_precip = (historical_precip - scaling_params_fit$mean_historical_precip) / scaling_params_fit$sd_historical_precip,
    distance = (distance - scaling_params_fit$mean_distance) / scaling_params_fit$sd_distance,
    log_gate_price = (log_gate_price - scaling_params_fit$mean_log_gate_price) / scaling_params_fit$sd_log_gate_price
  ) %>%
  select(X1, lat, sq_lat, historical_temp, sq_temp, historical_precip, distance, log_gate_price, log_slaughter, weights, group_id)

df_fit_scaled <- df_fit %>%
  mutate(
    lat = (lat - scaling_params_fit$mean_lat) / scaling_params_fit$sd_lat,
    sq_lat = (sq_lat - scaling_params_fit$mean_sq_lat) / scaling_params_fit$sd_sq_lat,
    historical_temp = (historical_temp - scaling_params_fit$mean_historical_temp) / scaling_params_fit$sd_historical_temp,
    sq_temp = (sq_temp - scaling_params_fit$mean_sq_temp) / scaling_params_fit$sd_sq_temp,
    historical_precip = (historical_precip - scaling_params_fit$mean_historical_precip) / scaling_params_fit$sd_historical_precip,
    distance = (distance - scaling_params_fit$mean_distance) / scaling_params_fit$sd_distance,
    log_gate_price = (log_gate_price - scaling_params_fit$mean_log_gate_price) / scaling_params_fit$sd_log_gate_price
  ) %>%
  select(X1, lat, sq_lat, historical_temp, sq_temp, historical_precip, distance, log_gate_price, group_id, muni_id)






st_write(df_reg_scaled,
  "data/calibration/theta_reg.geojson",
  driver = "GeoJSON",
  delete_dsn = TRUE
)




id_sfdata <- calib_df %>%
  select(id)

site_theta <- df_fit_scaled %>%
  st_intersection(id_sfdata)

site_theta$muni_site_area <-
  st_area(site_theta) %>%
  units::set_units(ha) %>%
  unclass()

st_write(
  site_theta,
  sprintf("data/calibration/theta_fit_78.geojson", n),
  driver = "GeoJSON",
  delete_dsn = TRUE
)




load("data/calibration/gamma_calibration_1043_sites.Rdata")

calib_1043 <- calib_1043 %>%
  st_transform(st_crs(muni_data))

id_sfdata <- calib_1043 %>%
  select(id)

site_theta <- df_fit_scaled %>%
  st_intersection(id_sfdata)

site_theta$muni_site_area <-
  st_area(site_theta) %>%
  units::set_units(ha) %>%
  unclass()

st_write(
  site_theta,
  sprintf("data/calibration/theta_fit_1043.geojson", n),
  driver = "GeoJSON",
  delete_dsn = TRUE
)
