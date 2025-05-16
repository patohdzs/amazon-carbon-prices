# Load necessary libraries
library(tidyverse)



# Load variables at the muni level to calibrate gamma
load("data/calibration/gamma_calibration_1043_sites.Rdata")

# Convert biomass into CO2e, add column of ones, take logs, and scale
df <- calib_1043 %>%
  mutate(
    X1 = 1,
    log_hist_precip = log(hist_precip),
    log_hist_temp = log(hist_temp),
    log_co2e_ha_2017 = log(co2e),
  ) %>%
  filter(!is.na(log_co2e_ha_2017))


# Step 1: Calculate the mean and standard deviation from the first dataset
scaling_params <- df %>%
  summarise(
    mean_log_hist_precip = mean(log_hist_precip, na.rm = TRUE),
    sd_log_hist_precip = sd(log_hist_precip, na.rm = TRUE),
    mean_log_hist_temp = mean(log_hist_temp, na.rm = TRUE),
    sd_log_hist_temp = sd(log_hist_temp, na.rm = TRUE),
    mean_lat = mean(lat, na.rm = TRUE),
    sd_lat = sd(lat, na.rm = TRUE),
    mean_lon = mean(lon, na.rm = TRUE),
    sd_lon = sd(lon, na.rm = TRUE)
  )

# Step 2: Use these parameters to scale the first dataset
df_scaled <- df %>%
  mutate(
    log_hist_precip = (log_hist_precip - scaling_params$mean_log_hist_precip) / scaling_params$sd_log_hist_precip,
    log_hist_temp = (log_hist_temp - scaling_params$mean_log_hist_temp) / scaling_params$sd_log_hist_temp,
    lat = (lat - scaling_params$mean_lat) / scaling_params$sd_lat,
    lon = (lon - scaling_params$mean_lon) / scaling_params$sd_lon,
    latlon = lat * lon
  ) %>%
  select(
    X1,
    log_hist_precip,
    log_hist_temp,
    lat,
    lon,
    latlon,
    log_co2e_ha_2017,
    id_group
  )


load("data/calibration/gamma_calibration_78_sites.Rdata")

df2 <- calib_df %>%
  mutate(
    X1 = 1,
    log_hist_precip = log(hist_precip),
    log_hist_temp = log(hist_temp),
    latlon = lat * lon,
    log_co2e_ha_2017 = log(co2e)
  ) 

# Step 3: Apply the same scaling to the second dataset
df2_scaled <- df2 %>%
  mutate(
    log_hist_precip = (log_hist_precip - scaling_params$mean_log_hist_precip) / scaling_params$sd_log_hist_precip,
    log_hist_temp = (log_hist_temp - scaling_params$mean_log_hist_temp) / scaling_params$sd_log_hist_temp,
    lat = (lat - scaling_params$mean_lat) / scaling_params$sd_lat,
    lon = (lon - scaling_params$mean_lon) / scaling_params$sd_lon,
    latlon = lat * lon,
    id_group = row_number(),
  ) %>%
  select(
    X1,
    log_hist_precip,
    log_hist_temp,
    lat,
    lon,
    latlon,
    log_co2e_ha_2017,
    id_group,
  ) 








# Output municipality-level regression data
st_write(df_scaled,
         "data/calibration/gamma_reg.geojson",
         driver = "GeoJSON",
         delete_dsn = TRUE
)

# Output municipality-level regression data
st_write(df2_scaled,
         "data/calibration/gamma_fit_78.geojson",
         driver = "GeoJSON",
         delete_dsn = TRUE
)




# Convert biomass into CO2e, add column of ones, take logs, and scale
df3 <- calib_1043 %>%
  mutate(
    X1 = 1,
    log_hist_precip = log(hist_precip),
    log_hist_temp = log(hist_temp),
    latlon = lat * lon,
    log_co2e_ha_2017 = log(co2e)
  ) %>%
  mutate(
    log_hist_precip = (log_hist_precip - scaling_params$mean_log_hist_precip) / scaling_params$sd_log_hist_precip,
    log_hist_temp = (log_hist_temp - scaling_params$mean_log_hist_temp) / scaling_params$sd_log_hist_temp,
    lat = (lat - scaling_params$mean_lat) / scaling_params$sd_lat,
    lon = (lon - scaling_params$mean_lon) / scaling_params$sd_lon,
    latlon = lat * lon
  ) %>%
  select(
    X1,
    log_hist_precip,
    log_hist_temp,
    lat,
    lon,
    latlon,
    log_co2e_ha_2017,
    id_group
  ) 



# Output municipality-level regression data
st_write(df3,
         "data/calibration/gamma_fit_1043.geojson",
         driver = "GeoJSON",
         delete_dsn = TRUE
)



