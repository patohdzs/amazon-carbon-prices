library(sf)
library(tictoc)
library(terra)
library(tidyverse)
library(conflicted)
library(sjlabelled)

# Resolve conflicts
conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)

# Start timer
tic(msg = "merge_agb_sec_veg.R script", log = TRUE)

# Load calibrated gamma
load("data/calibration/gamma_calibration_1043_sites.Rdata")

# Load clean mapbiomas raster
clean_mapbiomas <- rast("data/clean/land_use_cover_2000.tif")

# Load secondary vegetation age data for 2017
sec_veg_age_rst <- rast(
  list.files(
    "data/raw/mapbiomas/secondary_vegetation_age/",
    pattern = "2017",
    full.names = TRUE
  )
)

# Load aboveground biomass raster data for 2017
agb_rst <- map(
  list.files(
    "data/raw/esa/above_ground_biomass/",
    pattern = "_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2017-fv3.0.tif",
    full.names = TRUE
  ),
  rast
)

pq_rst <- rast(
  list.files(
    "data/raw/mapbiomas/pasture_quality/",
    full.names = TRUE
  )
)

# Load percipitation data
precip_rsts <- rast("data/clean/precipitation.tif")

# Load GHI data
ghi_rst <- rast("data/raw/GHI.tif")

# Combine AGB zonal rasters into a single raster
agb_rst <- do.call(merge, agb_rst)

# Crop the the Amazonia subset
sec_veg_age_rst <- sec_veg_age_rst |>
  crop(clean_mapbiomas) |>
  mask(clean_mapbiomas)

# Resample percipitation data into sites
mean_precip_rst <- precip_rsts |> mean()

# Get gamma and site number as rasters
gamma_rst <- calib_df |>
  rasterize(sec_veg_age_rst, field = "site_reg_gamma")

site_rst <- calib_df |>
  rasterize(sec_veg_age_rst, field = "id")

# Set pixels outside the  pasture quality time range to NA
sec_veg_age_rst[sec_veg_age_rst < 2 | sec_veg_age_rst > 17] <- NA

# Take random sample of pixels
n_samples <- 120000000
sec_veg_age <- spatSample(
  sec_veg_age_rst,
  n_samples,
  method = "regular",
  na.rm = TRUE,
  xy = TRUE
)

# Extract values from rasters for the sampled pixels
pq <- extract(pq_rst, sec_veg_age[, c("x", "y")])
agb <- extract(agb_rst, sec_veg_age[, c("x", "y")])
gamma <- extract(gamma_rst, sec_veg_age[, c("x", "y")])
site <- extract(site_rst, sec_veg_age[, c("x", "y")])
mean_precip <- extract(mean_precip_rst, sec_veg_age[, c("x", "y")])
ghi <- extract(ghi_rst, sec_veg_age[, c("x", "y")])

# Combine the extracted values
df <- data.frame(sec_veg_age, pq, agb, gamma, site, mean_precip, ghi) |>
  as_tibble()

# Rename some columns
df <- df |>
  rename(agb = !!names(df)[ncol(df) - 8]) |>
  rename(
    site = id,
    gamma = site_reg_gamma,
    age = sec_veg_age_2017,
    percip = mean,
    ghi = GHI
  ) |>
  select(-matches("^ID"))

# Get last pasture quality
df <- df |>
  mutate(last_year_pasture = 2017 - age) |>
  mutate(last_pq = case_when(
    last_year_pasture == 2000 ~ pasture_quality_2000,
    last_year_pasture == 2001 ~ pasture_quality_2001,
    last_year_pasture == 2002 ~ pasture_quality_2002,
    last_year_pasture == 2003 ~ pasture_quality_2003,
    last_year_pasture == 2004 ~ pasture_quality_2004,
    last_year_pasture == 2005 ~ pasture_quality_2005,
    last_year_pasture == 2006 ~ pasture_quality_2006,
    last_year_pasture == 2007 ~ pasture_quality_2007,
    last_year_pasture == 2008 ~ pasture_quality_2008,
    last_year_pasture == 2009 ~ pasture_quality_2009,
    last_year_pasture == 2010 ~ pasture_quality_2010,
    last_year_pasture == 2011 ~ pasture_quality_2011,
    last_year_pasture == 2012 ~ pasture_quality_2012,
    last_year_pasture == 2013 ~ pasture_quality_2013,
    last_year_pasture == 2014 ~ pasture_quality_2014,
    last_year_pasture == 2015 ~ pasture_quality_2015,
    last_year_pasture == 2016 ~ pasture_quality_2016,
    last_year_pasture == 2017 ~ pasture_quality_2017
  ))

# Convert AGB -> CO2e and get accumulation ratio
df <- df |>
  mutate(co2e = (agb / 2) * (44 / 12)) |>
  mutate(ratio = co2e / gamma)

# Save
save(df, file = "data/calibration/carbon_accumulation.Rdata")
