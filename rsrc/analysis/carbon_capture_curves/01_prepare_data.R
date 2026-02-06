# ============================================================================
# Script: 01_prepare_data.R
# Purpose: Prepare and align raster data for AGB-secondary vegetation analysis
# ============================================================================

# Load libraries ----------------------------------------------------------
library(sf)
library(stars)
library(tictoc)
library(terra)
library(tidyverse)
library(conflicted)

# Resolve conflicts
conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)

# Start timer
tic(msg = "Data preparation script", log = TRUE)

# Load raster data --------------------------------------------------------

# Secondary vegetation age for 2017
sec_veg_age_rst <- rast(
  list.files(
    "data/raw/mapbiomas/secondary_vegetation_age/",
    pattern = "2017",
    full.names = TRUE
  )
)

# Aboveground biomass rasters
agb_rasters <- map(
  list.files(
    "data/raw/esa/above_ground_biomass/",
    pattern = "_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2017-fv3.0.tif",
    full.names = TRUE
  ),
  rast
)

# Pasture quality for 2014
pq_rst <- rast(
  list.files(
    "data/raw/mapbiomas/pasture_quality/",
    pattern = "2014",
    full.names = TRUE
  )
)

# Gamma calibration data
load("data/calibration/gamma_calibration_1043_sites.Rdata")

# Process rasters ---------------------------------------------------------

# Aggregate secondary vegetation raster by factor of 4 (increases pixel size)
sec_veg_age_agg <- aggregate(
  sec_veg_age_rst,
  fact = 4, fun = mean, na.rm = TRUE
)

# Merge AGB tiles into single raster
combined_agb_raster <- do.call(merge, agb_rasters)

# Resample all rasters to match sec_veg resolution
resampled_agb_raster <- resample(combined_agb_raster, sec_veg_age_agg, method = "near")
resampled_pq_rst <- resample(pq_rst, sec_veg_age_agg, method = "near")

# Convert gamma data to raster and resample
raster_gamma <- as(
  st_rasterize(calib_df %>% dplyr::select(site_reg_gamma, geometry)),
  "SpatRaster"
)
resampled_gamma <- resample(raster_gamma, sec_veg_age_agg, method = "near")

# Save processed rasters --------------------------------------------------
writeRaster(sec_veg_age_agg, "data/calibration/sec_veg_age_agg.tif", overwrite = TRUE)
writeRaster(resampled_agb_raster, "data/calibration/resampled_agb_aggregated.tif", overwrite = TRUE)
writeRaster(resampled_pq_rst, "data/calibration/resampled_pq_rst.tif", overwrite = TRUE)
writeRaster(resampled_gamma, "data/calibration/resampled_gamma.tif", overwrite = TRUE)

# Create analysis mask and calculate area ---------------------------------

# Create mask for valid pixels (non-NA and non-zero for agb and sec)
mask_layer <- !is.na(resampled_agb_raster) &
  !is.na(sec_veg_age_agg) &
  !is.na(resampled_gamma) &
  resampled_agb_raster != 0 &
  sec_veg_age_agg != 0

# Calculate number of valid pixels
n_pixels <- global(mask_layer, fun = "sum", na.rm = TRUE)
cat("Number of valid pixels:", n_pixels$sum, "\n")

# Calculate total area in hectares
cell_areas <- cellSize(mask_layer, unit = "ha")
total_area_ha <- global(cell_areas * mask_layer, fun = "sum", na.rm = TRUE)
cat("Total area (hectares):", total_area_ha$sum, "\n")

# Convert rasters to dataframe --------------------------------------------

# Extract values and convert AGB from Mg/ha to Mg C/ha
combined_df <- data.frame(
  agb = values(resampled_agb_raster)[, 1] / 2 * (44 / 12), # Convert to carbon
  sec = values(sec_veg_age_agg)[, 1],
  gamma = values(resampled_gamma)[, 1],
  pq = values(resampled_pq_rst)[, 1]
) %>%
  filter(agb != 0, sec != 0)

# Save final dataframe
save(combined_df, file = "data/calibration/combined_df.Rdata")

# End timer
toc(log = TRUE)

cat("\nData preparation complete!\n")
cat("Output saved to: data/calibration/combined_df.Rdata\n")
