library(sf)
library(terra)
library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)

in_file <- "data/raw/mapbiomas/land_use_cover/COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-2017.tif"
land_use_rst <- rast(in_file)

# Read in secondary vegetation raster
sec_veg_area_rst <- rast(list.files(
  "data/raw/mapbiomas/secondary_vegetation_area/",
  full.names = TRUE
))

# Crop the the Amazonia subset
sec_veg_area_rst <- sec_veg_area_rst |>
  crop(land_use_rst) |>
  mask(land_use_rst)

# Create dummy for secondary vegetation class
sec_veg_area_rst[sec_veg_area_rst != 303] <- 0
sec_veg_area_rst[sec_veg_area_rst == 303] <- 1

# Aggregate into pixel shares
sec_veg_area_rst <- aggregate(
  sec_veg_area_rst,
  fact = 2250,
  fun = sum,
  na.rm = TRUE
) / (2250^2)

# Rename raster layers
names(sec_veg_area_rst) <- sapply(
  1:nlyr(sec_veg_area_rst),
  function(x) glue("share_sec_veg_{2005 + x}")
)

# Write raster
out_file <- "data/processed/share_sec_veg.tif"
writeRaster(sec_veg_area_rst, out_file, overwrite = TRUE)
