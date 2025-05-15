library(sf)
library(terra)
library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)

in_file <- "data/raw/mapbiomas/land_use_cover/COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-2017.tif"
land_use_rst <- rast(in_file)

# Read in secondary vegetation age rasters
sec_veg_age_rst <- rast(
  list.files(
    "data/raw/mapbiomas/secondary_vegetation_age/",
    full.names = TRUE
  )
)

# Crop the the Amazonia subset
sec_veg_age_rst <- sec_veg_age_rst |>
  crop(land_use_rst) |>
  mask(land_use_rst)


# Aggregate into pixels
sec_veg_age_rst <- aggregate(
  sec_veg_age_rst,
  fact = 2250,
  fun = median,
  na.rm = TRUE
)

# Rename raster layer
names(sec_veg_area_rst) <- sapply(
  1:nlyr(sec_veg_area_rst),
  function(x) glue("sec_veg_age_{2005 + x}")
)

# Write raster
out_file <- "data/processed/sec_veg_age.tif"
writeRaster(sec_veg_age_rst, out_file, overwrite = TRUE)
