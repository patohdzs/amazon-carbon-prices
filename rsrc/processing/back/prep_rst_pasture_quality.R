library(sf)
library(glue)
library(tictoc)
library(terra)
library(conflicted)

conflicts_prefer(terra::extract())

# START TIMER
tic(msg = "prep_biome_class_rasters.R script", log = TRUE)

in_files <- "data/raw/mapbiomas/pq-export/mapbiomas-brazil-collection-80-pasture-quality-amazonia-2017.tif"

for (level in c(1, 2, 3)) {
  # Read in pasture quality raster
  pasture_quality_rst <- rast(in_files)

  # Create dummy for pasture quality class
  pasture_quality_rst[pasture_quality_rst != level] <- 0
  pasture_quality_rst[pasture_quality_rst == level] <- 1

  # Aggregate into pixel shares
  pasture_quality_rst <- aggregate(
    pasture_quality_rst,
    fact = 2250,
    fun = sum,
    na.rm = TRUE
  ) / (2250^2)

  # Rename raster layers
  names(pasture_quality_rst) <- sapply(
    1:nlyr(pasture_quality_rst),
    function(x) glue("share_pasture_quality_{level}_2017")
  )

  # Write rasters
  out_file <- glue("data/processed/share_pq_{level}_1043_sites.tif")
  writeRaster(pasture_quality_rst, out_file, overwrite = TRUE)
}
