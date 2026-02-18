# Plot valid pixels colored by secondary vegetation age

library(terra)
library(viridis)
library(sf)
library(rnaturalearth)
library(geobr)

load("data/processed/gamma_calibration_1043_sites.Rdata")

# Dissolve site boundaries (project to meters to avoid s2 degenerate artifacts)
# This is an expedient workaround because without it the plot was coloring a subset of
# the gridlines in a weird way.

calib_proj <- st_transform(calib_1043, 5880) |>
  st_set_precision(1) |>
  st_make_valid() |>
  st_collection_extract("POLYGON")
study_boundary <- st_union(calib_proj) |> st_transform(4326)

# Base layers
rivers <- ne_download(scale = 10, type = "rivers_lake_centerlines",
  category = "physical", returnclass = "sf")
rivers_clip <- st_intersection(rivers, study_boundary)

states <- read_state(year = 2000, showProgress = FALSE)
states <- st_transform(states, st_crs(study_boundary))
states_clip <- st_intersection(st_geometry(states), study_boundary)

# Load rasters and build mask (valid = non-NA, non-zero, age <= 30)
cc_dir <- "data/calibration/carbon_capture_curves"
sec_r <- rast(file.path(cc_dir, "sec_veg_age_agg.tif"))
agb_r <- rast(file.path(cc_dir, "resampled_agb_aggregated.tif"))
gam_r <- rast(file.path(cc_dir, "resampled_gamma.tif"))

mask_layer <- !is.na(agb_r) & !is.na(sec_r) & !is.na(gam_r) &
  agb_r != 0 & sec_r != 0 & sec_r <= 30

sec_masked <- terra::mask(sec_r, mask_layer, maskvalues = 0)
sec_masked <- crop(sec_masked, vect(study_boundary))
sec_int <- ceiling(sec_masked)

# Second panel: age >= 6 only
sec_age6 <- terra::mask(sec_int, sec_masked >= 6, maskvalues = 0)
cell_area_ha <- cellSize(sec_int, unit = "ha")
area_all <- global(mask(cell_area_ha, sec_int), fun = "sum", na.rm = TRUE)$sum
area_age6 <- global(mask(cell_area_ha, sec_age6), fun = "sum", na.rm = TRUE)$sum

# Plot
dir.create("output/figures/carbon_capture", recursive = TRUE, showWarnings = FALSE)
png("output/figures/carbon_capture/valid_pixels_age_map.png",
  width = 6000, height = 2400, res = 300)
par(mfrow = c(1, 2))

plot(sec_int,
  main = paste0("All Valid (", format(round(area_all), big.mark = ","), " ha)"),
  col = viridis(100, direction = -1),
  plg = list(title = "Age (years)"))
plot(st_geometry(rivers_clip), col = "steelblue", lwd = 1.2, add = TRUE)
plot(states_clip, border = "black", lwd = 1.5, add = TRUE)
plot(study_boundary, border = "black", lwd = 1.5, add = TRUE)

plot(sec_age6,
  main = paste0("Age >= 6 (", format(round(area_age6), big.mark = ","), " ha)"),
  col = viridis(100, direction = -1),
  plg = list(title = "Age (years)"))
plot(st_geometry(rivers_clip), col = "steelblue", lwd = 1.2, add = TRUE)
plot(states_clip, border = "gray40", lwd = 0.8, add = TRUE)
plot(study_boundary, border = "black", lwd = 1.5, add = TRUE)

dev.off()
