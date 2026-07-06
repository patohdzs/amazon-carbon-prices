# =============================================================================
# extract_z_zbar_forest_all_years.R
#
# STANDALONE EXTRACTOR — raw MapBiomas tiles + the saved 1043-site grid  ->
#   z_<year>, zbar_<year>, area_forest_<year>
# for EVERY year present in the MapBiomas annual land-use/cover tiles
# (the in-repo pipeline hard-codes only 1995 / 2008 / 2017).
#
# It takes the project's already-computed site grid (grid_1043_sites.geojson)
# as a raw input and computes the three land-use variables onto those exact
# cells for all years. Because the grid is given, NONE of the theta / biome
# machinery is needed: the only inputs are the MapBiomas tiles and the geojson,
# and the output is exactly the 1043 calibration sites (with matching `id`).
#
# Faithful to the pipeline's numbers:
#   - shares  = aggregate(30 m class dummy, fact = 2250) / 2250^2   (prep_rst_land_use.R)
#   - area_ha = cellSize() of the aggregated ~67.5 km cell          (prep_rst_biome.R)
#   - z = agricultural-use area ; zbar = forest area + z            (calibrate_*_model.R)
#
# INPUTS (raw):
#   data/raw/mapbiomas/land_use_cover/COLECAO_5_..._AMAZONIA-<year>.tif  (1 per year)
#   data/calibration/grid_1043_sites.geojson   (site polygons + `id`, EPSG:5880)
#
# OUTPUTS:
#   <out_dir>/z_zbar_area_forest_all_years.csv      (one row per site)
#   <out_dir>/grid_all_years_sites.geojson          (site polygons + values)
#
# USAGE:  run from the project root, e.g.
#   Rscript extract_z_zbar_forest_all_years.R
#
# Swap `grid_geojson` for grid_78_sites.geojson to get the 78-site version.
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
})

# ------------------------------- CONFIG --------------------------------------
proj_root    <- getwd()                          # run from repo root
raw_dir      <- file.path(proj_root, "data", "raw")
lu_dir       <- file.path(raw_dir, "mapbiomas", "land_use_cover")
grid_geojson <- file.path(proj_root, "data", "calibration", "grid_1043_sites.geojson")

out_dir      <- file.path(proj_root, "data", "calibration")
out_csv      <- file.path(out_dir, "z_zbar_area_forest_all_years.csv")
out_geojson  <- file.path(out_dir, "grid_all_years_sites.geojson")

agg_factor   <- 2250     # 2250 * 30 m = 67.5 km site cells (matches pipeline)

# MapBiomas Collection 5 class codes (see docs/.../mapbiomasClass_id_legend.pdf).
# These reproduce the pipeline exactly. z = agricultural_use area; the third
# pipeline class ("other") is not needed for z / zbar / area_forest.
class_codes <- list(
  forest           = 3,                  # Forest Formation
  agricultural_use = c(15, 20, 39, 41)   # Pasture, Sugar cane, Soybean, Other temp. crops
)

# ---------------------- 1. Discover MapBiomas year tiles ---------------------
if (!dir.exists(lu_dir))      stop("MapBiomas dir not found: ", lu_dir)
if (!file.exists(grid_geojson)) stop("Site grid not found: ", grid_geojson)

tile_files <- list.files(
  lu_dir,
  pattern    = "COLECAO_5_DOWNLOADS_COLECOES_ANUAL_AMAZONIA_AMAZONIA-\\d{4}\\.tif$",
  full.names = TRUE
)
if (length(tile_files) == 0L) stop("No MapBiomas annual tiles found in: ", lu_dir)

years <- as.integer(regmatches(
  basename(tile_files),
  regexpr("(?<=AMAZONIA-)\\d{4}", basename(tile_files), perl = TRUE)
))
o <- order(years); years <- years[o]; tile_files <- tile_files[o]
message(sprintf("Found %d MapBiomas tiles spanning %d-%d.",
                length(years), min(years), max(years)))

# ----------------------------- 2. Read site grid -----------------------------
sites <- st_read(grid_geojson, quiet = TRUE)     # `id` + polygons (EPSG:5880)
message(sprintf("Site grid: %d sites from %s", nrow(sites), basename(grid_geojson)))

# ---------------- 3. Per-year, per-class shares (forest, ag use) -------------
# For each year: reclassify 30 m pixels to a class dummy (1 = class, 0 = not),
# then aggregate to the ~67.5 km grid -> class share. Numerically identical to
# prep_rst_land_use.R (sum of dummies / agg_factor^2).
share_layers <- list()
for (i in seq_along(years)) {
  y <- years[i]
  r <- rast(tile_files[[i]])
  for (cl in names(class_codes)) {
    dummy <- r %in% class_codes[[cl]]            # 1 = in class, 0 = otherwise
    agg   <- aggregate(dummy, fact = agg_factor, fun = sum, na.rm = TRUE) / (agg_factor^2)
    # Defensive: keep every layer on one common grid (tiles are co-registered).
    if (length(share_layers) &&
        !isTRUE(compareGeom(agg, share_layers[[1]], stopOnError = FALSE)))
      agg <- resample(agg, share_layers[[1]])
    names(agg) <- sprintf("share_%s_%d", cl, y)
    share_layers[[length(share_layers) + 1L]] <- agg
  }
  rm(r); gc()
  message(sprintf("  [%d/%d] processed %d", i, length(years), y))
}

# True geographic area (ha) of each ~67.5 km site cell (pipeline uses cellSize).
site_area <- cellSize(share_layers[[1]], unit = "ha")
names(site_area) <- "pixel_area_ha"

# ---------------- 4. Sample shares + area onto the 1043 cells ----------------
stk <- rast(c(share_layers, list(site_area)))
# Each site polygon is exactly one aggregated cell; sample at an interior point.
pts  <- vect(st_transform(st_point_on_surface(sites), crs(stk)))
vals <- terra::extract(stk, pts, ID = FALSE)     # one row per site, layers as cols
sites <- cbind(sites, vals)
sites$site_area_ha <- sites$pixel_area_ha

# ---------------- 5. Shares -> areas (ha); then z, area_forest, zbar ---------
for (y in years) {
  for (cl in names(class_codes)) {
    sc <- sprintf("share_%s_%d", cl, y)
    ac <- sprintf("area_%s_%d", cl, y)
    sites[[ac]] <- sites[[sc]] * sites$site_area_ha
  }
  sites[[sprintf("z_%d", y)]]    <- sites[[sprintf("area_agricultural_use_%d", y)]]
  sites[[sprintf("zbar_%d", y)]] <- sites[[sprintf("area_forest_%d", y)]] +
                                    sites[[sprintf("z_%d", y)]]
}

# ----------------------------- 6. Write outputs ------------------------------
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Centroid lon/lat (WGS84) for convenience in downstream projects.
cent <- st_coordinates(st_centroid(st_transform(sites, 4326)))
sites$lon <- cent[, 1]
sites$lat <- cent[, 2]

year_cols <- as.vector(vapply(years, function(y) c(
  sprintf("area_forest_%d", y),
  sprintf("z_%d", y),
  sprintf("zbar_%d", y)
), character(3)))

csv_cols <- c("id", "lon", "lat", "site_area_ha", year_cols)

st_write(sites[, c("id", year_cols)], out_geojson,
         driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
write.csv(st_drop_geometry(sites[, csv_cols]), out_csv, row.names = FALSE)

message("Wrote:\n  ", out_csv, "\n  ", out_geojson)

# =============================================================================
# CHECK — the 1995 / 2008 / 2017 columns produced here should match
# data/calibration/calibration_1043_sites.csv (same `id`), since both use the
# same shares (aggregate fact = 2250) and the same cellSize areas. The only
# addition is the remaining years.
# =============================================================================
