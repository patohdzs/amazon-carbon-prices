# Prepare raster data for AGB-secondary vegetation analysis

library(sf)
library(terra)
library(tidyverse)


# Study extent from calibration
calibration_path <- "data/calibration/gamma_calibration_1043_sites.Rdata"
load(calibration_path)
calib_vect <- vect(calib_1043)

# Load and crop rasters
sec_veg_age_candidates <- list.files(
  "data/raw/mapbiomas/secondary_vegetation_age/",
  pattern = "2017.*\\.tif$",
  full.names = TRUE
)
sec_veg_age_file <- sec_veg_age_candidates[
  grepl("mapbiomas_brazil_sv_age_2017\\.tif$", sec_veg_age_candidates)
]
if (length(sec_veg_age_file) == 0 && length(sec_veg_age_candidates) == 1) {
  sec_veg_age_file <- sec_veg_age_candidates
}
if (length(sec_veg_age_file) != 1) {
  stop(
    sprintf(
      "Expected exactly one secondary vegetation age raster, got %d: %s",
      length(sec_veg_age_candidates),
      paste(sec_veg_age_candidates, collapse = ", ")
    )
  )
}
sec_veg_age_rst <- rast(sec_veg_age_file)
sec_veg_age_rst <- crop(sec_veg_age_rst, ext(project(calib_vect, sec_veg_age_rst)))

# Crop each AGB tile first so merge only handles tiles intersecting the study area
ext_overlaps <- function(a, b) {
  av <- as.vector(a)
  bv <- as.vector(b)
  !(av[2] < bv[1] || av[1] > bv[2] || av[4] < bv[3] || av[3] > bv[4])
}

study_extent_ll <- ext(sec_veg_age_rst)
agb_files <- list.files(
  "data/raw/esa/above_ground_biomass/",
  pattern = "_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2017-fv3.0.tif",
  full.names = TRUE
)
agb_rasters <- map(agb_files, function(path) {
  tile <- rast(path)
  if (!ext_overlaps(ext(tile), study_extent_ll)) {
    return(NULL)
  }
  crop(tile, study_extent_ll)
}) %>%
  compact()

if (length(agb_rasters) == 0) {
  stop("No AGB tiles overlap the study extent.")
}

combined_agb <- do.call(merge, agb_rasters)
rm(agb_rasters)

# Aggregate sec veg 4x, resample AGB and gamma to match; write direct to disk
out_dir <- "data/calibration/carbon_capture_curves"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sec_veg_age_agg <- aggregate(sec_veg_age_rst, fact = 4, fun = mean, na.rm = TRUE,
  filename = file.path(out_dir, "sec_veg_age_agg.tif"), overwrite = TRUE)
rm(sec_veg_age_rst)

agb_resamp <- resample(combined_agb, sec_veg_age_agg, method = "near",
  filename = file.path(out_dir, "resampled_agb_aggregated.tif"), overwrite = TRUE)
rm(combined_agb)

gamma_params_path <- "data/calibration/productivity_params_1043.csv"
gamma_params <- read.csv(gamma_params_path)
gamma_params$id <- seq_len(nrow(gamma_params))
calib_bayesian <- calib_1043 %>%
  select(id, geometry) %>%
  left_join(gamma_params %>% select(id, gamma_fit), by = "id")

gamma_resamp <- rasterize(
  project(vect(calib_bayesian), sec_veg_age_agg),
  sec_veg_age_agg,
  field = "gamma_fit",
  touches = TRUE,
  filename = file.path(out_dir, "resampled_gamma.tif"),
  overwrite = TRUE
)

# Chunked extraction (read from disk)
n_rows <- nrow(sec_veg_age_agg)
chunk_size <- 500
chunks <- split(seq_len(n_rows), ceiling(seq_len(n_rows) / chunk_size))

combined_list <- vector("list", length(chunks))
for (i in seq_along(chunks)) {
  rows <- chunks[[i]]
  n <- max(rows) - min(rows) + 1
  agb_vals <- values(agb_resamp, row = min(rows), nrows = n)[, 1]
  sec_vals <- values(sec_veg_age_agg, row = min(rows), nrows = n)[, 1]
  gam_vals <- values(gamma_resamp, row = min(rows), nrows = n)[, 1]

  valid <- !is.na(agb_vals) & !is.na(sec_vals) & !is.na(gam_vals) &
    agb_vals != 0 & sec_vals != 0

  if (any(valid)) {
    combined_list[[i]] <- data.frame(
      agb = agb_vals[valid] / 2 * (44 / 12),  # Mg/ha biomass -> Mg CO2e/ha
      sec = sec_vals[valid],
      gamma = gam_vals[valid]
    )
  }
  rm(agb_vals, sec_vals, gam_vals) # to free memory, very heavy job!
  if (i %% 100 == 0) gc()
}

combined_list <- combined_list[!sapply(combined_list, is.null)]
combined_df <- do.call(rbind, combined_list)
rm(combined_list)
gc()

save(combined_df, file = file.path(out_dir, "combined_df.Rdata"))
