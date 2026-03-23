# Same as 01_prepare_data.R but using mapbiomas_brazil_sv_age_2017.tif
# Output: combined_df_2.Rdata

library(sf)
library(stars)
library(terra)
library(tidyverse)

load("data/processed/gamma_calibration_1043_sites.Rdata")
study_extent <- ext(vect(calib_1043))

# Brazil raster is already 4x4 aggregated natively - no aggregation needed
sec_veg_age_rst <- rast("data/raw/mapbiomas/secondary_vegetation_age/mapbiomas_brazil_sv_age_2017.tif")
sec_veg_age_rst <- crop(sec_veg_age_rst, study_extent)

agb_rasters <- map(
  list.files("data/raw/esa/above_ground_biomass/",
    pattern = "_ESACCI-BIOMASS-L4-AGB-MERGED-100m-2017-fv3.0.tif", full.names = TRUE),
  rast
)
combined_agb <- do.call(merge, agb_rasters)
combined_agb <- crop(combined_agb, study_extent)
rm(agb_rasters)

out_dir <- "data/calibration/carbon_capture_curves"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# No aggregation - Brazil raster already 4x4 native
sec_veg_age_agg <- writeRaster(sec_veg_age_rst,
  filename = file.path(out_dir, "sec_veg_age_agg_2.tif"), overwrite = TRUE)
rm(sec_veg_age_rst)

agb_resamp <- resample(combined_agb, sec_veg_age_agg, method = "near",
  filename = file.path(out_dir, "resampled_agb_aggregated_2.tif"), overwrite = TRUE)
rm(combined_agb)

gamma_params <- read.csv("data/calibration/hmc/productivity_params_1043.csv")
gamma_params$id <- seq_len(nrow(gamma_params))
calib_bayesian <- calib_1043 %>%
  select(id, geometry) %>%
  left_join(gamma_params %>% select(id, gamma_fit), by = "id")

raster_gamma <- as(st_rasterize(calib_bayesian %>% select(gamma_fit, geometry)), "SpatRaster")
gamma_resamp <- resample(raster_gamma, sec_veg_age_agg, method = "near",
  filename = file.path(out_dir, "resampled_gamma_2.tif"), overwrite = TRUE)
rm(raster_gamma)

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
      agb = agb_vals[valid] / 2 * (44 / 12),
      sec = sec_vals[valid],
      gamma = gam_vals[valid]
    )
  }
  rm(agb_vals, sec_vals, gam_vals)
  if (i %% 100 == 0) gc()
}

combined_list <- combined_list[!sapply(combined_list, is.null)]
combined_df_2 <- do.call(rbind, combined_list)
rm(combined_list)
gc()

save(combined_df_2, file = file.path(out_dir, "combined_df_2.Rdata"))
