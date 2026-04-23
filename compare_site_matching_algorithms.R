#!/usr/bin/env Rscript

# Compare two ways of assigning each of the 1043 fine sites to one of the 78 coarse sites:
# 1) Current algorithm in rsrc/calibration/site_gamma_reg_1043.R line 79:
#    nearest coarse centroid to each fine centroid.
# 2) Alternative algorithm:
#    coarse polygon with the largest overlap area with each fine polygon.
#
# Outputs:
# - site_assignment_algorithm_comparison.png
# - output/site_assignment_algorithm_comparison.csv

OVERLAP_SHARE_THRESHOLD <- 0.01

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
})

message("Loading calibration polygons...")
load("data/calibration/calibration_1043_sites.Rdata")
fine_sites <- calib_df
load("data/calibration/calibration_78_sites.Rdata")
coarse_sites <- calib_df

if (!inherits(fine_sites, "sf") || !inherits(coarse_sites, "sf")) {
  stop("Expected sf objects in both calibration files.")
}

fine_sites <- fine_sites %>% mutate(fine_id = dplyr::row_number())
coarse_sites <- coarse_sites %>% mutate(coarse_id = dplyr::row_number())

message("Aligning CRS (same as line 75 in site_gamma_reg_1043.R)...")
fine_sites <- st_transform(fine_sites, st_crs(coarse_sites))

message("Running main method: largest polygon overlap area...")
fine_valid <- st_make_valid(fine_sites)
coarse_valid <- st_make_valid(coarse_sites)
fine_areas <- as.numeric(st_area(fine_valid))

intersects_list <- st_intersects(fine_valid, coarse_valid)
overlap_assignment <- rep(NA_integer_, nrow(fine_valid))
best_overlap_area <- rep(NA_real_, nrow(fine_valid))
best_overlap_share <- rep(NA_real_, nrow(fine_valid))

for (i in seq_len(nrow(fine_valid))) {
  candidates <- intersects_list[[i]]
  if (length(candidates) == 0) {
    next
  }

  overlap_geoms <- st_intersection(
    fine_valid[i, "fine_id"],
    coarse_valid[candidates, "coarse_id"]
  )

  if (nrow(overlap_geoms) == 0) {
    next
  }

  overlap_areas <- as.numeric(st_area(overlap_geoms))
  overlap_shares <- overlap_areas / fine_areas[i]
  best_local_idx <- which.max(overlap_areas)

  best_overlap_area[i] <- overlap_areas[best_local_idx]
  best_overlap_share[i] <- overlap_shares[best_local_idx]

  if (!is.na(best_overlap_share[i]) && best_overlap_share[i] >= OVERLAP_SHARE_THRESHOLD) {
    overlap_assignment[i] <- overlap_geoms$coarse_id[best_local_idx]
  }
}

message("Running fallback method: nearest centroid...")
nearest_centroid_assignment <- st_nearest_feature(
  st_centroid(fine_sites),
  st_centroid(coarse_sites)
)

hybrid_assignment <- dplyr::coalesce(overlap_assignment, nearest_centroid_assignment)

comparison <- tibble(
  fine_id = fine_sites$fine_id,
  nearest_centroid = nearest_centroid_assignment,
  largest_overlap = overlap_assignment,
  overlap_then_nearest = hybrid_assignment,
  max_overlap_area = best_overlap_area,
  max_overlap_share = best_overlap_share
) %>%
  mutate(
    status = case_when(
      is.na(largest_overlap) ~ "No qualifying overlap",
      nearest_centroid == largest_overlap ~ "Same assignment",
      TRUE ~ "Different assignment"
    )
  )

n_total <- nrow(comparison)
n_same <- sum(comparison$status == "Same assignment")
n_diff <- sum(comparison$status == "Different assignment")
n_none <- sum(comparison$status == "No qualifying overlap")

message(sprintf("Total fine sites: %d", n_total))
message(sprintf("Same assignment: %d (%.1f%%)", n_same, 100 * n_same / n_total))
message(sprintf("Different assignment: %d (%.1f%%)", n_diff, 100 * n_diff / n_total))
message(sprintf(
  "No qualifying overlap (max overlap share < %.0f%% or no intersection): %d",
  100 * OVERLAP_SHARE_THRESHOLD,
  n_none
))

if (n_diff > 0) {
  message("\nFirst 20 differing assignments (fine_id, nearest_centroid, largest_overlap):")
  print(
    comparison %>%
      filter(status == "Different assignment") %>%
      select(fine_id, nearest_centroid, largest_overlap) %>%
      head(20)
  )
}

if (!dir.exists("output")) {
  dir.create("output", recursive = TRUE)
}

write.csv(
  comparison,
  file = "output/site_assignment_algorithm_comparison.csv",
  row.names = FALSE
)
message("Saved: output/site_assignment_algorithm_comparison.csv")

comparison_with_geom <- fine_sites %>%
  select(fine_id) %>%
  left_join(comparison, by = "fine_id")

plot_data <- bind_rows(
  comparison_with_geom %>%
    mutate(
      method = "Nearest centroid",
      coarse_assignment = as.factor(nearest_centroid)
    ) %>%
    select(method, coarse_assignment, status, geometry),
  comparison_with_geom %>%
    mutate(
      method = sprintf("Largest overlap area (>= %.0f%%)", 100 * OVERLAP_SHARE_THRESHOLD),
      coarse_assignment = as.factor(largest_overlap)
    ) %>%
    select(method, coarse_assignment, status, geometry),
  comparison_with_geom %>%
    mutate(
      method = "Largest overlap, then nearest centroid fallback",
      coarse_assignment = as.factor(overlap_then_nearest),
      status = NA_character_
    ) %>%
    select(method, coarse_assignment, status, geometry)
)

subtitle_text <- sprintf(
  "Same: %d (%.1f%%) | Different: %d (%.1f%%) | No overlap: %d",
  n_same,
  100 * n_same / n_total,
  n_diff,
  100 * n_diff / n_total,
  n_none
)

p <- ggplot(plot_data) +
  geom_sf(aes(fill = coarse_assignment), color = NA) +
  geom_sf(
    data = comparison_with_geom %>% filter(status == "Different assignment"),
    fill = NA,
    color = "red",
    linewidth = 0.18
  ) +
  geom_sf(
    data = coarse_sites,
    fill = NA,
    color = "black",
    linewidth = 0.12
  ) +
  facet_wrap(~method, ncol = 3) +
  scale_fill_viridis_d(option = "turbo", na.value = "grey85", guide = "none") +
  labs(
    title = "Fine-to-coarse assignment across three mappings",
    subtitle = subtitle_text,
    caption = sprintf(
      "Black lines: coarse-site boundaries. Red outlines: fine sites where nearest-centroid and largest-overlap methods differ. Largest-overlap match requires >= %.0f%% overlap of fine polygon area.",
      100 * OVERLAP_SHARE_THRESHOLD
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

figure_path <- "site_assignment_algorithm_comparison.png"
ggsave(figure_path, plot = p, width = 18, height = 7.5, dpi = 180, bg = "white")
message("Saved: ", figure_path)
