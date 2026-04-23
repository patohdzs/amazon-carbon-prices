#!/usr/bin/env Rscript
# investigate_site_matching.R
#
# Compares two algorithms for assigning each of 1043 fine sites to one of 78
# coarse sites:
#   Method A (current): nearest coarse-site centroid         [line 79 of site_gamma_reg_1043.R]
#   Method B (alt):     coarse site with largest polygon overlap with the fine site
#
# Produces site_matching_comparison.png in the project root.

library(sf)
library(tidyverse)
library(ggplot2)

# ── 1. Load data ───────────────────────────────────────────────────────────────
message("Loading fine sites (1043)...")
load("data/calibration/calibration_1043_sites.Rdata")   # → calib_df
fine_sites <- calib_df %>% mutate(fine_idx = row_number())

message("Loading coarse sites (78)...")
load("data/calibration/calibration_78_sites.Rdata")     # → calib_df
coarse_sites <- calib_df %>% mutate(coarse_idx = row_number())

# Align CRS (reproduce the transform on line 75 of site_gamma_reg_1043.R)
fine_sites <- st_transform(fine_sites, st_crs(coarse_sites))

stopifnot(nrow(fine_sites) == 1043, nrow(coarse_sites) == 78)

# ── 2. Method A – centroid distance (current algorithm) ───────────────────────
message("Method A: nearest centroid...")
match_centroid <- st_nearest_feature(
  st_centroid(fine_sites),
  st_centroid(coarse_sites)
)

# ── 3. Method B – largest area overlap ────────────────────────────────────────
message("Method B: largest area overlap (st_join largest=TRUE)...")
overlap_join <- st_join(
  fine_sites  %>% select(fine_idx),
  coarse_sites %>% select(coarse_idx),
  join    = st_intersects,
  largest = TRUE          # keep only the overlapping coarse site with the largest area
)
match_overlap <- overlap_join$coarse_idx   # NA if no coarse site overlaps at all

# ── 4. Compare ─────────────────────────────────────────────────────────────────
comparison <- tibble(
  fine_idx = fine_sites$fine_idx,
  method_A = match_centroid,
  method_B = match_overlap,
  status   = case_when(
    is.na(method_B)      ~ "No overlap found",
    method_A == method_B ~ "Agree",
    TRUE                 ~ "Disagree"
  )
)

n_agree    <- sum(comparison$status == "Agree")
n_disagree <- sum(comparison$status == "Disagree")
n_no_ovlp  <- sum(comparison$status == "No overlap found")

message(sprintf(
  "\n--- Results ---\nFine sites : %d\nAgree      : %d (%.1f%%)\nDisagree   : %d (%.1f%%)\nNo overlap : %d",
  nrow(comparison),
  n_agree,    100 * n_agree    / nrow(comparison),
  n_disagree, 100 * n_disagree / nrow(comparison),
  n_no_ovlp
))

# Print the sites that disagree
if (n_disagree > 0) {
  cat("\nDisagreeing fine sites (fine_idx | method_A coarse_idx | method_B coarse_idx):\n")
  print(comparison %>% filter(status == "Disagree"))
}

# ── 5. Plot ────────────────────────────────────────────────────────────────────
message("Building plot...")

plot_sf <- fine_sites %>%
  left_join(comparison %>% select(fine_idx, status), by = "fine_idx") %>%
  mutate(status = factor(status, levels = c("Agree", "Disagree", "No overlap found")))

status_colours <- c(
  "Agree"            = "#4daf4a",
  "Disagree"         = "#e41a1c",
  "No overlap found" = "#ff7f00"
)

subtitle_text <- sprintf(
  "%d / %d fine sites (%.1f%%) assigned differently — %d have no overlap with any coarse site",
  n_disagree, nrow(comparison), 100 * n_disagree / nrow(comparison), n_no_ovlp
)

p <- ggplot() +
  geom_sf(data = plot_sf,    aes(fill = status), colour = NA) +
  geom_sf(data = coarse_sites, fill = NA, colour = "black", linewidth = 0.35) +
  scale_fill_manual(values = status_colours, name = NULL, drop = FALSE) +
  labs(
    title    = "Centroid-distance vs. largest-area-overlap site matching",
    subtitle = subtitle_text,
    caption  = "Black outlines = 78 coarse sites"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "bottom",
    plot.subtitle    = element_text(size = 10, colour = "grey40"),
    plot.caption     = element_text(size = 9,  colour = "grey50")
  )

out_path <- "site_matching_comparison.png"
ggsave(out_path, p, width = 11, height = 9, dpi = 150)
message("Saved: ", out_path)
