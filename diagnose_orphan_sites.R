#!/usr/bin/env Rscript
# Diagnose the two fine sites that method B couldn't assign.

library(sf)
library(tidyverse)

load("data/calibration/calibration_1043_sites.Rdata"); fine_sites   <- calib_df %>% mutate(fine_idx = row_number())
load("data/calibration/calibration_78_sites.Rdata");   coarse_sites <- calib_df %>% mutate(coarse_idx = row_number())
fine_sites <- st_transform(fine_sites, st_crs(coarse_sites))

# Method B
overlap_join  <- st_join(fine_sites %>% select(fine_idx),
                         coarse_sites %>% select(coarse_idx),
                         join = st_intersects, largest = TRUE)
orphan_idx    <- which(is.na(overlap_join$coarse_idx))

# Method A assignment for the orphans
match_centroid <- st_nearest_feature(st_centroid(fine_sites), st_centroid(coarse_sites))

# Centroids
fine_ctr   <- st_centroid(fine_sites)   %>% st_coordinates()
coarse_ctr <- st_centroid(coarse_sites) %>% st_coordinates()

# Distance from each orphan fine centroid to its assigned coarse centroid (method A),
# compared to distance to the nearest coarse-site polygon BOUNDARY (in degrees -> metres).
coarse_union <- st_union(coarse_sites)

for (i in orphan_idx) {
  f_geom <- st_geometry(fine_sites)[i]
  f_ctr  <- st_centroid(f_geom)
  assigned_c <- match_centroid[i]
  c_ctr  <- st_centroid(st_geometry(coarse_sites)[assigned_c])

  # Is the fine polygon disjoint from the union of coarse polygons?
  is_outside <- as.logical(st_disjoint(f_geom, coarse_union, sparse = FALSE))

  # Distance from fine polygon to nearest point on any coarse polygon
  gap_to_coarse <- as.numeric(st_distance(f_geom, coarse_union))
  # Distance from fine centroid to assigned coarse centroid
  ctr_to_ctr    <- as.numeric(st_distance(f_ctr, c_ctr))

  cat(sprintf(
    "\nfine_idx = %d\n  fine centroid    : (%.3f, %.3f)\n  assigned coarse  : coarse_idx=%d, centroid (%.3f, %.3f)\n  disjoint from all coarse polys? %s\n  gap fine-poly -> nearest coarse-poly : %.4f (CRS units)\n  distance fine-ctr -> assigned coarse-ctr: %.4f\n",
    i,
    fine_ctr[i, "X"], fine_ctr[i, "Y"],
    assigned_c,
    coarse_ctr[assigned_c, "X"], coarse_ctr[assigned_c, "Y"],
    is_outside,
    gap_to_coarse,
    ctr_to_ctr
  ))
}

cat(sprintf("\nCRS: %s\n", st_crs(coarse_sites)$input))
cat(sprintf("Union of coarse sites vs union of fine sites — same extent?\n"))
print(st_bbox(coarse_sites))
print(st_bbox(fine_sites))
