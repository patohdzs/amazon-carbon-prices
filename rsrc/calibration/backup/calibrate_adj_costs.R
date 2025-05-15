library(sf)
library(ggplot2)
library(terra)
library(tidyverse)
library(conflicted)
library(stargazer)

conflicts_prefer(dplyr::filter)
conflicts_prefer(terra::extract)

# Load raster files
sec_veg_age_rst <- rast(
  list.files(
    "data/raw/mapbiomas/secondary_vegetation_age/",
    pattern = "2017",
    full.names = TRUE
  )
)

clean_mapbiomas <- rast("data/clean/land_use_cover_2000.tif")

# Select Amazonia subset
sec_veg_age_rst <- sec_veg_age_rst |>
  crop(clean_mapbiomas) |>
  mask(clean_mapbiomas)

# Compute pixel areas
pixel_areas <- cellSize(sec_veg_age_rst, unit = "ha")

# Calculate the total area for each secondary vegetation age
ref <- zonal(pixel_areas, sec_veg_age_rst, sum, na.rm = TRUE) %>%
  slice(-1) %>%
  mutate(prop = area / sum(area, na.rm = TRUE))

print(ref)

# Get mean yearly area reforested
avg_ref <- mean(ref[1:10, 2])

# Compute implied adjustment cost params (MC = MR condition)
price_v <- 52
zeta_v <- price_v / avg_ref

print("Zeta (reforestation)")
print(zeta_v)

# Deforestation cost calibration
# We transform to dollars using an FX rate of 4.14 (December 2019)
deforestation_cost <- 1614.54 / 4.14

# The forest area in 2008 represented 72% of the Legal Amazon area,
#   which covers 501,506,775 ha, so the transition in hectares is
#   0.065*0.72*501,506,775 = 23,470,517,
#   resulting in an annual average of 2,347,052 ha.
transition_area <- (0.065 * 0.72 * 501506775) / (2017 - 2008 + 1)
zeta_u <- deforestation_cost / transition_area
print("Zeta (deforestation)")
print(zeta_u)

# Alternative value based on https://shorturl.at/jgp9i
zeta_alt <- 483 / transition_area
