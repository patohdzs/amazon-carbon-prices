library(sf)
library(terra)
library(tidyverse)

stability_rst <- rast("data/raw/bio/stability_regions/STAB2050_MB.tif")

# Convert the raster to a dataframe
raster_df <- as.data.frame(stability_rst, xy = TRUE) %>% as_tibble()

raster_df <- raster_df %>%
  mutate(
    region = cut(
      STAB2050_MB.tif,
      breaks = c(-Inf, 0.5, 1.5, 2.5),
      labels = c("Stable Forest", "Bistable", "Stable Savanna"),
      include.lowest = TRUE
    )
  )


# Basic ggplot2 plot
fig <- ggplot(data = raster_df, aes(x = x, y = y, fill = region)) +
  geom_raster() +
  scale_fill_viridis_d() +
  coord_equal() +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Region"
  ) +
  theme(legend.position = "bottom")

ggsave(filename = "results/bio/stability_regions.pdf", plot = fig)


clean_mapbiomas_rst <- rast("data/clean/land_use_cover_2000.tif")


stability_rst <- crop(stability_rst, vect(calib_df))
ovr <- mask(stability_rst, vect(calib_df))
