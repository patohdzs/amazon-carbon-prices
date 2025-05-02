# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: MAPS OF 1043 SITES MODEL CALIBRATED VALUES
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

# SETUP

# Load all necessary 
library(sf)           
library(dplyr)          
library(ggplot2)         
library(RColorBrewer)    
library(tictoc)          
library(here)           


# START TIMER
tictoc::tic(msg = "mapsCalibration_1043SitesModel.R script", log = T)

# DATA INPUT

map_basin <- st_read(
  dsn = "data/raw/mapbiomas/basin",
  layer = "BASIN_LEVEL_2_PNRH"
)


# 1043 SITES MODEL CALIBRATION VARIABLES
load(here::here("data/calibration/", "calibration_1043_sites.Rdata"))

# AMAZON BIOME VECTOR DATA
load(here::here("data/clean/amazon_biome.Rdata"))


productivity_params <- read.csv(here::here("data/calibration/", "productivity_params_1043.csv"))%>%
  mutate(id = row_number())

calib_df <- calib_df %>%
  left_join(productivity_params, by = "id")
calib_df <- calib_df %>%
  mutate(
    x_1995 = calib_df$gamma_fit * (zbar_1995 - z_1995),
    x_2008 = calib_df$gamma_fit * (zbar_2008 - z_2008),
    x_2017 = calib_df$gamma_fit * (zbar_2017 - z_2017)
  )

# # GLOBAL MODEL CALIBRATION VARIABLES
# load(here::here("data/calibration/globalModel", "calibration_globalModel.Rdata"))

# stop()
# # extract average price
# p_mean <- rownames(matrixTransition.2prices) %>% as.numeric() %>% mean()

# # clear unnecessary objects
# rm(matrixTransition.2prices, calibration.globalModel)

# DATASET CLEANUP AND PREP

# adjust projection
#calib_df <- sf::st_transform(calib_df, sf::st_crs(4326))
#amazon_biome <- sf::st_transform(amazon_biome, sf::st_crs(calib_df))


calib_df <- calib_df %>%
  st_transform(st_crs(map_basin))

amazon_biome <- amazon_biome %>%
  st_transform(st_crs(map_basin))


# PLOT MAPS

# z_2017
plot <- ggplot2::ggplot(data = calib_df %>%
  dplyr::mutate(z_2017 = ggplot2::cut_number(round(z_2017 / 1000, digits = 1), n = 5, dig.lab = 3))) +
  ggplot2::geom_sf(aes(fill = z_2017)) +
  ggplot2::scale_fill_brewer(name = expression(paste("Z"[2017]^"i", ~"(thousand hectares)")), palette = "YlOrRd") +
  ggplot2::geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::theme(
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    legend.position = "bottom", legend.key.width = unit(3, "cm"), legend.margin = margin(t = -0.5, r = 0, b = 0.2, l = 0, unit = "cm"),
    legend.text = element_text(size = 20, face = "bold")
  )

ggplot2::ggsave(filename = here::here("plots/calibration/1043SitesModel/map_z2017_1043Sites.png"), plot = plot, width = 8, height = 6)



# zbar_2017
plot <- ggplot2::ggplot(data = calib_df %>%
  dplyr::mutate(zbar_2017 = ggplot2::cut_number(round(zbar_2017 / 1000, digits = 1), n = 5, dig.lab = 3))) +
  ggplot2::geom_sf(aes(fill = zbar_2017)) +
  ggplot2::scale_fill_brewer(name = expression(paste(bar(z)^"i", ~"(thousand hectares)")), palette = "YlOrRd") +
  ggplot2::geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::theme(
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    legend.position = "bottom", legend.key.width = unit(3, "cm"), legend.margin = margin(t = -0.5, r = 0, b = 0.2, l = 0, unit = "cm"),
    legend.text = element_text(size = 20, face = "bold")
  )

ggplot2::ggsave(filename = here::here("plots/calibration/1043SitesModel/map_zbar_2017.png"), plot = plot, width = 8, height = 6)

# share_z_2017
plot <- ggplot2::ggplot(data = calib_df %>%
  dplyr::mutate(
    share_z_2017 = (z_2017 / (zbar_2017)) * 100,
    share_z_2017 = cut(round(share_z_2017),
      breaks = c(0, 0.5, 20, 40, 60, 80, 100), include.lowest = T, dig.lab = 3,
      labels = c("[0]", "(0-20]", "(20-40]", "(40-60]", "(60-80]", "(80-100]")
    )
  )) +
  ggplot2::geom_sf(aes(fill = share_z_2017)) +
  ggplot2::scale_fill_manual(name = expression(paste("Z"[2017]^"i", ~"(% of ", bar(z)^"i", ")")), values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
  ggplot2::geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::theme(
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    legend.position = "bottom", legend.key.width = unit(2, "cm"), legend.margin = margin(t = -0.5, r = 0, b = 0.2, l = 0, unit = "cm"),
    legend.text = element_text(size = 20, face = "bold")
  )

ggplot2::ggsave(filename = here::here("plots/calibration/1043SitesModel/map_z2017Share_1043Sites.png"), plot = plot, width = 8, height = 6)

# # A_2017
# ggplot2::ggplot(data = calib_df %>%
#                   dplyr::mutate(A_2017 = (z_2017*theta_fit*p_mean)/1000,
#                                 A_2017 = ggplot2::cut_number(round(A_2017),  n = 5, dig.lab = 5))) +
#   ggplot2::geom_sf(aes(fill = A_2017)) +
#   ggplot2::scale_fill_brewer(name = expression(paste("A"[2017]^"i", "(thousand USD)")), palette = "YlOrRd") +
#   ggplot2::geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
#   ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
#   ggplot2::theme(panel.grid.major = element_line(colour = "white"),
#                  panel.grid.minor = element_line(colour = "white"),
#                  panel.background = element_blank(),
#                  strip.background = element_rect(fill = NA),
#                  axis.line = element_blank(), axis.ticks = element_blank(),
#                  axis.title = element_blank(), axis.text = element_blank(),
#                  legend.title = element_text(hjust = 0.5, size = 22, face = "bold"),
#                  legend.position = "bottom", legend.key.width = unit(3, "cm"), legend.margin=margin(t=-0.5, r=0, b=0.2, l=0, unit="cm"),
#                  legend.text = element_text(size = 18, face = "bold"))

# ggplot2::ggsave(filename = here::here("plots/calibration/1043SitesModel/map_A2017_1043Sites.png"), plot= plot, width = 8, height = 6)

# x_2017
plot <- ggplot2::ggplot(data = calib_df %>%
  dplyr::mutate(x_2017 = ggplot2::cut_number(round(x_2017 / 1000000), n = 5, dig.lab = 3))) +
  ggplot2::geom_sf(aes(fill = x_2017)) +
  ggplot2::scale_fill_brewer(name = expression(paste("x"[2017]^"i", ~"(CO2e million Mg)")), palette = "YlOrRd") +
  ggplot2::geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::theme(
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    legend.position = "bottom", legend.key.width = unit(3, "cm"), legend.margin = margin(t = -0.5, r = 0, b = 0.2, l = 0, unit = "cm"),
    legend.text = element_text(size = 20, face = "bold")
  )

ggplot2::ggsave(filename = here::here("plots/calibration/1043SitesModel/map_x2017_1043Sites.png"), plot = plot, width = 8, height = 6)

# gamma_fit
# Define breaks and labels
breaks <- c(200, 380, 425, 530, 730, 1100)
labels <- c("200-380", "380-425", "425-530", "530-730", "730-1100")

# Modify your plotting code
plot <- ggplot2::ggplot(data = calib_df) +
  ggplot2::geom_sf(aes(fill = cut(gamma_fit, breaks = breaks, labels = labels))) +
  ggplot2::scale_fill_brewer(
    name = expression(paste(gamma^"i", ~"(CO2e Mg/ha)")),
    palette = "YlOrRd",
    breaks = labels
  ) + # Using the labels directly here, since they represent the bins now.
  ggplot2::geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::theme(
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    legend.position = "bottom", legend.key.width = unit(2, "cm"),
    legend.margin = margin(t = -0.5, r = 0, b = 0.2, l = 0, unit = "cm"),
    legend.text = element_text(size = 20, face = "bold")
  )

ggplot2::ggsave(filename = here::here("plots/calibration/1043SitesModel/map_gamma_fit.png"), plot = plot, width = 8, height = 6)

# theta_fit
plot <- ggplot2::ggplot(data = calib_df %>%
  dplyr::mutate(theta_fit = ggplot2::cut_number(round(theta_fit, 1), n = 5, dig.lab = 2))) +
  ggplot2::geom_sf(aes(fill = theta_fit)) +
  ggplot2::scale_fill_brewer(name = expression(paste(theta^"i")), palette = "YlOrRd") +
  ggplot2::geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::theme(
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    legend.position = "bottom", legend.key.width = unit(3, "cm"), legend.margin = margin(t = -0.5, r = 0, b = 0.2, l = 0, unit = "cm"),
    legend.text = element_text(size = 20, face = "bold")
  )

ggplot2::ggsave(filename = here::here("plots/calibration/1043SitesModel/map_theta_fit.png"), plot = plot, width = 8, height = 6)

# END TIMER
tictoc::toc(log = T)
