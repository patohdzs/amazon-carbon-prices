# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: MAPS OF 1043 SITES MODEL PREDICTED VALUES
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

# SETUP
library(scales)
library(sf)
library(ggplot2)
library(readr)
library(dplyr)

# START TIMER
tictoc::tic(msg = "mapsPrediction_1043SitesModel.R script", log = T)

# OPTIONS
options(scipen = 999)

# DATA INPUT

# extract high and low price
p_high <- 41.11


# 1043 SITES MODEL CALIBRATION VARIABLES
load(here::here("data/calibration", "calibration_1043_sites.Rdata"))


# 1043 SITES AGGREGATE PREDICTION
aux.prices <- c(4.7, 14.7, 19.7, 24.7, 29.7)





prediction.1043SitesModel <- purrr::map_df(
  .x = aux.prices,
  .f = ~ {
    file_path <- here::here("output/optimization/hmc/gurobi/1043sites/xi_1.0/pa_41.11",
                            paste0("pe_", .x, "/Z.txt"))
    
    # Read the .txt file directly
    data <- readr::read_delim(file_path, 
                              delim = ",", 
                              col_names = FALSE,  # No column names
                              col_types = cols(.default = "n"))  # All numeric
    
    # Add a time column (T/R) ranging from 0 to 201
    time_period <- 0:(nrow(data) - 1)  # Assuming 202 time steps
    
    # Rename the columns to numeric ids (1, 2, ..., 1043)
    colnames(data) <- as.character(1:ncol(data))
    
    # Add the additional variables (T/R, p_e and p_a)
    data %>%
      dplyr::mutate(`T/R` = time_period, 
                    p_e = .x, 
                    p_a = p_high) %>%
      dplyr::relocate(`T/R`, .before = 1)  # Place T/R as the first column
  }
)



# AMAZON BIOME VECTOR DATA
load(here::here("data/clean/amazon_biome.Rdata"))


calib_df <- calib_df %>%
  st_transform(,crs=4326)

amazon_biome <- amazon_biome %>%
  st_transform(,crs=4326)


relative_entropy <- read.csv(here::here("output/figures/entropy/site_1043/xi1.0", "kl_divergences_theta_gamma.csv"))%>%
  mutate(id = row_number())


calib_df <- calib_df %>%
  left_join(relative_entropy, by = "id")

  

# select relevant calibrated variables
calib_df <-
  calib_df %>%
  dplyr::select(id, z_2017, zbar_2017,theta_b0,theta_b15,gamma_b0,gamma_b15)



# adjust predicted data to site-year panel + add calibrated variables
prediction.1043SitesModel <-
  prediction.1043SitesModel %>%
  tidyr::pivot_longer(cols = -c("T/R", "p_a", "p_e"), names_to = "id", values_to = "z_t") %>%
  dplyr::mutate(id = as.numeric(stringr::str_trim(id))) %>%
  dplyr::rename(time = "T/R") %>%
  dplyr::mutate(across(.cols = everything(), .fns = as.numeric)) %>%
  dplyr::right_join(calib_df, by = c("id" = "id")) %>% # match by id, guarantee that prediction data uses the same id than the calibrated data
  sf::st_as_sf() %>%
  dplyr::mutate(
    z_t =  1e11 * z_t / zbar_2017,
    z_2017_1043Sites =  100* z_2017 / zbar_2017
  ) # transform share to %

# clean environment
rm(calib_df)

# adjust projection
prediction.1043SitesModel <- sf::st_as_sf(prediction.1043SitesModel)
prediction.1043SitesModel <- sf::st_transform(prediction.1043SitesModel, sf::st_crs(4326))
amazon_biome <- sf::st_transform(amazon_biome, sf::st_crs(prediction.1043SitesModel))







# GENERATE MAPS


circle_theta_b0 <- prediction.1043SitesModel %>%
  dplyr::filter(id %in% c(986)) %>%
  sf::st_centroid() %>%
  dplyr::mutate(
    x = sf::st_coordinates(.)[,1],  # Extract x coordinate
    y = sf::st_coordinates(.)[,2]   # Extract y coordinate
  )

circle_theta_b15 <- prediction.1043SitesModel %>%
  dplyr::filter(id %in% c(1040)) %>%
  sf::st_centroid() %>%
  dplyr::mutate(
    x = sf::st_coordinates(.)[,1],  # Extract x coordinate
    y = sf::st_coordinates(.)[,2]   # Extract y coordinate
  )

circle_gamma_b0 <- prediction.1043SitesModel %>%
  dplyr::filter(id %in% c(691)) %>%
  sf::st_centroid() %>%
  dplyr::mutate(
    x = sf::st_coordinates(.)[,1],  # Extract x coordinate
    y = sf::st_coordinates(.)[,2]   # Extract y coordinate
  )

circle_gamma_b15 <- prediction.1043SitesModel %>%
  dplyr::filter(id %in% c(904)) %>%
  sf::st_centroid() %>%
  dplyr::mutate(
    x = sf::st_coordinates(.)[,1],  # Extract x coordinate
    y = sf::st_coordinates(.)[,2]   # Extract y coordinate
  )



# relative entropy b0
plot_theta_b0 <- prediction.1043SitesModel %>%
  filter(time == 0, p_e == aux.prices[3]) %>%
  arrange(theta_b0) %>%
  mutate(
    rank = row_number(),
    bin = ntile(rank, 6)
  ) %>%
  group_by(bin) %>%
  mutate(
    bin_label = paste0(min(rank), "–", max(rank))
  ) %>%
  ungroup() %>%
  mutate(
    bin_label = factor(bin_label, levels = unique(bin_label))  # ensure correct order
  ) %>%
  ggplot() +
  geom_sf(aes(fill = bin_label)) +
  scale_fill_brewer(name = NULL, palette = "YlOrRd", drop = FALSE) +
  geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
  geom_point(data = circle_theta_b0, aes(x = x, y = y),
             color = "blue", size = 30, shape = 21, fill = NA, stroke = 4) +
  guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  theme(
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
    legend.position = "bottom",
    legend.margin = margin(t = -1, r = 0, b = 0.3, l = 0, unit = "cm"),
    legend.text = element_text(size = 50, face = "bold"),
    plot.margin = unit(c(-0.5, -1, 0, -1), "cm")
  )



plot_theta_b15 <- prediction.1043SitesModel %>%
  filter(time == 0, p_e == aux.prices[3]) %>%
  arrange(theta_b15) %>%
  mutate(
    rank = row_number(),
    bin = ntile(rank, 6)
  ) %>%
  group_by(bin) %>%
  mutate(
    bin_label = paste0(min(rank), "–", max(rank))
  ) %>%
  ungroup() %>%
  mutate(
    bin_label = factor(bin_label, levels = unique(bin_label))  # preserve order
  ) %>%
  ggplot() +
  geom_sf(aes(fill = bin_label)) +
  scale_fill_brewer(name = NULL, palette = "YlOrRd", drop = FALSE) +
  geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
  geom_point(data = circle_theta_b15, aes(x = x, y = y),
             color = "blue", size = 30, shape = 21, fill = NA, stroke = 4) +
  guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  theme(
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
    legend.position = "bottom",
    legend.margin = margin(t = -1, r = 0, b = 0.3, l = 0, unit = "cm"),
    legend.text = element_text(size = 50, face = "bold"),
    plot.margin = unit(c(-0.5, -1, 0, -1), "cm")
  )




plot_gamma_b0 <- prediction.1043SitesModel %>%
  filter(time == 0, p_e == aux.prices[3]) %>%
  arrange(gamma_b0) %>%
  mutate(
    rank = row_number(),
    bin = ntile(rank, 6)
  ) %>%
  group_by(bin) %>%
  mutate(
    bin_label = paste0(min(rank), "–", max(rank))
  ) %>%
  ungroup() %>%
  mutate(
    bin_label = factor(bin_label, levels = unique(bin_label))
  ) %>%
  ggplot() +
  geom_sf(aes(fill = bin_label)) +
  scale_fill_brewer(name = NULL, palette = "YlOrRd", drop = FALSE) +
  geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
  geom_point(data = circle_gamma_b0, aes(x = x, y = y),
             color = "blue", size = 30, shape = 21, fill = NA, stroke = 4) +
  guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  theme(
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
    legend.position = "bottom",
    legend.margin = margin(t = -1, r = 0, b = 0.3, l = 0, unit = "cm"),
    legend.text = element_text(size = 50, face = "bold"),
    plot.margin = unit(c(-0.5, -1, 0, -1), "cm")
  )

plot_gamma_b15 <- prediction.1043SitesModel %>%
  filter(time == 0, p_e == aux.prices[3]) %>%
  arrange(gamma_b15) %>%
  mutate(
    rank = row_number(),
    bin = ntile(rank, 6)
  ) %>%
  group_by(bin) %>%
  mutate(
    bin_label = paste0(min(rank), "–", max(rank))
  ) %>%
  ungroup() %>%
  mutate(
    bin_label = factor(bin_label, levels = unique(bin_label))
  ) %>%
  ggplot() +
  geom_sf(aes(fill = bin_label)) +
  scale_fill_brewer(name = NULL, palette = "YlOrRd", drop = FALSE) +
  geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
  geom_point(data = circle_gamma_b15, aes(x = x, y = y),
             color = "blue", size = 30, shape = 21, fill = NA, stroke = 4) +
  guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  theme(
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
    legend.position = "bottom",
    legend.margin = margin(t = -1, r = 0, b = 0.3, l = 0, unit = "cm"),
    legend.text = element_text(size = 50, face = "bold"),
    plot.margin = unit(c(-0.5, -1, 0, -1), "cm")
  )




dir.create(here::here("plots/1043-hmc"), recursive = TRUE, showWarnings = FALSE)

ggpubr::ggexport(
  plot = plot_theta_b0,   
  filename = here::here(glue::glue("plots/1043-hmc/re_theta_b0.png")),  
  width = 2400,   
  height = 1500   
)

ggpubr::ggexport(
  plot = plot_theta_b15,   
  filename = here::here(glue::glue("plots/1043-hmc/re_theta_b15.png")),  
  width = 2400,   
  height = 1500   
)

ggpubr::ggexport(
  plot = plot_gamma_b0,   
  filename = here::here(glue::glue("plots/1043-hmc/re_gamma_b0.png")),  
  width = 2400,   
  height = 1500   
)

ggpubr::ggexport(
  plot = plot_gamma_b15,   
  filename = here::here(glue::glue("plots/1043-hmc/re_gamma_b15.png")),  
  width = 2400,   
  height = 1500   
)


# END TIMER
tictoc::toc(log = T)
