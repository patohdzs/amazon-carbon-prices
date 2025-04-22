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
load(here::here("data/calibration/hmc", "calibration_1043_sites.Rdata"))


# 1043 SITES AGGREGATE PREDICTION
aux.prices <- c(4.5, 14.5, 19.5, 24.5, 29.5)





prediction.1043SitesModel <- purrr::map_df(
  .x = aux.prices,
  .f = ~ {
    file_path <- here::here("output/optimization/hmc/gams/1043sites/xi_5.0/pa_41.11",
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


relative_entropy <- read.csv(here::here("output/figures/entropy/site_1043/xi5.0", "kl_divergences_theta_gamma.csv"))%>%
  mutate(id = row_number())

gamma_fit <- read.csv(here::here("data/calibration/hmc", "gamma_fit_1043.csv"))%>%
  mutate(id = row_number())
theta_fit <- read.csv(here::here("data/calibration/hmc", "theta_fit_1043.csv"))%>%
  mutate(id = row_number())
calib_df <- calib_df %>%
  left_join(gamma_fit, by = "id")%>%
  left_join(theta_fit, by="id")%>%
  left_join(relative_entropy, by = "id")

  
calib_df <- calib_df %>%
  mutate(
    x_1995 = calib_df$gamma_fit * (zbar_1995 - z_1995),
    x_2008 = calib_df$gamma_fit * (zbar_2008 - z_2008),
    x_2017 = calib_df$gamma_fit * (zbar_2017 - z_2017)
  )



# DATASET CLEANUP AND PREP

# select relevant calibrated variables
calib_df <-
  calib_df %>%
  dplyr::mutate(
    rank_theta_1043Sites = dense_rank(desc(theta_fit)), # rank values
    rank_gamma_1043Sites = dense_rank(desc(gamma_fit))
  ) %>% # rank values
  dplyr::select(id, z_2017, theta, rank_theta_1043Sites, rank_gamma_1043Sites, x_2017, zbar_2017,theta_b0,theta_b15,gamma_b0,gamma_b15)



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
  dplyr::filter(id %in% c(1021)) %>%
  sf::st_centroid() %>%
  dplyr::mutate(
    x = sf::st_coordinates(.)[,1],  # Extract x coordinate
    y = sf::st_coordinates(.)[,2]   # Extract y coordinate
  )

circle_gamma_b0 <- prediction.1043SitesModel %>%
  dplyr::filter(id %in% c(1015)) %>%
  sf::st_centroid() %>%
  dplyr::mutate(
    x = sf::st_coordinates(.)[,1],  # Extract x coordinate
    y = sf::st_coordinates(.)[,2]   # Extract y coordinate
  )

circle_gamma_b15 <- prediction.1043SitesModel %>%
  dplyr::filter(id %in% c(938)) %>%
  sf::st_centroid() %>%
  dplyr::mutate(
    x = sf::st_coordinates(.)[,1],  # Extract x coordinate
    y = sf::st_coordinates(.)[,2]   # Extract y coordinate
  )



# relative entropy b0
plot_theta_b0 <-
  ggplot2::ggplot(data = prediction.1043SitesModel %>%
                    dplyr::filter(time == 0, p_e == aux.prices[3]) %>%
                    dplyr::mutate(theta_b0 = cut(theta_b0,
                                            breaks = c(0, 0.1,1, 2, 4, 6, 12),
                                            include.lowest = T,
                                            dig.lab = 3,
                                            labels = c("[0~ 0.1", " ~ 1", " ~ 2", " ~ 4", " ~ 6", " ~ 12]")
                    ))) +
  ggplot2::geom_sf(aes(fill = theta_b0)) +
  ggplot2::scale_fill_manual(name = NULL, values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
  ggplot2::geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
  ggplot2::geom_point(data = circle_theta_b0, aes(x = x, y = y),
                      color = "blue", size = 30, shape = 21, fill = NA, stroke = 4) +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::theme(
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
    legend.position = "bottom", legend.margin = margin(t = -1, r = -0, b = 0.3, l = -0, unit = "cm"),
    legend.text = element_text(size = 50, face = "bold"),
    plot.margin = unit(c(t = -0.5, r = -1, b = -0, l = -1), "cm")
  )



plot_theta_b15 <-
  ggplot2::ggplot(data = prediction.1043SitesModel %>%
                    dplyr::filter(time == 0, p_e == aux.prices[3]) %>%
                    dplyr::mutate(theta_b15 = cut(theta_b15,
                                            breaks = c(0,0.01, 0.04, 0.08, 0.12,0.16,0.205),
                                            include.lowest = T,
                                            dig.lab = 3,
                                            labels = c("[0~ 0.01", " ~ 0.04", " ~ 0.08", " ~ 0.12", " ~ 0.16", " ~ 0.20]")
                    ))) +
  ggplot2::geom_sf(aes(fill = theta_b15)) +
  ggplot2::scale_fill_manual(name = NULL, values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
  ggplot2::geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
  ggplot2::geom_point(data = circle_theta_b15, aes(x = x, y = y),
                      color = "blue", size = 30, shape = 21, fill = NA, stroke = 4) +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::theme(
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
    legend.position = "bottom", legend.margin = margin(t = -1, r = -0, b = 0.3, l = -0, unit = "cm"),
    legend.text = element_text(size = 50, face = "bold"),
    plot.margin = unit(c(t = -0.5, r = -1, b = -0, l = -1), "cm")
  )



plot_gamma_b0 <-
  ggplot2::ggplot(data = prediction.1043SitesModel %>%
                    dplyr::filter(time == 0, p_e == aux.prices[3]) %>%
                    dplyr::mutate(gamma_b0 = cut(gamma_b0,
                                            breaks = c(0, 0.002,0.005, 0.01, 0.02, 0.03, 0.041),
                                            include.lowest = T,
                                            dig.lab = 3,
                                            labels = c("[0~ 0.002", " ~ 0.005", " ~ 0.01", " ~ 0.02", " ~ 0.03", " ~ 0.04]")
                    ))) +
  ggplot2::geom_sf(aes(fill = gamma_b0)) +
  ggplot2::scale_fill_manual(name = NULL, values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
  ggplot2::geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
  ggplot2::geom_point(data = circle_gamma_b0, aes(x = x, y = y),
                      color = "blue", size = 30, shape = 21, fill = NA, stroke = 4) +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::theme(
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
    legend.position = "bottom", legend.margin = margin(t = -1, r = -0, b = 0.3, l = -0, unit = "cm"),
    legend.text = element_text(size = 50, face = "bold"),
    plot.margin = unit(c(t = -0.5, r = -1, b = -0, l = -1), "cm")
  )


plot_gamma_b15 <-
  ggplot2::ggplot(data = prediction.1043SitesModel %>%
                    dplyr::filter(time == 0, p_e == aux.prices[3]) %>%
                    dplyr::mutate(gamma_b15 = cut(gamma_b15,
                                            breaks = c(0, 0.1,0.2, 0.3, 0.4, 0.5, 0.6),
                                            include.lowest = T,
                                            dig.lab = 3,
                                            labels = c("[0~ 0.1", " ~ 0.2", " ~ 0.3", " ~ 0.4", " ~ 0.5", " ~ 0.6]")
                    ))) +
  ggplot2::geom_sf(aes(fill = gamma_b15)) +
  ggplot2::scale_fill_manual(name = NULL, values = c("white", RColorBrewer::brewer.pal(5, "YlOrRd")), drop = FALSE) +
  ggplot2::geom_sf(data = amazon_biome, fill = NA, color = "darkgreen", size = 1.2) +
  ggplot2::geom_point(data = circle_gamma_b15, aes(x = x, y = y),
                      color = "blue", size = 30, shape = 21, fill = NA, stroke = 4) +
  ggplot2::guides(fill = guide_legend(label.position = "bottom", title.position = "top", nrow = 1)) +
  ggplot2::theme(
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white"),
    panel.background = element_blank(),
    strip.background = element_rect(fill = NA),
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank(),
    legend.title = element_text(hjust = 0.5, size = 40, face = "bold"),
    legend.position = "bottom", legend.margin = margin(t = -1, r = -0, b = 0.3, l = -0, unit = "cm"),
    legend.text = element_text(size = 50, face = "bold"),
    plot.margin = unit(c(t = -0.5, r = -1, b = -0, l = -1), "cm")
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
