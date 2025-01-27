library(sf)
library(glue)
library(tidyverse)
library(ggplot2)


load("data/calibration/calibration_1043_sites.Rdata")

hist_bio_loss <- read_csv("results/bio/bio_hist.csv")

calib_df$pct_change_bio <- hist_bio_loss$pct_change_bio_SAR

calib_df <- calib_df %>%
  st_transform(4326)

pct_breaks <- quantile(
  calib_df$pct_change_bio,
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5),
  na.rm = TRUE
)

pct_breaks[[5]] <- 0.00

calib_df <- calib_df %>%
  mutate(
    percentile_change_bio = cut(
      pct_change_bio,
      breaks = c(-Inf, pct_breaks, Inf),
      labels = c(
        glue("< {round(pct_breaks[[1]], 2)}"),
        glue("[{round(pct_breaks[[1]], 2)}, {round(pct_breaks[[2]], 2)})"),
        glue("[{round(pct_breaks[[2]], 2)}, {round(pct_breaks[[3]], 2)})"),
        glue("[{round(pct_breaks[[3]], 2)}, {round(pct_breaks[[4]], 2)})"),
        glue("[{round(pct_breaks[[4]], 3)}, {round(pct_breaks[[5]], 2)})"),
        glue("> {round(pct_breaks[[5]], 2)}")
      ),
      include.lowest = TRUE
    )
  )


fig <- calib_df %>%
  ggplot() +
  geom_sf(aes(fill = percentile_change_bio)) +
  scale_fill_brewer(
    name = "Change in number of species (%)",
    palette = "YlOrRd"
  ) +
  guides(fill = guide_legend(
    label.position = "bottom",
    title.position = "top", nrow = 1
  )) +
  guides(
    fill = guide_legend(
      label.position = "bottom",
      title.position = "top",
      nrow = 1
    )
  ) +
  theme(
    legend.title = element_text(face = "bold"),
    legend.position = "bottom",
  )

ggsave(filename = "results/bio/hist_bio_loss.pdf", plot = fig)


fig <- calib_df %>%
  ggplot() +
  geom_sf(aes(fill = pct_change_bio)) +
  scale_fill_viridis_c(name = "% change in number of species")

ggsave(filename = "results/bio/hist_bio_loss_pct.pdf", plot = fig)
