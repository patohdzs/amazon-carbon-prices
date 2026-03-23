# Carbon capture curve: df1 with empirical shifted 5y right + theoretical alpha=0.046 (no t0)

library(tidyverse)
library(ggplot2)

load("data/calibration/carbon_capture_curves/combined_df.Rdata")

df <- combined_df %>%
  filter(!is.na(agb), !is.na(sec), !is.na(gamma)) %>%
  filter(sec > 0, sec <= 30) %>%
  mutate(
    ratio   = agb / gamma,
    age_bin = pmin(30L, pmax(1L, ceiling(sec)))
  )

rm(combined_df)
gc()

m_dummy <- lm(ratio ~ 0 + factor(age_bin), data = df)
coefs <- coef(m_dummy)
ses   <- summary(m_dummy)$coefficients[, "Std. Error"]

emp_df <- tibble(
  age_bin     = 1:30,
  mean_ratio  = as.numeric(coefs[paste0("factor(age_bin)", 1:30)]),
  se          = as.numeric(ses[paste0("factor(age_bin)", 1:30)]),
  lower       = mean_ratio - se,
  upper       = mean_ratio + se
)

# Empirical with no shift, +5, +7
emp_plot <- bind_rows(
  emp_df %>% mutate(age_x = age_bin, type = "Estimates (no shift)"),
  emp_df %>% mutate(age_x = age_bin + 5, type = "Estimates (shift +5)"),
  emp_df %>% mutate(age_x = age_bin + 7, type = "Estimates (shift +7)")
)

# Theoretical: 1 - exp(-alpha * t), no t0, alpha = 0.046
theo_046 <- tibble(
  age = 0:38,
  value = 1 - exp(-0.046 * age),
  type = "Theoretical (alpha=0.046)"
)

dir.create("output/figures/carbon_capture", recursive = TRUE, showWarnings = FALSE)

p_main <- ggplot() +
  geom_point(data = emp_plot, aes(x = age_x, y = mean_ratio, color = type), size = 1.8) +
  geom_line(data = emp_plot, aes(x = age_x, y = mean_ratio, color = type), linewidth = 0.6) +
  geom_line(
    data = theo_046,
    aes(x = age, y = value, color = type),
    linewidth = 1
  ) +
  labs(
    x = "Age of secondary vegetation (years)",
    y = expression(ratio == AGB/gamma),
    title = "Carbon capture curve"
  ) +
  scale_x_continuous(limits = c(0, 30), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  scale_color_manual(
    values = c(
      "Estimates (no shift)" = "black",
      "Estimates (shift +5)" = "purple",
      "Estimates (shift +7)" = "darkgreen",
      "Theoretical (alpha=0.046)" = "blue"
    ),
    breaks = c("Estimates (no shift)", "Estimates (shift +5)", "Estimates (shift +7)", "Theoretical (alpha=0.046)"),
    drop = FALSE
  )

ggsave("output/figures/carbon_capture/gamma_secondary_vegetation_shifted.png",
       p_main, width = 8.5, height = 6)

rm(df, m_dummy, emp_df, emp_plot, theo_046, p_main)
gc()
