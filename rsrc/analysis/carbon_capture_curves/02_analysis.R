# Carbon capture curve: ratio = AGB/gamma ~ 1 - exp(-alpha * t)

library(tidyverse)
library(ggplot2)

# ---- Load data ----
load("data/calibration/carbon_capture_curves/combined_df.Rdata")

df <- combined_df %>%
  filter(!is.na(agb), !is.na(sec), !is.na(gamma)) %>%
  filter(sec > 0, sec <= 30) %>%
  mutate(
    ratio   = agb / gamma,
    t       = sec,
    age_bin = pmin(30L, pmax(1L, ceiling(sec)))
  )

rm(combined_df)
gc()

# ---- Empirical ratio evolution (for plot) ----
m_dummy <- lm(ratio ~ 0 + factor(age_bin), data = df)
m_dummy_summary <- summary(m_dummy)

coefs <- coef(m_dummy)
ses   <- m_dummy_summary$coefficients[, "Std. Error"]

emp_df <- tibble(
  age_bin     = 1:30,
  mean_ratio  = as.numeric(coefs[paste0("factor(age_bin)", 1:30)]),
  se          = as.numeric(ses[paste0("factor(age_bin)", 1:30)]),
  lower       = mean_ratio - se,
  upper       = mean_ratio + se
)


# Theoretical curve
theo_045 <- tibble(age = 0:30, value = 1 - exp(-0.045 * (0:30)),
  type = "Theoretical (alpha=0.045)")

# Empirical: use age_bin 1:30, add age 0 with NA or omit
emp_plot <- emp_df %>% mutate(type = "Estimates from data")

# ---- Combined plot (gamma_secondary_vegetation) ----
dir.create("output/figures/carbon_capture", recursive = TRUE, showWarnings = FALSE)

p_main <- ggplot() +
  geom_line(data = emp_plot, aes(x = age_bin, y = mean_ratio, group = 1, color = "Estimates from data"), linewidth = 0.8) +
  geom_point(data = emp_plot, aes(x = age_bin, y = mean_ratio, color = "Estimates from data"), size = 2) +
  geom_line(
    data = theo_045,
    aes(x = age, y = value, color = "Theoretical Function"),
    linetype = "dashed",
    linewidth = 1
  ) +
  labs(
    x = "Age of secondary vegetation (years)",
    y = expression("Ratio of maximum carbon density" == X[p]/(z[p] * gamma[p])),
    title = "Carbon capture curve"
  ) +
  scale_color_manual(
    values = c("Estimates from data" = "black", "Theoretical Function" = "blue"),
    breaks = c("Estimates from data", "Theoretical Function")
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "right")

ggsave("output/figures/carbon_capture/gamma_secondary_vegetation.png",
       p_main, width = 8.5, height = 6)

cat("R2 =", m_dummy_summary$r.squared, "\n")
cat("Adjusted R2 =", m_dummy_summary$adj.r.squared, "\n")

# Cleanup
rm(df, m_dummy, m_dummy_summary, emp_df, emp_plot, theo_045, p_main)
gc()
