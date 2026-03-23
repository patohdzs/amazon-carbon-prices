# Same as 02_analysis.R but using combined_df_2 (Brazil raster)

library(tidyverse)
library(ggplot2)

# ---- Load data (df2) ----
load("data/calibration/carbon_capture_curves/combined_df_2.Rdata")

df <- combined_df_2 %>%
  filter(!is.na(agb), !is.na(sec), !is.na(gamma)) %>%
  filter(sec > 0, sec <= 30) %>%
  mutate(
    ratio   = agb / gamma,
    t       = sec,
    age_bin = pmin(30L, pmax(1L, ceiling(sec)))
  )

rm(combined_df_2)
gc()

# ---- Empirical ratio evolution (for plot) ----
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

# ---- Parametric fit on binned data (mean ratio by integer age) ----
alpha0 <- 0.06
t0_0   <- 5
alpha0 <- max(1e-6, alpha0)

m_nls <- nls(
  mean_ratio ~ 1 - exp(-alpha * (t0 + age_bin)),
  data      = emp_df,
  start     = list(alpha = alpha0, t0 = t0_0),
  algorithm = "port",
  lower     = c(alpha = 1e-8, t0 = 0),
  control   = nls.control(maxiter = 2000, warnOnly = TRUE)
)

alpha_hat <- coef(m_nls)[["alpha"]]
t0_hat    <- coef(m_nls)[["t0"]]
beta_hat  <- 1

cat(sprintf("\nNLS estimates (df2, beta=1):\nalpha_hat = %.5f\nt0_hat = %.2f\n",
            alpha_hat, t0_hat))

# ---- Curves for combined plot ----
nls_curve <- tibble(
  age = 0:30,
  value = beta_hat * (1 - exp(-alpha_hat * (t0_hat + age))),
  type = sprintf("NLS (alpha=%.3f, t0=%.1f)", alpha_hat, t0_hat)
)
theo_045 <- tibble(age = 0:30, value = 1 - exp(-0.045 * (0:30)),
  type = "Theoretical (alpha=0.045)")
curve_all <- bind_rows(nls_curve, theo_045)
emp_plot <- emp_df %>% mutate(type = "Estimates from data")

# ---- Plot ----
dir.create("output/figures/carbon_capture", recursive = TRUE, showWarnings = FALSE)

p_main <- ggplot() +
  geom_point(data = emp_plot, aes(x = age_bin, y = mean_ratio, color = type), size = 2) +
  geom_line(data = emp_plot, aes(x = age_bin, y = mean_ratio, color = type), linewidth = 0.7) +
  geom_errorbar(
    data = emp_plot,
    aes(x = age_bin, ymin = lower, ymax = upper),
    width = 0.2,
    alpha = 0.6,
    color = "black"
  ) +
  geom_line(
    data = curve_all,
    aes(x = age, y = value, color = type),
    linewidth = 1
  ) +
  labs(
    x = "Age of secondary vegetation (years)",
    y = expression(ratio == AGB/gamma),
    title = "Carbon capture curve (Brazil raster, df2)"
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom")

nls_label <- sprintf("NLS (alpha=%.3f, t0=%.1f)", alpha_hat, t0_hat)
p_main <- p_main +
  scale_color_manual(
    values = c(
      "Estimates from data" = "black",
      setNames("purple", nls_label),
      "Theoretical (alpha=0.045)" = "blue"
    ),
    breaks = c("Estimates from data", nls_label, "Theoretical (alpha=0.045)")
  )

ggsave("output/figures/carbon_capture/gamma_secondary_vegetation_df2.png",
       p_main, width = 8.5, height = 6)

rm(df, m_dummy, emp_df, emp_plot, m_nls, nls_curve, theo_045, curve_all, p_main)
gc()
