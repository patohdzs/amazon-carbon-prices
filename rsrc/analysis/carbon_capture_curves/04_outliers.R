# Outliers: pixels with ratio = agb/gamma > 1

library(tidyverse)

# 1) Percentage of pixels with ratio > 1
load("data/calibration/carbon_capture_curves/combined_df.Rdata")
df <- combined_df %>%
  filter(!is.na(agb) & !is.na(sec) & !is.na(gamma), sec <= 30) %>%
  mutate(ratio = agb / gamma)

n_total <- nrow(df)
n_outlier <- sum(df$ratio > 1, na.rm = TRUE)
pct_outlier <- 100 * n_outlier / n_total
cat(sprintf("Pixels with ratio > 1: %s / %s (%.2f%%)\n", n_outlier, n_total, pct_outlier))

# 2) Scatterplot of x = -ln(1 - ratio) vs t (age)
ols_df <- df %>% filter(ratio < 1) %>% mutate(x = -log(1 - ratio))
m_no_int <- lm(x ~ sec - 1, data = ols_df)
m_int <- lm(x ~ sec, data = ols_df)
set.seed(42)
plot_df <- slice_sample(ols_df, n = min(5e4, nrow(ols_df)))
dir.create("output/figures/carbon_capture", recursive = TRUE, showWarnings = FALSE)
png("output/figures/carbon_capture/x_vs_t_scatterplot.png", width = 800, height = 600, res = 120)
plot(plot_df$sec, plot_df$x, pch = ".", col = rgb(0, 0, 0, 0.3),
  xlab = "t (age)", ylab = "x = -ln(1 - ratio)", main = "x vs t")
abline(coef(m_int)[1], coef(m_int)[2], col = "blue", lwd = 2)
abline(0, coef(m_no_int)["sec"], col = "red", lwd = 2)
legend("topleft", legend = c("With intercept", "No intercept"), col = c("blue", "red"), lwd = 2)
dev.off()

# 3) Percentage of each age cohort with ratio > 1
age_cohort <- pmin(30, pmax(1, ceiling(df$sec)))
cohort_stats <- df %>%
  mutate(age_cohort = age_cohort) %>%
  group_by(age_cohort) %>%
  summarise(
    n = n(),
    n_outlier = sum(ratio > 1, na.rm = TRUE),
    pct_outlier = 100 * n_outlier / n,
    .groups = "drop"
  )
cat("\n% ratio > 1 by age cohort:\n")
print(cohort_stats)
out_path <- "output/outliers_by_cohort.csv"
dir.create("output", showWarnings = FALSE)
write.csv(cohort_stats, out_path, row.names = FALSE)
cat(sprintf("\nSaved to %s\n", out_path))

# 4) Bar graph of % outliers by age cohort (stacked: blue = outlier %, grey = rest to 100%)
cohort_long <- cohort_stats %>%
  mutate(non_outlier = 100 - pct_outlier) %>%
  pivot_longer(c(pct_outlier, non_outlier), names_to = "type", values_to = "pct") %>%
  mutate(type = factor(type, levels = c("pct_outlier", "non_outlier"),
    labels = c("ratio > 1", "ratio ≤ 1")))
p_cohort <- ggplot(cohort_long, aes(x = age_cohort, y = pct, fill = type)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c("ratio ≤ 1" = "lightgrey", "ratio > 1" = "steelblue")) +
  labs(x = "Age cohort", y = expression("% " * W/gamma * " > 1"), title = "Outlier share by age cohort", fill = NULL) +
  theme_minimal() +
  theme(legend.position = "top")
ggsave("output/figures/carbon_capture/outliers_by_cohort.png", p_cohort, width = 8, height = 5)
