# Carbon capture curve: ratio = agb/gamma by age

library(tidyverse)
library(ggplot2)

load("data/calibration/carbon_capture_curves/combined_df.Rdata")

combined_df <- combined_df %>%
  filter(!is.na(agb) & !is.na(sec) & !is.na(gamma), sec <= 30) %>%
  mutate(ratio = agb / gamma)

generate_dummies <- function(df, sec_var = "sec", max_sec = 32) {
  for (i in 1:max_sec) {
    dummy_name <- paste0("dummy_", i)
    df <- df %>%
      mutate(!!dummy_name := ifelse(.data[[sec_var]] > (i - 1) & .data[[sec_var]] <= i, 1, 0))
  }
  return(df)
}
combined_df <- generate_dummies(combined_df, max_sec = 30)

no_intercept_formula <- as.formula(paste("ratio ~ -1 +", paste(paste0("dummy_", 1:30), collapse = " + ")))
no_intercept_model <- lm(no_intercept_formula, data = combined_df)
summary(no_intercept_model)

coef_no_intercept <- coef(no_intercept_model)[paste0("dummy_", 1:30)]
stderr_no_intercept <- summary(no_intercept_model)$coefficients[paste0("dummy_", 1:30), "Std. Error"]

no_intercept_coefficients_df <- data.frame(
  dummy = 1:30,
  coefficient = coef_no_intercept,
  lower_bound = coef_no_intercept - stderr_no_intercept,
  upper_bound = coef_no_intercept + stderr_no_intercept,
  type = "Estimates from data"
)

theoretical_df <- data.frame(
  dummy = 0:30,
  coefficient = 1 - exp(-0.045 * (0:30)),
  type = "Theoretical Function"
)

no_intercept_plot_df <- bind_rows(no_intercept_coefficients_df, theoretical_df)

p1 <- ggplot(no_intercept_plot_df, aes(x = dummy, y = coefficient, color = type, linetype = type)) +
  geom_point(data = no_intercept_coefficients_df) +
  geom_line(data = no_intercept_coefficients_df) +
  geom_errorbar(
    data = no_intercept_coefficients_df,
    aes(ymin = lower_bound, ymax = upper_bound),
    width = 0.2
  ) +
  geom_line(data = theoretical_df, linewidth = 1) +
  labs(x = "Age of Secondary Forest", y = "Percentage") +
  scale_color_manual(values = c("Estimates from data" = "black", "Theoretical Function" = "blue")) +
  scale_linetype_manual(values = c("Estimates from data" = "solid", "Theoretical Function" = "dashed")) +
  guides(color = guide_legend(title = NULL), linetype = guide_legend(title = NULL)) +
  theme_minimal() +
  theme(legend.position.inside = c(0.85, 0.15))

ggsave("output/figures/carbon_capture/gamma_secondary_vegetation.png", plot = p1, width = 8, height = 6)
