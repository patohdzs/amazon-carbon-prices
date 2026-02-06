# ============================================================================
# Script: 02_analysis.R
# Purpose: Run regression analyses and create visualizations for AGB ratios
# ============================================================================

# Load libraries ----------------------------------------------------------
library(tidyverse)
library(ggplot2)

# Load data ---------------------------------------------------------------
load("data/calibration/combined_df.Rdata")

# Prepare data ------------------------------------------------------------

# Filter and create ratio variable
combined_df <- combined_df %>%
  filter(!is.na(agb) & !is.na(sec) & !is.na(gamma)) %>%
  mutate(ratio = agb / gamma)

# Create version with pasture quality filter (for potential future use)
combined_df_pq <- combined_df %>%
  filter(pq != 0)

# Helper function to generate dummy variables -----------------------------
generate_dummies <- function(df, sec_var = "sec", max_sec = 32) {
  for (i in 1:max_sec) {
    dummy_name <- paste0("dummy_", i)
    df <- df %>%
      mutate(!!dummy_name := ifelse(.data[[sec_var]] > (i - 1) & .data[[sec_var]] <= i, 1, 0))
  }
  return(df)
}

# Generate dummy variables
combined_df <- generate_dummies(combined_df, max_sec = 32)

# Model 1: Without intercept (coefficients for all dummies 1-32) ---------
cat("\n=== MODEL 1: WITHOUT INTERCEPT ===\n")

no_intercept_formula <- as.formula(paste("ratio ~ -1 +", paste(paste0("dummy_", 1:32), collapse = " + ")))
no_intercept_model <- lm(no_intercept_formula, data = combined_df)
print(summary(no_intercept_model))

# Extract coefficients for plotting
coef_no_intercept <- coef(no_intercept_model)[paste0("dummy_", 1:32)]
stderr_no_intercept <- summary(no_intercept_model)$coefficients[paste0("dummy_", 1:32), "Std. Error"]

no_intercept_coefficients_df <- data.frame(
  dummy = 1:32,
  coefficient = coef_no_intercept,
  lower_bound = coef_no_intercept - stderr_no_intercept,
  upper_bound = coef_no_intercept + stderr_no_intercept,
  type = "Coefficients"
)

# Model 2: With intercept (coefficients for dummies 2-32) ----------------
cat("\n=== MODEL 2: WITH INTERCEPT ===\n")

with_intercept_formula <- as.formula(paste("ratio ~", paste(paste0("dummy_", 2:32), collapse = " + ")))
with_intercept_model <- lm(with_intercept_formula, data = combined_df)
print(summary(with_intercept_model))

# Extract coefficients for plotting
coef_with_intercept <- coef(with_intercept_model)[paste0("dummy_", 2:32)]
stderr_with_intercept <- summary(with_intercept_model)$coefficients[paste0("dummy_", 2:32), "Std. Error"]

with_intercept_coefficients_df <- data.frame(
  dummy = 2:32,
  coefficient = coef_with_intercept,
  lower_bound = coef_with_intercept - stderr_with_intercept,
  upper_bound = coef_with_intercept + stderr_with_intercept,
  type = "Coefficients"
)

# Create theoretical function for comparison ------------------------------
x_values <- 0:32
theoretical_values <- 1 - exp(-0.045 * x_values)

theoretical_df <- data.frame(
  dummy = x_values,
  coefficient = theoretical_values,
  type = "Theoretical Function"
)

# Plot 1: Without intercept -----------------------------------------------
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
  labs(
    title = "Coefficients of Dummy Variables with Theoretical Function (No Intercept)",
    x = "Dummy Variable Index (age)",
    y = "Coefficient Value",
    color = "Legend",
    linetype = "Legend"
  ) +
  scale_color_manual(values = c("Coefficients" = "black", "Theoretical Function" = "blue")) +
  scale_linetype_manual(values = c("Coefficients" = "solid", "Theoretical Function" = "dashed")) +
  theme_minimal() +
  theme(legend.position.inside = c(0.85, 0.15))

ggsave("results/carbon_capture/coefficients_plot.png",
  plot = p1, width = 8, height = 6
)
cat("\nPlot saved: coefficients_plot.png\n")

# Plot 2: With intercept --------------------------------------------------
with_intercept_plot_df <- bind_rows(with_intercept_coefficients_df, theoretical_df)

p2 <- ggplot(with_intercept_plot_df, aes(x = dummy, y = coefficient, color = type, linetype = type)) +
  geom_point(data = with_intercept_coefficients_df) +
  geom_line(data = with_intercept_coefficients_df) +
  geom_errorbar(
    data = with_intercept_coefficients_df,
    aes(ymin = lower_bound, ymax = upper_bound),
    width = 0.2
  ) +
  geom_line(data = theoretical_df, linewidth = 1) +
  labs(
    title = "Coefficients of Dummy Variables with Theoretical Function (With Intercept)",
    x = "Dummy Variable Index (age)",
    y = "Coefficient Value",
    color = "Legend",
    linetype = "Legend"
  ) +
  scale_color_manual(values = c("Coefficients" = "black", "Theoretical Function" = "blue")) +
  scale_linetype_manual(values = c("Coefficients" = "solid", "Theoretical Function" = "dashed")) +
  theme_minimal() +
  theme(legend.position.inside = c(0.85, 0.15))

ggsave("results/carbon_capture/coefficients_plot_with_intercept.png",
  plot = p2, width = 8, height = 6
)
cat("Plot saved: coefficients_plot_with_intercept.png\n")

# Parametric models using theoretical function ---------------------------

# Add theoretical predictor variable
combined_df <- combined_df %>%
  mutate(Tp = 1 - exp(-0.045 * sec))

# Parametric model 1: With intercept
cat("\n=== PARAMETRIC MODEL 1: WITH INTERCEPT ===\n")
parametric_intercept_model <- lm(ratio ~ Tp, data = combined_df)
print(summary(parametric_intercept_model))

# Parametric model 2: Without intercept
cat("\n=== PARAMETRIC MODEL 2: WITHOUT INTERCEPT ===\n")
parametric_no_intercept_model <- lm(ratio ~ -1 + Tp, data = combined_df)
print(summary(parametric_no_intercept_model))

cat("\nAnalysis complete!\n")
