library(tidyverse)
library(stargazer)
library(ggplot2)
library(broom)
library(fixest)

load("data/calibration/carbon_accumulation.Rdata")

# Fit regression from approach 1
model_1 <- lm(ratio ~ 0 + factor(age), data = df)

# Fit regression interacting with last pasture quality
model_2 <- df %>%
  filter(last_pq != 0) %>%
  lm(ratio ~ 0 + factor(age):factor(last_pq), data = .)

stargazer(model_1, model_2, out = "plots/gamma_alpha/table_1.tex")
stargazer(model_1, model_2, out = "plots/gamma_alpha/table_1.txt")

# Get coefficients
coefs <- coef(model_1)[paste0("factor(age)", 2:32)]
coefs_df <- data.frame(
  time = 2:32,
  coef = coefs,
  type = "Coefficients"
)

# Generate the theoretical values from t=0 to 32
time <- 0:max(df$age)
theoretical_values <- 1 - exp(-0.045 * time)

# Create a data frame for plotting the theoretical function
theoretical_df <- data.frame(
  time = time,
  coef = theoretical_values,
  type = "Theoretical Function"
)

# Combine the two data frames
plot_df <- bind_rows(coefs_df, theoretical_df)

# Plot the coefs and overlay the theoretical function
p <- plot_df |> ggplot(aes(x = time, y = coef, color = type, linetype = type)) +
  geom_point(data = coefs_df) +
  geom_line(data = coefs_df) +
  geom_line(data = theoretical_df) +
  labs(
    x = "Time (years)",
    y = "Coefficient",
    color = "Legend",
    linetype = "Legend"
  ) +
  scale_color_manual(values = c("Coefficients" = "black", "Theoretical Function" = "blue")) +
  scale_linetype_manual(values = c("Coefficients" = "solid", "Theoretical Function" = "dashed"))

# Save the plot
ggsave("plots/gamma_alpha/coefs.png", plot = p)

# Get coefficients from model with interactions
coefs_df <- coef(model_2) |>
  as.data.frame() |>
  rownames_to_column("term") |>
  as_tibble() |>
  rename(coef = `coef(model_2)`) |>
  separate(term, into = c("factor_age", "factor_last_pq"), sep = ":") |>
  mutate(
    age = gsub("factor\\(age\\)", "", factor_age),
    last_pq = gsub("factor\\(last_pq\\)", "", factor_last_pq)
  )

# Add theoretical curve
coefs_df <- coefs_df |>
  select(age, last_pq, coef) |>
  mutate(time = as.numeric(age)) |>
  mutate(theory = 1 - exp(-0.045 * time))

p <- coefs_df |> ggplot() +
  geom_line(aes(x = time, y = coef, color = last_pq)) +
  geom_line(aes(x = time, y = theory, color = "Theory")) +
  geom_point(aes(x = time, y = coef, color = last_pq)) +
  geom_point(aes(x = time, y = theory, color = "Theory")) +
  labs(
    x = "Time (years)",
    y = "Coefficient",
    color = "Legend",
    linetype = "Legend"
  )

# Save the plot
ggsave("plots/gamma_alpha/coefs_by_pq.png", plot = p)


# Perform approach 2 regression
df <- df %>%
  mutate(Tp = 1 - exp(-0.045 * age))

model_2 <- lm(ratio ~ 0 + Tp, data = df)
model_3 <- lm(ratio ~ Tp, data = df)
model_4 <- lm(ratio ~ Tp + mp, data = df)
stargazer(model_2, model_3, out = "plots/gamma_alpha/table_2.tex")
stargazer(model_2, model_3, out = "plots/gamma_alpha/table_2.txt")
