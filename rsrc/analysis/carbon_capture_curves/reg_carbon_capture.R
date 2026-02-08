library(tidyverse)
library(stargazer)
library(ggplot2)
library(broom)
library(fixest)
library(dplyr)

load("data/calibration/carbon_accumulation.Rdata")

# Fit regression on age dummies
model_1 <- lm(ratio ~ 0 + factor(age), data = df)

# Fit regression interacting with last pasture quality
model_2 <- df %>%
  filter(last_pq != 0) %>%
  lm(ratio ~ 0 + factor(age):factor(last_pq), data = .)

stargazer(model_1, model_2, out = "plots/carbon_capture/table_1.tex")
stargazer(model_1, model_2, out = "plots/carbon_capture/table_1.txt")

# Fit regression on theoretical law of motion
df <- df %>%
  mutate(Tp = 1 - exp(-0.045 * age))

model_2 <- lm(ratio ~ 0 + Tp, data = df)
model_3 <- lm(ratio ~ Tp, data = df)
stargazer(model_2, model_3, out = "plots/carbon_capture/table_2.tex")
stargazer(model_2, model_3, out = "plots/carbon_capture/table_2.txt")

# Fit regression interacting with precipitation and GHI
model_4 <- df %>%
  lm(ratio ~ 0 + percip + factor(age), data = .)

# Fit regression with GHI levels
model_5 <- df %>%
  lm(ratio ~ 0 + ghi + factor(age), data = .)

# Fit regression interacting with both GHI and mean percipitation
model_6 <- df %>%
  lm(ratio ~ 0 + percip + ghi + factor(age), data = .)

# Fit the above on logs
model_7 <- df %>%
  filter(!(ratio == 0)) %>%
  lm(log(ratio) ~ 0 + percip + ghi + factor(age), data = .)

stargazer(model_4, model_5, model_6, model_7, out = "plots/carbon_capture/table_4.tex")
stargazer(model_4, model_5, model_6, model_7, out = "plots/carbon_capture/table_4.txt")

# Drop rows with NAs
df <- df %>%
  filter(!is.na(ghi) & !is.na(percip) & !is.na(gamma))

# Regress gamma on ghi and percip
first_stage <- df %>%
  lm(log(gamma) ~ ghi + percip, data = .)

# Get residuals
df <- df %>%
  mutate(res = residuals(first_stage))

# Run full regression model with orthogonalized ratio
full_model_res <- df %>%
  filter(!(co2e == 0)) %>%
  lm(log(co2e) - res ~ percip + ghi + factor(age), data = .)

stargazer(first_stage, full_model_res, out = "plots/carbon_capture/table_5.txt")
stargazer(first_stage, full_model_res, out = "plots/carbon_capture/table_5.tex")

# correlation matrix between ghi and x and y
cor_matrix = cor(df[, c('ghi','x','y')])
write.table(cor_matrix, file = "plots/gamma_alpha/correlation_matrix_ghi_coordinates.txt", sep = "\t", col.names = NA)
latex_table <- xtable(cor_matrix)
print(latex_table, file = "plots/gamma_alpha/correlation_matrix.tex", type = "latex")

# Regression on only GHI without precip

log_gamma_from_ghi <- df %>%
  mutate(log_gamma = log(gamma)) %>%
  lm(log_gamma ~ ghi , data = .)

df <- df %>%
  mutate(
    residuals_model = residuals(log_gamma_from_ghi)
  )

df <- df %>%
mutate(log_co2e = log(co2e))%>%
  mutate(ratio_co2e_residuals = log_co2e - residuals_model)

# eliminate zero and infinity in ratio_co2e_residuals
df<- df[!is.na(df$ratio_co2e_residuals) & is.finite(df$ratio_co2e_residuals), ]

log_model_residual_ghi_age <- lm(ratio_co2e_residuals ~ 0 + ghi + factor(age), data = df)

stargazer(log_gamma_from_ghi, out = "plots/gamma_alpha/table_ghi.tex")
stargazer(log_gamma_from_ghi, out = "plots/gamma_alpha/table_ghi.txt")

stargazer(log_model_residual_ghi_age, out = "plots/gamma_alpha/table_co2e_ghi.tex")
stargazer(log_model_residual_ghi_age, out = "plots/gamma_alpha/table_co2e_ghi.txt")
