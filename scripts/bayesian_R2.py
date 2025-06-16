

import pandas as pd

import numpy as np
import matplotlib.pyplot as plt
from pysrc.services.data_service import (
    load_productivity_params,
    # load_reg_data,
    load_site_data,
    load_gamma_calib,
    load_theta_calib,
)
from pysrc.services.file_service import get_path

num_sites=1043


baseline_distribution = pd.read_csv(get_path("data","calibration")/f'distribution_parameters_all_{num_sites}.csv')
beta_gamma_cols = [col for col in baseline_distribution.columns if col.startswith('beta_gamma')]
beta_theta_cols = [col for col in baseline_distribution.columns if col.startswith('beta_theta')]
nu_gamma_cols = [col for col in baseline_distribution.columns if col.startswith('nu_gamma')]
nu_theta_cols = [col for col in baseline_distribution.columns if col.startswith('nu_theta')]



beta_gamma = np.array(baseline_distribution[beta_gamma_cols])
beta_theta = np.array(baseline_distribution[beta_theta_cols])
nu_gamma = np.array(baseline_distribution[nu_gamma_cols])
nu_theta = np.array(baseline_distribution[nu_theta_cols])
sigma_u_gamma = np.array(baseline_distribution["sigma_u_gamma"])
sigma_u_theta = np.array(baseline_distribution["sigma_u_theta"])



gamma_df=load_gamma_calib(num_sites, "reg")
theta_df=load_theta_calib(num_sites, "reg")


gamma_R2_draws = []
for s in range(len(sigma_u_gamma)):
    beta_s = beta_gamma[s, :]    # (K,)
    nu_s = nu_gamma[s, :]        # (B,)
    sigma_s = sigma_u_gamma[s]               # scalar

    yhat_s = gamma_df["X_gamma"] @ beta_s.T + nu_s[gamma_df["m_gamma"]-1]
    residual_s = np.array(gamma_df["y_gamma"]) - yhat_s

    var_fit = np.var(yhat_s, ddof=1)
    var_res = np.var(residual_s, ddof=1)

    
    # var_fit = np.var(yhat_s, ddof=1)
    # var_res = sigma_s ** 2

    R2 = var_fit / (var_fit + var_res)
    gamma_R2_draws.append(R2)

gamma_R2_draws = np.array(gamma_R2_draws)

print("Gamma Bayesian R^2 summary:")
print(f"Mean:   {np.mean(gamma_R2_draws):.4f}")
print(f"Median: {np.median(gamma_R2_draws):.4f}")
print(f"Std:    {np.std(gamma_R2_draws):.4f}")

plt.hist(gamma_R2_draws, bins=30, edgecolor='black',density=True)
plt.title(r"Bayesian $R^2$ - log(CO2e)", fontsize=16)
plt.xlabel("$R^2$", fontsize=16)
plt.ylabel("density", fontsize=16)
plt.savefig(get_path("output", "figures") / f"bayesian_r2_gamma_{num_sites}.png")
plt.close()



theta_R2_draws = []
for s in range(len(sigma_u_theta)):
    beta_s = beta_theta[s, :]    # (K_theta,)
    nu_s = nu_theta[s, :]        # (B_theta,)
    sigma_s = sigma_u_theta[s]   # scalar

    
    
    yhat_s = theta_df["X_theta"] @ beta_s.T + nu_s[theta_df["m_theta"] - 1]
    residual_s = np.array(theta_df["y_theta"]) - yhat_s

    var_fit = np.var(yhat_s, ddof=1)
    var_res = np.var(residual_s, ddof=1)

    


    R2 = var_fit / (var_fit + var_res)
    theta_R2_draws.append(R2)

theta_R2_draws = np.array(theta_R2_draws)

# Summary
print("Theta Bayesian R^2 summary:")
print(f"Mean:   {np.mean(theta_R2_draws):.4f}")
print(f"Median: {np.median(theta_R2_draws):.4f}")
print(f"Std:    {np.std(theta_R2_draws):.4f}")

# Plot
plt.hist(theta_R2_draws, bins=30, edgecolor='black',density=True)
plt.title(r"Bayesian $R^2$ - log(Slaughter Value)", fontsize=16)
plt.xlabel("$R^2$", fontsize=16)
plt.ylabel("density", fontsize=16)
plt.savefig(get_path("output", "figures") / f"bayesian_r2_theta_{num_sites}.png")
plt.close()

