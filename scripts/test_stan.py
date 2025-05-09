import time

import matplotlib.pyplot as plt
import pandas as pd
from cmdstanpy import CmdStanModel

from pysrc.services.data_service import load_gamma_calib, load_theta_calib
from pysrc.services.file_service import get_path

output_dir = get_path("data", "calibration")

# Compile Stan code
sampler = CmdStanModel(
    stan_file=get_path("stan_model") / "baseline.stan",
    cpp_options={"STAN_THREADS": "true"},
    force_compile=True,
)


# Organize input data for Stan
num_sites = 1043
data = dict(
    num_sites=num_sites,
    **load_gamma_calib(num_sites, "reg"),
    **load_gamma_calib(num_sites, "fit"),
    **load_theta_calib(num_sites, "reg"),
    **load_theta_calib(num_sites, "fit"),
)

# Set sampling params
stan_kwargs = dict(
    iter_sampling=1500,
    iter_warmup=500,
    show_progress=True,
    seed=1,
    inits=0.2,
    chains=4,
)


# Sampling from adjusted distribution
start_time = time.time()
fit = sampler.sample(
    data=data,
    **stan_kwargs,
)
end_time = time.time()
print(f"Finished sampling! Time: {end_time - start_time:.4f} seconds")
print(fit.diagnose())


gamma_mean = fit.stan_variable("gamma").mean(axis=0)
theta_mean = fit.stan_variable("theta").mean(axis=0)

df = pd.DataFrame({"gamma_fit": gamma_mean, "theta_fit": theta_mean})

# Save to CSV
df.to_csv(output_dir / f"productivity_params_{num_sites}.csv", index=False)


## Test the model

beta_gamma = fit.stan_variable("beta_gamma")
beta_theta = fit.stan_variable("beta_theta")
nu_gamma = fit.stan_variable("nu_gamma")
nu_theta = fit.stan_variable("nu_theta")

sigma_u_gamma = fit.stan_variable("sigma_u_gamma")
sigma_u_theta = fit.stan_variable("sigma_u_theta")
sigma_v_gamma = fit.stan_variable("sigma_v_gamma")
sigma_v_theta = fit.stan_variable("sigma_v_theta")


df_beta_gamma = pd.DataFrame(
    beta_gamma, columns=[f"beta_gamma_{i}" for i in range(beta_gamma.shape[1])]
)
df_beta_theta = pd.DataFrame(
    beta_theta, columns=[f"beta_theta_{i}" for i in range(beta_theta.shape[1])]
)
df_nu_gamma = pd.DataFrame(
    nu_gamma, columns=[f"nu_gamma_{i}" for i in range(nu_gamma.shape[1])]
)
df_nu_theta = pd.DataFrame(
    nu_theta, columns=[f"nu_theta_{i}" for i in range(nu_theta.shape[1])]
)
df_sigma_u_gamma = pd.DataFrame({"sigma_u_gamma": sigma_u_gamma})
df_sigma_u_theta = pd.DataFrame({"sigma_u_theta": sigma_u_theta})
df_sigma_v_gamma = pd.DataFrame({"sigma_v_gamma": sigma_v_gamma})
df_sigma_v_theta = pd.DataFrame({"sigma_v_theta": sigma_v_theta})

df_all = pd.concat(
    [
        df_beta_gamma,
        df_beta_theta,
        df_nu_gamma,
        df_nu_theta,
        df_sigma_u_gamma,
        df_sigma_u_theta,
        df_sigma_v_gamma,
        df_sigma_v_theta,
    ],
    axis=1,
)

df_all.to_csv(
    output_dir / f"distribution_parameters_all_{num_sites}.csv",
    index=False,
)


plt.figure(figsize=(8, 5))
plt.hist(fit.stan_variable("sigma_u_gamma"), bins=50)  # 30 bins, nice edges
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.title(r"Histogram of $\sigma_u$")
# plt.grid(True)
plt.tight_layout()
plt.savefig(f"results/{num_sites}/histogram_sigma_u_gamma.png", dpi=300)
plt.show()


plt.figure(figsize=(8, 5))
plt.hist(fit.stan_variable("sigma_v_gamma"), bins=50)  # 30 bins, nice edges
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.title(r"Histogram of $\sigma_nu$")
# plt.grid(True)
plt.tight_layout()
plt.savefig(f"results/{num_sites}/histogram_sigma_nu_gamma.png", dpi=300)
plt.show()


plt.figure(figsize=(8, 5))
plt.hist(1 / fit.stan_variable("sigma_u_gamma") ** 2, bins=50)  # 30 bins, nice edges
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.title(r"Histogram of $\eta$")
# plt.grid(True)
plt.tight_layout()
plt.savefig(f"results/{num_sites}/histogram_eta_gamma.png", dpi=300)
plt.show()


plt.figure(figsize=(8, 5))
plt.hist(1 / fit.stan_variable("sigma_v_gamma") ** 2, bins=50)  # 30 bins, nice edges
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.title(r"Histogram of $\zeta$")
# plt.grid(True)
plt.tight_layout()
plt.savefig(f"results/{num_sites}/histogram_zeta_gamma.png", dpi=300)
plt.show()
