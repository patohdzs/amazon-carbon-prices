import argparse
import os
import time

import numpy as np
import pandas as pd
from cmdstanpy import CmdStanModel
import pickle
from pysrc.services.data_service import load_gamma_calib, load_theta_calib
from pysrc.services.file_service import get_path

# Define script arguments
parser = argparse.ArgumentParser(description="baseline sampling")
parser.add_argument("--sites", type=int, default=1043)
args = parser.parse_args()
num_sites = args.sites

# Define output paths
output_dir = get_path("data", "calibration")
os.makedirs(get_path("output", "tables"), exist_ok=True)

# Compile Stan code
sampler = CmdStanModel(
    stan_file=get_path("stan_model") / "baseline.stan",
    cpp_options={"STAN_THREADS": "true"},
    force_compile=True,
)

# Collect input data for Stan
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

# Sampling from baseline distribution
start_time = time.time()
fit = sampler.sample(
    data=data,
    **stan_kwargs,
)
end_time = time.time()
print(f"Finished sampling! Time: {end_time - start_time:.4f} seconds")
print(fit.diagnose())

# Compute average baseline gamma and theta
gamma_mean = fit.stan_variable("gamma").mean(axis=0)
theta_mean = fit.stan_variable("theta").mean(axis=0)

# Save productivity params
df = pd.DataFrame({"gamma_fit": gamma_mean, "theta_fit": theta_mean})

# Save to CSV
df.to_csv(output_dir / f"productivity_params_{num_sites}.csv", index=False)



gamma_samples = fit.stan_variable("gamma")  # shape: (draws, ...)
theta_samples = fit.stan_variable("theta")


output_path=os.path.join(
    str(get_path("output")),
    "sampling",
    "gurobi",
    f"{num_sites}sites",
    "baseline",
)
if not os.path.exists(output_path):
    os.makedirs(output_path)
outfile_path = output_path + "/results.pcl"
with open(outfile_path, "wb") as outfile:
    pickle.dump({"gamma": gamma_samples, "theta": theta_samples}, outfile)
    print(f"Results saved to {outfile_path}")



beta_gamma = fit.stan_variable("beta_gamma")
beta_theta = fit.stan_variable("beta_theta")
nu_gamma = fit.stan_variable("nu_gamma")
nu_theta = fit.stan_variable("nu_theta")

sigma_u_gamma = fit.stan_variable("sigma_u_gamma")
sigma_u_theta = fit.stan_variable("sigma_u_theta")
sigma_v_gamma = fit.stan_variable("sigma_v_gamma")
sigma_v_theta = fit.stan_variable("sigma_v_theta")


percentiles_gamma = np.percentile(beta_gamma, [10, 50, 90], axis=0)
percentiles_theta = np.percentile(beta_theta, [10, 50, 90], axis=0)
percentiles_gamma_sigma_u = np.percentile(sigma_u_gamma, [10, 50, 90], axis=0)
percentiles_theta_sigma_u = np.percentile(sigma_u_theta, [10, 50, 90], axis=0)
percentiles_gamma_sigma_v = np.percentile(sigma_v_gamma, [10, 50, 90], axis=0)
percentiles_theta_sigma_v = np.percentile(sigma_v_theta, [10, 50, 90], axis=0)

gamma_quan_table = pd.DataFrame(
    {
        "10th_percentile": percentiles_gamma[0],
        "50th_percentile": percentiles_gamma[1],
        "90th_percentile": percentiles_gamma[2],
    }
)
gamma_quan_table.to_csv(
    get_path("output", "tables") / f"gamma_percentiles_{num_sites}.csv", index=False
)

theta_quan_table = pd.DataFrame(
    {
        "10th_percentile": percentiles_theta[0],
        "50th_percentile": percentiles_theta[1],
        "90th_percentile": percentiles_theta[2],
    }
)
theta_quan_table.to_csv(
    get_path("output", "tables") / f"theta_percentiles_{num_sites}.csv", index=False
)


sigma_quan_table = pd.DataFrame(
    {
        "gamma_sigma_u_10th": [percentiles_gamma_sigma_u[0]],
        "gamma_sigma_v_10th": [percentiles_gamma_sigma_v[0]],
        "theta_sigma_u_10th": [percentiles_theta_sigma_u[0]],
        "theta_sigma_v_10th": [percentiles_theta_sigma_v[0]],
        "gamma_sigma_u_50th": [percentiles_gamma_sigma_u[1]],
        "gamma_sigma_v_50th": [percentiles_gamma_sigma_v[1]],
        "theta_sigma_u_50th": [percentiles_theta_sigma_u[1]],
        "theta_sigma_v_50th": [percentiles_theta_sigma_v[1]],
        "gamma_sigma_u_90th": [percentiles_gamma_sigma_u[2]],
        "gamma_sigma_v_90th": [percentiles_gamma_sigma_v[2]],
        "theta_sigma_u_90th": [percentiles_theta_sigma_u[2]],
        "theta_sigma_v_90th": [percentiles_theta_sigma_v[2]],
    }
)

sigma_quan_table.to_csv(
    get_path("output", "tables") / f"sigma_percentiles_{num_sites}.csv", index=False
)


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
