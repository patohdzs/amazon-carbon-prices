import geopandas as gpd
import numpy as np
import pandas as pd
from cmdstanpy import CmdStanModel

from pysrc.services.file_service import get_path

# Compile Stan code
sampler = CmdStanModel(
    stan_file=get_path("stan_model") / "baseline.stan",
    cpp_options={"STAN_THREADS": "true"},
    force_compile=True,
)


def load_gamma_calib(num_sites: int, type: str = "reg"):
    data_dir = get_path("data", "calibration")
    if type == "fit":
        # Used for gamma projection onto fitted values
        df = gpd.read_file(data_dir / f"gamma_fit_{num_sites}.geojson")

        # Get design matrix and its dimensions
        X = df.iloc[:, :6].to_numpy()
        N, K = X.shape

        # Large group indicator
        m = df["id_group"].astype(int)

        return {
            "X_gamma_fit": X,
            "m_gamma_fit": m,
        }
    else:
        # Used for gamma regression
        df = gpd.read_file(data_dir / "gamma_reg.geojson")

        # Get design matrix and its dimensions
        M = df["id_group"].unique().size
        y = df["log_co2e_ha_2017"]
        X = df.iloc[:, :6].to_numpy()
        N, K = X.shape

        # Large group indicator
        m = df["id_group"].astype(int)

        return {
            "N_gamma": N,
            "M_gamma": M,
            "K_gamma": K,
            "y_gamma": y,
            "X_gamma": X,
            "m_gamma": m,
        }


def load_theta_calib(num_sites: int, type: str = "reg"):
    data_dir = get_path("data", "calibration")
    if type == "fit":
        df = gpd.read_file(data_dir / f"theta_fit_{1043}.geojson")

        # Create the projection matrix
        G = (
            df.pivot(index="id", columns="muni_id", values="muni_site_area")
            .fillna(0)
            .to_numpy()
        )

        # Normalize to make row-stochastic
        G = G / G.sum(axis=1, keepdims=True)

        # Collapse data set to the municipality level
        df = df.sort_values("muni_id")

        # Keep first observation per municipality
        df = df.drop_duplicates(subset="muni_id", keep="first")

        # Get design matrix
        X = df.iloc[:, :8].to_numpy()
        C, _ = X.shape

        # Large group indicator
        m = df["group_id"].astype(int)

        # Cattle price in 2017
        pa_2017 = 44.9736197781184
        
        
        
        return {
            "C_theta_fit": C,
            "X_theta_fit": X,
            "m_theta_fit": m,
            "G_theta_fit": G,
            "pa_2017": pa_2017,
        }

    else:
        df = gpd.read_file(data_dir / "theta_reg.geojson")
        # Get number of groups
        M = df["group_id"].unique().size

        # Get design matrix and its dimensions
        y = df["log_slaughter"]
        X = df.iloc[:, :8].to_numpy()
        N, K = X.shape

        # Large group indicator
        m = df["group_id"].astype(int)

        W = df["weights"].values
        W= np.sqrt(W/np.std(W))

        return {
            "N_theta": N,
            "M_theta": M,
            "K_theta": K,
            "y_theta": y,
            "X_theta": X,
            "m_theta": m,
            "W_theta": W,
        }


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
    iter_sampling=5000,
    iter_warmup=500,
    show_progress=True,
    seed=1,
    inits=0.2,
)

# Sampling from adjusted distribution
fit = sampler.sample(
    data=data,
    **stan_kwargs,
)
print("Finished sampling!")
print(fit.diagnose())




gamma_fit_mean = fit.stan_variable("gamma").mean(axis=0)
theta_fit_mean = fit.stan_variable("theta").mean(axis=0)

df = pd.DataFrame(
    {
        "gamma_fit_mean": gamma_fit_mean,
        "theta_fit_mean": theta_fit_mean,
    }
)

# Save to CSV
df.to_csv(get_path("data","calibration") / "productivity_params.csv", index=False)



## Test the model

sigma_u_gamma=fit.stan_variable("sigma_u_gamma").mean(axis=0)
sigma_u_theta=fit.stan_variable("sigma_u_theta").mean(axis=0)
print("sigma_u_gamma",sigma_u_gamma)
print("sigma_u_theta",sigma_u_theta)
sigma_v_gamma=fit.stan_variable("sigma_v_gamma").mean(axis=0)
sigma_v_theta=fit.stan_variable("sigma_v_theta").mean(axis=0)
print("sigma_v_gamma",sigma_v_gamma)
print("sigma_v_theta",sigma_v_theta)


import matplotlib.pyplot as plt
plt.figure(figsize=(8, 5))
plt.hist(fit.stan_variable("sigma_u_gamma"), bins=50)  # 30 bins, nice edges
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title(r'Histogram of $\sigma_u$')
# plt.grid(True)
plt.tight_layout()
plt.savefig("test/histogram_sigma_u_gamma.png", dpi=300)
plt.show()


plt.figure(figsize=(8, 5))
plt.hist(fit.stan_variable("sigma_v_gamma"), bins=50)  # 30 bins, nice edges
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title(r'Histogram of $\sigma_nu$')
# plt.grid(True)
plt.tight_layout()
plt.savefig("test/histogram_sigma_nu_gamma.png", dpi=300)
plt.show()



plt.figure(figsize=(8, 5))
plt.hist(1/fit.stan_variable("sigma_u_gamma")**2, bins=50)  # 30 bins, nice edges
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title(r'Histogram of $\eta$')
# plt.grid(True)
plt.tight_layout()
plt.savefig("test/histogram_eta_gamma.png", dpi=300)
plt.show()


plt.figure(figsize=(8, 5))
plt.hist(1/fit.stan_variable("sigma_v_gamma")**2, bins=50)  # 30 bins, nice edges
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title(r'Histogram of $\zeta$')
# plt.grid(True)
plt.tight_layout()
plt.savefig("test/histogram_zeta_gamma.png", dpi=300)
plt.show()