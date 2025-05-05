import geopandas as gpd
import numpy as np
import pandas as pd
from cmdstanpy import CmdStanModel
from scipy.sparse import coo_matrix
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
        df = gpd.read_file(data_dir / f"theta_fit_{num_sites}.geojson")

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
        
        
        G_sparse = coo_matrix(G)
        row_G_theta = G_sparse.row + 1  # Stan uses 1-based indexing
        col_G_theta = G_sparse.col + 1
        val_G_theta = G_sparse.data
        N_nonzero_G_theta = len(val_G_theta)
        
        return {
            "C_theta_fit": C,
            "X_theta_fit": X,
            "m_theta_fit": m,
            # "G_theta_fit": G,
            "N_nonzero_G_theta": N_nonzero_G_theta,
            "row_G_theta": row_G_theta,
            "col_G_theta": col_G_theta,
            "val_G_theta": val_G_theta,
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
    iter_sampling=1500,
    iter_warmup=500,
    show_progress=True,
    seed=1, 
    inits=0.2,
    chains=8,
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
        "gamma_fit": gamma_fit_mean,
        "theta_fit": theta_fit_mean,
    }
)

# Save to CSV
df.to_csv(get_path("data","calibration") / f"productivity_params_{num_sites}.csv", index=False)



## Test the model

beta_gamma = fit.stan_variable("beta_gamma")
beta_theta = fit.stan_variable("beta_theta")
nu_gamma = fit.stan_variable("nu_gamma")
nu_theta = fit.stan_variable("nu_theta")

sigma_u_gamma=fit.stan_variable("sigma_u_gamma")
sigma_u_theta=fit.stan_variable("sigma_u_theta")
sigma_v_gamma=fit.stan_variable("sigma_v_gamma")
sigma_v_theta=fit.stan_variable("sigma_v_theta")


df_beta_gamma = pd.DataFrame(beta_gamma, columns=[f"beta_gamma_{i}" for i in range(beta_gamma.shape[1])])
df_beta_theta = pd.DataFrame(beta_theta, columns=[f"beta_theta_{i}" for i in range(beta_theta.shape[1])])
df_nu_gamma = pd.DataFrame(nu_gamma, columns=[f"nu_gamma_{i}" for i in range(nu_gamma.shape[1])])
df_nu_theta = pd.DataFrame(nu_theta, columns=[f"nu_theta_{i}" for i in range(nu_theta.shape[1])])
df_sigma_u_gamma = pd.DataFrame({'sigma_u_gamma': sigma_u_gamma})
df_sigma_u_theta = pd.DataFrame({'sigma_u_theta': sigma_u_theta})
df_sigma_v_gamma = pd.DataFrame({'sigma_v_gamma': sigma_v_gamma})
df_sigma_v_theta = pd.DataFrame({'sigma_v_theta': sigma_v_theta})

df_all = pd.concat([
    df_beta_gamma,
    df_beta_theta,
    df_nu_gamma,
    df_nu_theta,
    df_sigma_u_gamma,
    df_sigma_u_theta,
    df_sigma_v_gamma,
    df_sigma_v_theta
], axis=1)

df_all.to_csv(get_path("data", "calibration") / f"distribution_parameters_all_{num_sites}.csv", index=False)





import matplotlib.pyplot as plt
plt.figure(figsize=(8, 5))
plt.hist(fit.stan_variable("sigma_u_gamma"), bins=50)  # 30 bins, nice edges
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title(r'Histogram of $\sigma_u$')
# plt.grid(True)
plt.tight_layout()
plt.savefig(f"test/{num_sites}/histogram_sigma_u_gamma.png", dpi=300)
plt.show()


plt.figure(figsize=(8, 5))
plt.hist(fit.stan_variable("sigma_v_gamma"), bins=50)  # 30 bins, nice edges
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title(r'Histogram of $\sigma_nu$')
# plt.grid(True)
plt.tight_layout()
plt.savefig(f"test/{num_sites}/histogram_sigma_nu_gamma.png", dpi=300)
plt.show()



plt.figure(figsize=(8, 5))
plt.hist(1/fit.stan_variable("sigma_u_gamma")**2, bins=50)  # 30 bins, nice edges
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title(r'Histogram of $\eta$')
# plt.grid(True)
plt.tight_layout()
plt.savefig(f"test/{num_sites}/histogram_eta_gamma.png", dpi=300)
plt.show()


plt.figure(figsize=(8, 5))
plt.hist(1/fit.stan_variable("sigma_v_gamma")**2, bins=50)  # 30 bins, nice edges
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title(r'Histogram of $\zeta$')
# plt.grid(True)
plt.tight_layout()
plt.savefig(f"test/{num_sites}/histogram_zeta_gamma.png", dpi=300)
plt.show()