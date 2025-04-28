import geopandas as gpd
import numpy as np
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
        X = df.iloc[:, 1:6].to_numpy()
        N, K = X.shape

        # Large group indicator
        m = df["id_group"].astype(int)

        return {
            "X_gamma_fit": X,
            "m_gamma_fit": m,
        }
    else:
        # Used for gamma regression
        df = gpd.read_file(data_dir / f"gamma_reg_{num_sites}.geojson")

        # Get design matrix and its dimensions
        M = df["id_group"].unique().size
        y = df["log_co2e_ha_2017"]
        X = df.iloc[:, 1:6].to_numpy()
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

        # Get design matrix
        X = df.iloc[:, 1:8].to_numpy()
        C, _ = X.shape

        # Large group indicator
        m = df["group_id"].astype(int)

        # Municipal level to site level projection matrix
        G = np.array(
            [(df["id"].to_numpy() == i).astype(int) for i in range(1, num_sites + 1)]
        )

        # Multiply by area overalp weights
        G = df["muni_site_area"].to_numpy() * G
        G = G / G.sum(axis=1, keepdims=True)

        # Cattle price in 2017
        pa_2017 = 44.9736197781184

        return {
            "C_theta": C,
            "X_theta_fit": X,
            "m_theta_fit": m,
            "G_theta_fit": G,
            "pa_2017": pa_2017,
        }

    else:
        df = gpd.read_file(data_dir / f"theta_reg_{num_sites}.geojson")
        # Get number of groups
        M = df["group_id"].unique().size

        # Get design matrix and its dimensions
        y = df["log_slaughter"]
        X = df.iloc[:, 1:8].to_numpy()
        N, K = X.shape

        # Large group indicator
        m = df["group_id"].astype(int)

        return {
            "N_theta": N,
            "M_theta": M,
            "K_theta": K,
            "y_theta": y,
            "X_theta": X,
            "m_theta": m,
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
    iter_sampling=10000,
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


fit.stan_variable("gamma").mean(axis=0)
