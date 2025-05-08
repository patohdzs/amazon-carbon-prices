import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from ..services.file_service import get_path


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

        G_sparse = csr_matrix(G)

        # Extract CSR components
        G_w = G_sparse.data  # Non-zero values
        G_v = G_sparse.indices + 1  # Column indices
        G_u = G_sparse.indptr + 1  # Row pointers
        G_nnz = len(G_w)

        return {
            "C_theta_fit": C,
            "X_theta_fit": X,
            "m_theta_fit": m,
            "G_nnz_theta": G_nnz,
            "G_w_theta": G_w,
            "G_v_theta": G_v,
            "G_u_theta": G_u,
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
        W = np.sqrt(W / np.std(W))

        return {
            "N_theta": N,
            "M_theta": M,
            "K_theta": K,
            "y_theta": y,
            "X_theta": X,
            "m_theta": m,
            "W_theta": W,
        }


def load_site_data(num_sites: int, year: int = 2017, norm_fac: float = 1e9):
    # Set data directory
    data_dir = get_path("data", "calibration")

    # Read data file
    file_path = data_dir / f"calibration_{num_sites}_sites.csv"
    file_path = data_dir / f"calibration_{num_sites}_sites.csv"
    df = pd.read_csv(file_path)

    # Extract information
    z = df[f"z_{year}"].to_numpy()
    zbar = df["zbar_2017"].to_numpy()
    forest_area = df[f"area_forest_{year}"].to_numpy()

    # Normalize Z and forest data
    z /= norm_fac
    zbar /= norm_fac
    forest_area /= norm_fac

    return (zbar, z, forest_area)


def load_productivity_params(num_sites: int):
    data_dir = get_path("data", "calibration")
    productivity_params = pd.read_csv(data_dir / f"productivity_params_{num_sites}.csv")

    theta = productivity_params["theta_fit"]

    gamma = productivity_params["gamma_fit"]

    return (theta.to_numpy()[:,].flatten(), gamma.to_numpy()[:,].flatten())


def load_price_data():
    # Read data file
    file_path = get_path("data", "calibration") / "seriesPriceCattle_prepared.csv"
    df = pd.read_csv(file_path)
    average_prices = df.groupby("year")["price_real_mon_cattle"].mean()
    p_a_list = np.array(average_prices)
    return p_a_list


def load_site_data_1995(num_sites: int, norm_fac: float = 1e9):
    # Set data directory
    data_dir = get_path("data", "calibration")

    # Read data file
    file_path = data_dir / f"calibration_{num_sites}_sites.csv"
    df = pd.read_csv(file_path)

    # Extract information
    z_1995 = df["z_1995"].to_numpy()
    z_2008 = df["z_2008"].to_numpy()
    zbar_1995 = df["zbar_1995"].to_numpy()
    forest_area_1995 = df["area_forest_1995"].to_numpy()

    # Normalize Z data
    zbar_1995 /= norm_fac
    z_1995 /= norm_fac
    forest_area_1995 /= norm_fac

    (theta, gamma) = load_productivity_params(num_sites)

    print(f"Data successfully loaded from {data_dir}")
    return (
        zbar_1995,
        z_1995,
        forest_area_1995,
        z_2008,
        theta,
        gamma,
    )
