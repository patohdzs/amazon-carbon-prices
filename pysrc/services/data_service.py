import geopandas as gpd
import numpy as np
import pandas as pd

from ..services.file_service import get_path


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
    productivity_parameters = pd.read_csv(data_dir / f"productivity_params_{num_sites}.csv")
    
    theta = productivity_parameters["theta_fit"]

    gamma = productivity_parameters["gamma_fit"]

    return (theta.to_numpy()[:,].flatten(), gamma.to_numpy()[:,].flatten())


def load_reg_data(num_sites: int):
    # Set data directory
    data_dir = get_path("data", "calibration", "hmc")

    # Read site level data
    site_theta_df = gpd.read_file(data_dir / f"theta_fit_{num_sites}.geojson")
    site_gamma_df = gpd.read_file(data_dir / f"gamma_data_site_{num_sites}.geojson")

    # Remove geometries
    site_theta_df = site_theta_df.iloc[:, :-1]
    site_gamma_df = site_gamma_df.iloc[:, :-1]

    print(f"Data successfully loaded from {data_dir}")
    return (
        site_theta_df,
        site_gamma_df,
    )


def load_price_data():
    # Read data file
    file_path = (
        get_path("data", "calibration", "hmc") / "seriesPriceCattle_prepared.csv"
    )
    df = pd.read_csv(file_path)
    average_prices = df.groupby("year")["price_real_mon_cattle"].mean()
    p_a_list = np.array(average_prices)
    return p_a_list


def load_site_data_1995(num_sites: int, norm_fac: float = 1e9):
    # Set data directory
    data_dir = get_path("data", "calibration", "hmc")

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

    theta = (
        pd.read_csv(
            get_path("data", "calibration", "hmc") / f"theta_fit_{num_sites}.csv"
        )
        .to_numpy()[:,]
        .flatten()
    )
    gamma = (
        pd.read_csv(
            get_path("data", "calibration", "hmc") / f"gamma_fit_{num_sites}.csv"
        )
        .to_numpy()[:,]
        .flatten()
    )

    print(f"Data successfully loaded from {data_dir}")
    return (
        zbar_1995,
        z_1995,
        forest_area_1995,
        z_2008,
        theta,
        gamma,
    )
