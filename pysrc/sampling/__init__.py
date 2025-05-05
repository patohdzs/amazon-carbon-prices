import numpy as np
import pandas as pd

from ..services.file_service import get_path


# def theta_adj_reg_data(num_sites, theta_df):
#     M = 62

#     # Get design matrix and its dimensions
#     X = theta_df.iloc[:, :8].to_numpy()
#     N, K = X.shape

#     # Large group indicator
#     G = theta_df.iloc[:, 8]
#     G = G.astype(int)

#     # Municipal level to site level projection matrix
#     SG = np.array(
#         [(theta_df["id"].to_numpy() == i).astype(int) for i in range(1, num_sites + 1)]
#     )

#     # Multiply by area overalp weights
#     SG = theta_df["muni_site_area"].to_numpy() * SG
#     SG = SG / SG.sum(axis=1, keepdims=True)

#     return {
#         "M_theta": M,
#         "X_theta": X,
#         "N_theta": N,
#         "K_theta": K,
#         "SG_theta": SG,
#         "G_theta": G,
#     }


# def gamma_adj_reg_data(num_sites, gamma_df):
#     # Get design matrix and its dimensions
#     M = 78
#     X = gamma_df.iloc[:, :6].to_numpy()
#     N, K = X.shape

#     # Large group indicator
#     G = gamma_df.iloc[:, 7]

#     return {
#         "M_gamma": M,
#         "X_gamma": X,
#         "N_gamma": N,
#         "K_gamma": K,
#         "G_gamma": G,
#     }


# def load_gamma_calib(num_sites: int, type: str = "reg"):
#     data_dir = get_path("data", "calibration")
#     if type == "fit":
#         # Used for gamma projection onto fitted values
#         df = gpd.read_file(data_dir / f"gamma_fit_{num_sites}.geojson")

#         # Get design matrix and its dimensions
#         X = df.iloc[:, :6].to_numpy()
#         N, K = X.shape

#         # Large group indicator
#         m = df["id_group"].astype(int)

#         return {
#             "X_gamma_fit": X,
#             "m_gamma_fit": m,
#         }
#     else:
#         # Used for gamma regression
#         df = gpd.read_file(data_dir / "gamma_reg.geojson")

#         # Get design matrix and its dimensions
#         M = df["id_group"].unique().size
#         y = df["log_co2e_ha_2017"]
#         X = df.iloc[:, :6].to_numpy()
#         N, K = X.shape

#         # Large group indicator
#         m = df["id_group"].astype(int)

#         return {
#             "N_gamma": N,
#             "M_gamma": M,
#             "K_gamma": K,
#             "y_gamma": y,
#             "X_gamma": X,
#             "m_gamma": m,
#         }


# def load_theta_calib(num_sites: int, type: str = "reg"):
#     data_dir = get_path("data", "calibration")
#     if type == "fit":
#         df = gpd.read_file(data_dir / f"theta_fit_{num_sites}.geojson")

#         # Create the projection matrix
#         G = (
#             df.pivot(index="id", columns="muni_id", values="muni_site_area")
#             .fillna(0)
#             .to_numpy()
#         )

#         # Normalize to make row-stochastic
#         G = G / G.sum(axis=1, keepdims=True)

#         # Collapse data set to the municipality level
#         df = df.sort_values("muni_id")

#         # Keep first observation per municipality
#         df = df.drop_duplicates(subset="muni_id", keep="first")

#         # Get design matrix
#         X = df.iloc[:, :8].to_numpy()
#         C, _ = X.shape

#         # Large group indicator
#         m = df["group_id"].astype(int)

#         # Cattle price in 2017
#         pa_2017 = 44.9736197781184
        
        
#         G_sparse = coo_matrix(G)
#         row_G_theta = G_sparse.row + 1  # Stan uses 1-based indexing
#         col_G_theta = G_sparse.col + 1
#         val_G_theta = G_sparse.data
#         N_nonzero_G_theta = len(val_G_theta)
        
#         return {
#             "C_theta_fit": C,
#             "X_theta_fit": X,
#             "m_theta_fit": m,
#             # "G_theta_fit": G,
#             "N_nonzero_G_theta": N_nonzero_G_theta,
#             "row_G_theta": row_G_theta,
#             "col_G_theta": col_G_theta,
#             "val_G_theta": val_G_theta,
#             "pa_2017": pa_2017,
#         }

#     else:
#         df = gpd.read_file(data_dir / "theta_reg.geojson")
#         # Get number of groups
#         M = df["group_id"].unique().size

#         # Get design matrix and its dimensions
#         y = df["log_slaughter"]
#         X = df.iloc[:, :8].to_numpy()
#         N, K = X.shape

#         # Large group indicator
#         m = df["group_id"].astype(int)

#         W = df["weights"].values
#         W= np.sqrt(W/np.std(W))

#         return {
#             "N_theta": N,
#             "M_theta": M,
#             "K_theta": K,
#             "y_theta": y,
#             "X_theta": X,
#             "m_theta": m,
#             "W_theta": W,
#         }



# # def baseline_hyperparams(var):
#     # Load baseline hyperparameters
#     data_dir = get_path("data", "calibration", "hmc")
#     df = pd.read_csv(data_dir / f"{var}_posterior_means_and_variances.csv")

#     beta_mean = df[df["Parameter"].str.startswith("beta_")]["Posterior Mean"].values
#     V_mean = df[df["Parameter"].str.startswith("V_")]["Posterior Mean"].values

#     # Construct coefficients vcov
#     beta_cov = np.zeros((len(beta_mean), len(beta_mean)))
#     for i in range(len(beta_mean)):
#         for j in range(len(beta_mean)):
#             param_name = f"var_beta_{i+1}_{j+1}"
#             beta_cov[i, j] = df[df["Parameter"] == param_name][
#                 "Posterior Variance"
#             ].values[0]

#     V_var = df[df["Parameter"].str.startswith("V_")]["Posterior Variance"].values

#     return {
#         f"beta_{var}_mean": beta_mean,
#         f"beta_{var}_vcov": beta_cov,
#         f"V_{var}_mean": V_mean,
#         f"V_{var}_var": V_var,
#     }
