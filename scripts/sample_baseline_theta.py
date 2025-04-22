import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.stats import gamma
from tqdm import tqdm

from pysrc.sampling import theta_adj_reg_data
from pysrc.services.file_service import get_path

data_dir = get_path("data", "calibration", "hmc")

theta_reg = gpd.read_file(data_dir / "theta_reg.geojson")

X = theta_reg.iloc[:, 0:8].values  # Columns 1 to 6 as X
Y = theta_reg.iloc[:, 8].values  # Column 7 as Y
site_ids = theta_reg.iloc[:, -2].values - 1  # Column 8 as site_ids
site_ids = site_ids.astype(int)
weights = theta_reg.iloc[:, 9].values
weights = weights / np.std(weights)


np.random.seed(45)


def gibbs_sampling_with_data(
    X,
    Y,
    site_ids,
    weights,
    n_iterations=500000,
    beta_prior_mean=None,
    beta_prior_cov=None,
    eta_prior_shape=2,
    eta_prior_scale=1,
    nu_prior_shape=2,
    nu_prior_scale=1,
):
    n_data, n_features = X.shape
    n_sites = len(np.unique(site_ids))

    # Initialize priors
    if beta_prior_mean is None:
        beta_prior_mean = np.zeros(n_features)
    if beta_prior_cov is None:
        beta_prior_cov = np.eye(n_features)  # Large prior variance

    # Initialize Gibbs sampling arrays
    beta_samples = np.zeros((n_iterations, n_features))
    eta_samples = np.zeros(n_iterations)
    nu_samples = np.zeros(n_iterations)
    V_samples = np.zeros((n_iterations, n_sites))

    # Initial values based on prior
    beta_current = np.random.multivariate_normal(beta_prior_mean, beta_prior_cov)
    eta_current = gamma.rvs(a=eta_prior_shape, scale=1 / eta_prior_scale)
    nu_current = gamma.rvs(a=nu_prior_shape, scale=1 / nu_prior_scale)
    V_current = np.zeros(n_sites) + 0.01

    sq_weight_matrix = np.diag(np.sqrt(weights))
    weight_X = sq_weight_matrix @ X
    weight_Y = sq_weight_matrix @ Y

    # Gibbs sampler
    for i in tqdm(range(n_iterations), desc="Gibbs Sampling Progress"):
        # Step i: Draw beta conditioned on V, eta, and data
        precision_beta = eta_current * (weight_X.T @ weight_X)
        mean_beta = np.linalg.inv(precision_beta) @ (
            eta_current
            * weight_X.T
            @ (weight_Y - sq_weight_matrix @ V_current[site_ids])
        )
        beta_current = np.random.multivariate_normal(
            mean_beta, np.linalg.inv(precision_beta)
        )
        beta_samples[i, :] = beta_current

        # Step ii: Draw eta conditioned on V and data
        residuals = (
            weight_Y - weight_X @ beta_current - sq_weight_matrix @ V_current[site_ids]
        )
        eta_shape_post = n_data / 2
        eta_scale_post = 0.5 * np.sum(residuals**2)
        eta_current = gamma.rvs(a=eta_shape_post, scale=1 / eta_scale_post)
        eta_samples[i] = eta_current

        # Step iii: Draw nu conditioned on V
        nu_shape_post = n_sites / 2
        nu_scale_post = 0.5 * np.sum(V_current**2)
        nu_current = gamma.rvs(a=nu_shape_post, scale=1 / nu_scale_post)
        # nu_current = 9.225
        nu_samples[i] = nu_current

        # nu_current = 1e10
        # nu_samples[i] = nu_current

        # Step iv: Draw V_j conditioned on beta, nu, eta, and data
        for j in range(n_sites):
            site_mask = site_ids == j
            # n_j = np.sum(site_mask)
            n_j = np.sum(weights[site_mask])
            weight_matrix_sub = np.diag(weights[site_mask])
            V_bar_j = np.sum(
                weight_matrix_sub @ Y[site_mask]
                - weight_matrix_sub @ X[site_mask] @ beta_current
            )
            mean_V_j = (eta_current * V_bar_j) / (n_j * eta_current + nu_current)
            precision_V_j = n_j * eta_current + nu_current
            V_current[j] = np.random.normal(mean_V_j, np.sqrt(1 / precision_V_j))
        V_samples[i, :] = V_current

    beta_samples_final = beta_samples[300000:, :]
    eta_samples_final = eta_samples[300000:]
    nu_samples_final = nu_samples[300000:]
    V_samples_final = V_samples[300000:, :]

    beta_posterior_mean = np.mean(beta_samples_final, axis=0)
    beta_posterior_variance = np.cov(beta_samples_final, rowvar=False, ddof=0)

    eta_posterior_mean = np.mean(eta_samples_final)
    eta_posterior_variance = np.var(eta_samples_final, ddof=0)

    nu_posterior_mean = np.mean(nu_samples_final)
    nu_posterior_variance = np.var(nu_samples_final, ddof=0)

    V_posterior_mean = np.mean(V_samples_final, axis=0)
    V_posterior_variance = np.var(V_samples_final, ddof=0, axis=0)

    percentiles = np.percentile(beta_samples_final, [10, 50, 90], axis=0)

    quan_table = pd.DataFrame({
        "10th_percentile": percentiles[0],
        "50th_percentile": percentiles[1],
        "90th_percentile": percentiles[2]
    })

    # Save to CSV
    quan_table.to_csv(get_path("output", "tables") /"theta_percentiles.csv", index=False)



    return (
        beta_posterior_mean,
        eta_posterior_mean,
        nu_posterior_mean,
        V_posterior_mean,
        beta_posterior_variance,
        eta_posterior_variance,
        nu_posterior_variance,
        V_posterior_variance,
    )


# Run Gibbs sampling with the extracted data
# Assuming X, Y, and site_ids are already defined
(
    beta_posterior_mean,
    eta_posterior_mean,
    nu_posterior_mean,
    V_posterior_mean,
    beta_posterior_variance,
    eta_posterior_variance,
    nu_posterior_variance,
    V_posterior_variance,
) = gibbs_sampling_with_data(X=X, Y=Y, site_ids=site_ids, weights=weights)

# Output the results
print(f"Posterior mean of beta: {beta_posterior_mean}")
print(f"Posterior mean of eta: {eta_posterior_mean}")
print(f"Posterior mean of nu: {nu_posterior_mean}")


data = []


for i in range(len(beta_posterior_mean)):
    data.append(
        {
            "Parameter": f"beta_{i+1}",
            "Posterior Mean": beta_posterior_mean[i],
            "Posterior Variance": None,
        }
    )


for i in range(len(beta_posterior_mean)):
    for j in range(len(beta_posterior_mean)):
        data.append(
            {
                "Parameter": f"var_beta_{i+1}_{j+1}",
                "Posterior Mean": None,
                "Posterior Variance": beta_posterior_variance[i, j],
            }
        )


data.append(
    {
        "Parameter": "eta",
        "Posterior Mean": eta_posterior_mean,
        "Posterior Variance": eta_posterior_variance,
    }
)


data.append(
    {
        "Parameter": "nu",
        "Posterior Mean": nu_posterior_mean,
        "Posterior Variance": nu_posterior_variance,
    }
)


for i in range(len(V_posterior_mean)):
    data.append(
        {
            "Parameter": f"V_{i+1}",
            "Posterior Mean": V_posterior_mean[i],
            "Posterior Variance": V_posterior_variance[i],
        }
    )


df = pd.DataFrame(data)


df.to_csv(data_dir / "theta_posterior_means_and_variances.csv", index=False)

## fit values

beta = np.random.multivariate_normal(
    beta_posterior_mean, beta_posterior_variance, 100000
)
V = np.random.normal(
    V_posterior_mean, np.sqrt(V_posterior_variance), (100000, len(np.unique(site_ids)))
)


mean_pa_2017 = 44.97362
theta_fit_df_78 = gpd.read_file(data_dir / "theta_fit_78.geojson")
X_fit = theta_adj_reg_data(78, theta_fit_df_78)["X_theta"]
weights = theta_adj_reg_data(78, theta_fit_df_78)["SG_theta"]
fit_ids = theta_fit_df_78["group_id"].values - 1  # Column 8 as site_ids
fit_ids = fit_ids.astype(int)
site_theta = np.mean(np.exp(X_fit @ beta.T + V[:, fit_ids].T), axis=1)
theta_78 = weights @ site_theta / 44.97362
theta_78_df = pd.DataFrame(theta_78, columns=["theta_fit"])
theta_78_df.to_csv(data_dir / "theta_fit_78.csv", index=False)

theta_fit_df_1043 = gpd.read_file(data_dir / "theta_fit_1043.geojson")
X_fit = theta_adj_reg_data(1043, theta_fit_df_1043)["X_theta"]
weights = theta_adj_reg_data(1043, theta_fit_df_1043)["SG_theta"]
fit_ids = theta_fit_df_1043["group_id"].values - 1  # Column 8 as site_ids
fit_ids = fit_ids.astype(int)
site_theta = np.mean(np.exp(X_fit @ beta.T + V[:, fit_ids].T), axis=1)
theta_1043 = weights @ site_theta / 44.97362
theta_1043_df = pd.DataFrame(theta_1043, columns=["theta_fit"])
theta_1043_df.to_csv(data_dir / "theta_fit_1043.csv", index=False)
