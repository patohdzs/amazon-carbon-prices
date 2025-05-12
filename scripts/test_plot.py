import matplotlib.pyplot as plt
import pandas as pd

from pysrc.services.file_service import get_path

num_sites = 1043

data_dir = get_path("data", "calibration")
df = pd.read_csv(data_dir / f"distribution_parameters_all_{num_sites}.csv")


df["precision_eta_gamma"] = 1.0 / df["sigma_u_gamma"] ** 2
df["precision_eta_theta"] = 1.0 / df["sigma_u_theta"] ** 2
df["precision_zeta_gamma"] = 1.0 / df["sigma_v_gamma"] ** 2
df["precision_zeta_theta"] = 1.0 / df["sigma_v_theta"] ** 2


plot_dir = get_path("test", "plots") / f"histograms_{num_sites}"
plot_dir.mkdir(parents=True, exist_ok=True)


# Suppose this is your column list:
columns = df.columns.tolist()

latex_labels = []

for col in columns:
    # beta_gamma
    if col.startswith("beta_gamma_"):
        idx = col.split("_")[-1]
        latex_labels.append(r"$\beta_{\gamma," + idx + r"}$")
    # beta_theta
    elif col.startswith("beta_theta_"):
        idx = col.split("_")[-1]
        latex_labels.append(r"$\beta_{\theta," + idx + r"}$")
    # nu_gamma
    elif col.startswith("nu_gamma_"):
        idx = col.split("_")[-1]
        latex_labels.append(r"$\nu_{\gamma," + idx + r"}$")
    # nu_theta
    elif col.startswith("nu_theta_"):
        idx = col.split("_")[-1]
        latex_labels.append(r"$\nu_{\theta," + idx + r"}$")
    # sigma variables
    elif col == "sigma_u_gamma":
        latex_labels.append(r"$\sigma_{u,\gamma}$")
    elif col == "sigma_u_theta":
        latex_labels.append(r"$\sigma_{u,\theta}$")
    elif col == "sigma_v_gamma":
        latex_labels.append(r"$\sigma_{v,\gamma}$")
    elif col == "sigma_v_theta":
        latex_labels.append(r"$\sigma_{v,\theta}$")
    elif col == "precision_eta_gamma":
        latex_labels.append(r"$\eta_{\gamma}$")
    elif col == "precision_eta_theta":
        latex_labels.append(r"$\eta_{\theta}$")
    elif col == "precision_zeta_gamma":
        latex_labels.append(r"$\zeta_{\gamma}$")
    elif col == "precision_zeta_theta":
        latex_labels.append(r"$\zeta_{\theta}$")
    else:
        # fallback if any unknown column
        latex_labels.append(col)


# Loop through each column and plot
for i, col in enumerate(df.columns):
    plt.figure(figsize=(6, 4))
    plt.hist(df[col], bins=50)
    # plt.title(f"Histogram of {col}")
    plt.xlabel(latex_labels[i])
    plt.ylabel("Frequency")
    # plt.tight_layout(rect=[0, 0.1, 1, 1])
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)
    # Save each histogram
    plt.savefig(plot_dir / f"{col}_histogram.png", dpi=300)
    plt.close()  # close figure after saving
