import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from pysrc.services.data_service import load_site_data
from pysrc.services.file_service import get_path

# Allow TeX for plotting
plt.rcParams["text.usetex"] = True
plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.size": 10,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
    }
)


def land_allocation(pee=7.6, num_sites=1043, opt="gurobi", pa=41.11, model="det", xi=1):
    # Set transfer levels
    b = [0, 10, 15, 20, 25]

    # Set corresponding emissions prices
    pe = [pee + bi for bi in b]

    # Get z_bar data
    output_folder = get_path("output") / "figures"

    # Load site data
    (zbar, _, _) = load_site_data(num_sites)

    variable_dict = {}
    for j in range(5):
        order = j
        results_dir = (
            get_path("output")
            / "optimization"
            / model
            / opt
            / f"{num_sites}sites"
            / f"pa_{pa}"
            / f"pe_{pe[order]}"
        )

        dfz = np.loadtxt(results_dir / "Z.txt", delimiter=",")
        dfx = np.sum(np.loadtxt(results_dir / "X.txt", delimiter=","), axis=1)

        if model == "det":
            result_folder = os.path.join(
                str(get_path("output")),
                "optimization",
                model,
                opt,
                f"{num_sites}sites",
                f"pa_{pa}",
                f"pe_{pe[order]}",
            )
        elif model == "hmc":
            result_folder = os.path.join(
                str(get_path("output")),
                "optimization",
                model,
                opt,
                f"{num_sites}sites",
                f"xi_{xi}",
                f"pa_{pa}",
                f"pe_{pe[order]}",
            )
        dfz = np.loadtxt(os.path.join(result_folder, "Z.txt"), delimiter=",")
        dfx = np.sum(
            np.loadtxt(os.path.join(result_folder, "X.txt"), delimiter=","), axis=1
        )

        variable_dict[f"results_zper{j}"] = []
        variable_dict[f"results_xagg{j}"] = dfx[:51]
        for i in range(51):
            result_zper = (np.sum(dfz[i]) / np.sum(zbar)) * 100
            variable_dict[f"results_zper{j}"].append(result_zper)

    time = list(range(len(variable_dict[f"results_zper{0}"])))
    plt.figure(figsize=(10, 6))
    custom_labels = [
        "$p^{{ee}}$={0}       $b$".format(pee),
        "0",
        "10",
        "15",
        "20",
        "25",
    ]

    plt.plot([], [], " ", label=custom_labels[0])
    for i in range(5):
        if i in [0, 2, 4]:
            if i == 0:
                color = "red"
            elif i == 2:
                color = "green"
            elif i == 4:
                color = "blue"
            plt.plot(
                time,
                variable_dict[f"results_zper{i}"],
                label=custom_labels[i + 1],
                linewidth=4,
                color=color,
            )

    plt.xlabel("years", fontsize=16)
    plt.ylabel("Z(%)", fontsize=16)
    plt.xlim(0, max(time) + 2)
    plt.yticks([0, 5, 10, 15, 20, 25], ["0", "5", "10", "15", "20", "25"])
    plt.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=5,
        frameon=False,
        fontsize=16,
    )
    plt.savefig(
        output_folder / f"pred_zshare_{num_sites}_sites_det.png",
        format="png",
        bbox_inches="tight",
    )
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.plot([], [], " ", label=custom_labels[0])
    for i in range(5):
        if i in [0, 2, 4]:
            if i == 0:
                color = "red"
            elif i == 2:
                color = "green"
            elif i == 4:
                color = "blue"
            plt.plot(
                time,
                variable_dict[f"results_xagg{i}"],
                label=custom_labels[i + 1],
                linewidth=4,
                color=color,
            )
    plt.xlabel("years", fontsize=18)
    plt.ylabel("X(billions CO2e)", fontsize=18)
    plt.xlim(0, max(time) + 2)
    plt.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=5,
        frameon=False,
        fontsize=18,
    )
    plt.savefig(
        output_folder + f"/plot_pred_x_{num_sites}_sites_det.png",
        format="png",
        bbox_inches="tight",
    )
    plt.show()


def density(pee=7.6, num_sites=78, solver="gurobi", pa=41.11, xi=1, model="det"):
    output_folder = (
        str(get_path("output")) + f"/figures/density/site_{num_sites}/xi{xi}/"
    )
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    result_folder = os.path.join(
        str(get_path("output")),
        "sampling",
        solver,
        f"{num_sites}sites",
        f"pa_{pa}",
        f"xi_{xi}",
    )
    prior_folder = os.path.join(
        str(get_path("output")),
        "sampling",
        solver,
        f"{num_sites}sites",
        f"pa_{pa}",
        "xi_10000",
    )

    with open(result_folder + f"/pe_{pee}/results.pcl", "rb") as f:
        b0 = pickle.load(f)

    with open(result_folder + f"/pe_{pee+15}/results.pcl", "rb") as f:
        b15 = pickle.load(f)

    with open(prior_folder + f"/pe_{pee+15}/results.pcl", "rb") as f:
        results_unadjusted = pickle.load(f)

    theta_unadjusted = results_unadjusted["final_sample"][:16000, :num_sites]
    gamma_unadjusted = results_unadjusted["final_sample"][:16000, num_sites:]
    theta_adjusted_b0 = b0["final_sample"][:16000, :num_sites]
    gamma_adjusted_b0 = b0["final_sample"][:16000, num_sites:]
    theta_adjusted_b15 = b15["final_sample"][:16000, :num_sites]
    gamma_adjusted_b15 = b15["final_sample"][:16000, num_sites:]

    for idx in range(num_sites):
        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        num_bins = 100

        global_min = min(
            gamma_unadjusted[:, idx].min(),
            gamma_adjusted_b0[:, idx].min(),
            gamma_adjusted_b15[:, idx].min(),
        )
        global_max = max(
            gamma_unadjusted[:, idx].max(),
            gamma_adjusted_b0[:, idx].max(),
            gamma_adjusted_b15[:, idx].max(),
        )
        np.linspace(global_min, global_max, num_bins + 1)

        sns.kdeplot(
            gamma_unadjusted[:, idx],
            label="baseline",
            color="black",
            fill=False,
            alpha=0.6,
            linewidth=4,
        )  # Bright blue
        sns.kdeplot(
            gamma_adjusted_b0[:, idx],
            label="b=0",
            color="red",
            fill=True,
            alpha=0.6,
            linewidth=4,
        )  # Bright red
        sns.kdeplot(
            gamma_adjusted_b15[:, idx],
            label="b=15",
            color="blue",
            fill=True,
            alpha=0.6,
            linewidth=4,
        )  # Bright red

        plt.title(rf"Probability density for $\gamma$ and site {idx+1}", fontsize=16)
        plt.xlabel("parameter value", fontsize=16)
        plt.ylabel("density", fontsize=16)
        plt.legend(fontsize=16)
        file_name = os.path.join(output_folder, f"gamma_distribution_{idx+1}.png")
        fig.savefig(file_name, format="png")
        plt.close()

    for idx in range(num_sites):
        fig, axes = plt.subplots(1, 1, figsize=(8, 6))
        num_bins = 100

        global_min = min(
            theta_unadjusted[:, idx].min(),
            theta_adjusted_b0[:, idx].min(),
            theta_adjusted_b15[:, idx].min(),
        )
        global_max = max(
            theta_unadjusted[:, idx].max(),
            theta_adjusted_b0[:, idx].max(),
            theta_adjusted_b15[:, idx].max(),
        )
        np.linspace(global_min, global_max, num_bins + 1)

        sns.kdeplot(
            theta_unadjusted[:, idx],
            label="baseline",
            color="black",
            fill=False,
            alpha=0.6,
            linewidth=4,
        )  # Bright blue
        sns.kdeplot(
            theta_adjusted_b0[:, idx],
            label="b=0",
            color="red",
            fill=True,
            alpha=0.6,
            linewidth=4,
        )  # Bright red
        sns.kdeplot(
            theta_adjusted_b15[:, idx],
            label="b=15",
            color="blue",
            fill=True,
            alpha=0.6,
            linewidth=4,
        )  # Bright red

        plt.title(rf"Probability density for $\Theta$ and site {idx+1}", fontsize=16)
        plt.xlabel("parameter value", fontsize=16)
        plt.ylabel("density", fontsize=16)
        plt.legend(fontsize=16)
        file_name = os.path.join(output_folder, f"theta_distribution_{idx+1}.png")
        fig.savefig(file_name, format="png")
        plt.close()
    return


def trajectory_diff(
    num_sites=78, pe_hmc=7.1, pe_det=5.3, b=0, solver="gams", pa=41.11, xi=1
):
    pe_hmc += b
    pe_det += b

    result_folder = os.path.join(str(get_path("output")), "optimization")
    output_folder = str(get_path("output")) + "/figures"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    df_ori = pd.read_csv(
        str(get_path("data")) + f"/calibration/hmc/calibration_{num_sites}_sites.csv"
    )
    dfz_bar = df_ori["zbar_2017"]
    dfz_bar_np = dfz_bar.to_numpy()

    dfz_hmc = np.loadtxt(
        os.path.join(
            result_folder
            + "/hmc/"
            + solver
            + f"/{num_sites}sites/xi_{xi}/pa_{pa}/pe_{pe_hmc}/",
            "Z.txt",
        ),
        delimiter=",",
    )

    dfz_zeronp_hmc = dfz_hmc
    result_zper_hmc = np.zeros((51, 1))
    result_zper_hmc = result_zper_hmc[:, 0]
    for i in range(51):
        result_zper_hmc[i] = (
            np.sum(dfz_zeronp_hmc[i]) / (np.sum(dfz_bar_np) / 1e9)
        ) * 100

    dfz_det = np.loadtxt(
        os.path.join(
            result_folder
            + "/det/"
            + solver
            + f"/{num_sites}sites/pa_{pa}/pe_{pe_det}/",
            "Z.txt",
        ),
        delimiter=",",
    )

    dfz_zeronp_det = dfz_det
    result_zper_det = np.zeros((51, 1))
    result_zper_det = result_zper_det[:, 0]
    for i in range(51):
        result_zper_det[i] = (
            np.sum(dfz_zeronp_det[i]) / (np.sum(dfz_bar_np) / 1e9)
        ) * 100

    time = list(range(0, len(result_zper_hmc)))
    plt.figure(figsize=(10, 6))

    plt.plot(time, result_zper_hmc, label=rf"$\xi$={xi}", linewidth=4, color="blue")
    plt.plot(time, result_zper_det, label=r"$\xi=\infty$", linewidth=4, color="red")
    plt.xlabel("years", fontsize=16)
    plt.ylabel("Z(%)", fontsize=16)
    plt.xlim(0, max(time) + 2)
    plt.legend(loc="upper left", ncol=5, frameon=False, fontsize=16)
    plt.savefig(
        output_folder
        + f"/aggregate_percentage_Z_b{b}_pehmc_{pe_hmc}_pedet_{pe_det}.png"
    )
    plt.show()

    return


def plot_agg_land_use(results_0, results_15, results_25, zbar):
    # Compute aggregate land use as a percentage of site area
    pct_Z_0 = 100 * (results_0.Z.sum(axis=1) / zbar.sum())
    pct_Z_15 = 100 * (results_15.Z.sum(axis=1) / zbar.sum())
    pct_Z_25 = 100 * (results_25.Z.sum(axis=1) / zbar.sum())

    # Plotting the agricultural land use trajectories
    fig = plt.figure(figsize=_set_size())
    plt.plot(pct_Z_0[:50], color="red", label=r"b=\$0")
    plt.plot(pct_Z_15[:50], color="green", label=r"b=\$15")
    plt.plot(pct_Z_25[:50], color="blue", label=r"b=\$25")

    # Adding legend
    plt.legend()

    # Adding labels and title
    plt.xlabel("Time (years)")
    plt.ylabel(r"$Z_t$ (\%)")
    return fig


def plot_agg_carbon_stock(results_0, results_15, results_25):
    # Plotting the carbon stock trajectories
    fig = plt.figure(figsize=_set_size())
    plt.plot(results_0.X.sum(axis=1)[:50], color="red", label=r"b=\$0")
    plt.plot(results_15.X.sum(axis=1)[:50], color="green", label=r"b=\$15")
    plt.plot(results_25.X.sum(axis=1)[:50], color="blue", label=r"b=\$25")

    # Adding legend
    plt.legend()

    # Adding labels and title
    plt.xlabel("Time (years)")
    plt.ylabel(r"$X_t$ (CO$_2$e Gt)")
    return fig


def plot_transfers(results_15, results_25, kappa=2.094215255):
    fig = plt.figure(figsize=_set_size())

    for b, res, color in zip([15, 25], [results_15, results_25], ["red", "blue"]):
        kappa = 2.094215255

        # Get X and Z trajectories
        X = res.X
        Z = res.Z

        # Compute X_dot
        X_dot = np.diff(X, axis=0)

        # Compute net transfers
        transfers = -b * (kappa * Z[1:] - X_dot).sum(axis=1)

        # Plot net transfers
        plt.plot(transfers[:50], color=color, label=rf"b=\${b}")

    # Adding legend
    plt.legend()

    # Adding labels and title
    plt.xlabel("Time (years)")
    plt.ylabel(r"Net Transfers (\$ billion)")

    return fig


def _set_size(width=397.485, scale=0.5):
    """
    Parameters
    ----------
    width: float
            Document textwidth or columnwidth in pts
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy

    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width * scale

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    golden_ratio = (5**0.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt

    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio

    fig_dim = (fig_width_in, fig_height_in)

    return fig_dim
