from pysrc.analysis.figures import (
    plot_agg_carbon_stock,
    plot_agg_land_use,
    plot_transfers,
)
from pysrc.analysis.tables import make_pvd_table, make_transfers_table
from pysrc.optimization import solve_planner_problem
from pysrc.services.data_service import load_productivity_params, load_site_data
from pysrc.services.file_service import get_path

# Model hyperparameters
solver = "gurobi"
pee = 6.1
pa = 41.11
num_sites = 78
T = 200

# Load site data
(zbar_2017, z_2017, forest_area_2017) = load_site_data(num_sites)

# Load productivity params
(theta, gamma) = load_productivity_params(num_sites)

# Computing carbon absorbed in start period
x_2017 = gamma * forest_area_2017

# Solve planner problem with different transfers
results = []
for b in range(0, 30, 5):
    results.append(
        solve_planner_problem(
            x0=x_2017,
            z0=z_2017,
            zbar=zbar_2017,
            gamma=gamma,
            theta=theta,
            time_horizon=T,
            price_cattle=pa,
            price_emissions=pee + b,
        )
    )

# Set output directory
output_dir = get_path("results")
suffix = f"sites{num_sites}_pa{pa}_pee{pee}"

# Plot aggregate carboon stock trajectories
fig = plot_agg_carbon_stock(results[0], results[3], results[5])
fig.savefig(output_dir / f"agg_X_det_{suffix}.pdf", bbox_inches="tight")

# Plot aggregate land use trajectories
fig = plot_agg_land_use(results[0], results[3], results[5], zbar_2017)
fig.savefig(output_dir / f"agg_Z_det_{suffix}.pdf", bbox_inches="tight")

# Plot net transfers
fig = plot_transfers(results[3], results[5])
fig.savefig(output_dir / f"net_transfers_det_{suffix}.pdf", bbox_inches="tight")

# Create PVD LaTex table
pvd_table = make_pvd_table(results, pee, pa, theta)
with open(output_dir / f"pvd_table_{suffix}.tex", "w") as f:
    f.write(pvd_table)

# Create transfers LaTeX table
transfers_table = make_transfers_table(results, pee)
with open(output_dir / f"transfers_table_{suffix}.tex", "w") as f:
    f.write(transfers_table)
