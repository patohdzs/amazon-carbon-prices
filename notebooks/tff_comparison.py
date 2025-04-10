from functools import partial

import matplotlib.pyplot as plt
import numpy as np

from pysrc.optimization import solve_planner_problem
from pysrc.services.data_service import load_productivity_params, load_site_data

solver = "gurobi"
pee = 6.6
pa = 41.11
num_sites = 1043
T = 200

(zbar, z_2017, forest_area_2017) = load_site_data(num_sites)

# Load baseline productivity params
(theta, gamma) = load_productivity_params(num_sites)

# Compute initial carbon stock
x_2017 = gamma * forest_area_2017

# Set fully forested Amazon
z0 = np.zeros_like(z_2017)

# Solve planner problem starting with fully forested Amazon
solve_tfff_problem = partial(
    solve_planner_problem,
    tfff_rent=0,
    x0=x_2017,
    z0=z0,
    zbar=zbar,
    gamma=gamma,
    theta=theta,
    time_horizon=T,
    price_cattle=pa,
    price_emissions=pee,
)


res_0 = solve_tfff_problem()

# Compute aggregate land use as a percentage of site area
pct_Z_0 = 100 * (res_0.Z.sum(axis=1) / zbar.sum())

# Plotting the agricultural land use trajectories
fig = plt.figure()
plt.plot(pct_Z_0[:50], color="red", label=r"b=\$0")

# Adding legend
plt.legend()

# Adding labels and title
plt.xlabel("Time (years)")
plt.ylabel(r"$Z_t$ (%)")
