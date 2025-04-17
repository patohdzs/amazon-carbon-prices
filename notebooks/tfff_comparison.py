from functools import partial

import matplotlib.pyplot as plt
import numpy as np

from pysrc.analysis import compute_planner_value
from pysrc.optimization import solve_planner_problem
from pysrc.services.data_service import load_productivity_params, load_site_data

solver = "gurobi"
pee = 6.6
pa = 41.11
num_sites = 1043
T = 200

(zbar, _, _) = load_site_data(num_sites)

# Load baseline productivity params
(theta, gamma) = load_productivity_params(num_sites)

# Set fully forested Amazon
z0 = np.zeros_like(zbar)

# Compute initial carbon stock
x0 = gamma * zbar

# Solve planner problem starting with fully forested Amazon
solve_tfff_problem = partial(
    solve_planner_problem,
    tfff_rent=0,
    x0=x0,
    z0=z0,
    zbar=zbar,
    gamma=gamma,
    theta=theta,
    time_horizon=T,
    price_cattle=pa,
    price_emissions=pee,
)

# Solve baseline problem
res_baseline = solve_tfff_problem()

# Compute aggregate land use as a percentage of site area
pct_Z_0 = 100 * (res_baseline.Z.sum(axis=1) / zbar.sum())

# Plotting the agricultural land use trajectories
fig = plt.figure()
plt.plot(pct_Z_0[:50], color="red", label=r"b=\$0, r=\$0")

# Adding legend
plt.legend()

# Adding labels and title
plt.xlabel("Time (years)")
plt.ylabel(r"$Z_t$ (%)")


# Solve TFFF problem
tfff_rent = 5
res = solve_tfff_problem(tfff_rent=tfff_rent)

# Compute aggregate land use as a percentage of site area
pct_Z = 100 * (res.Z.sum(axis=1) / zbar.sum())

# Plotting the agricultural land use trajectories
fig = plt.figure()
plt.plot(pct_Z[:50], color="red", label=rf"b=\$0, r=\${tfff_rent}")

# Adding legend
plt.legend()

# Adding labels and title
plt.xlabel("Time (years)")
plt.ylabel(r"$Z_t$ (%)")

# Set y-axis limits
plt.ylim(0, 5)

value_baseline = compute_planner_value(
    pee=pee,
    pa=pa,
    b=0,
    r=0,
    theta=theta,
    zbar=zbar,
    solution=res_baseline,
)

value_tfff = compute_planner_value(
    pee=pee,
    pa=pa,
    b=0,
    r=tfff_rent,
    theta=theta,
    zbar=zbar,
    solution=res,
)
print(f"value baseline: {value_baseline}")
print(f"value tfff: {value_tfff}")


f = tfff_rent * np.sum(zbar)
delta = 0.02

# Sum of geometric series: f * (1 - r^n) / (1 - r)
total_payoff = sum(f / ((1 + delta) ** t) for t in range(T))

print(f"Total discounted payoff over {T} periods: {total_payoff}")
