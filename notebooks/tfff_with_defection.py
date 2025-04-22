import numpy as np
from matplotlib import pyplot as plt

from pysrc.analysis import compute_planner_value
from pysrc.optimization import PlannerSolution, solve_planner_problem
from pysrc.services.data_service import load_productivity_params, load_site_data

b = 25
delta = 0.02
kappa = 2.094215255

# Short-term horizon
H = 50

solver = "gurobi"
pee = 6.6
pa = 41.11
num_sites = 1043
T = 200

# Load site data
(zbar, z_2017, forest_area_2017) = load_site_data(num_sites)

# Load productivity params
(theta, gamma) = load_productivity_params(num_sites)

# Computing carbon absorbed in start period
x_2017 = gamma * forest_area_2017

# Solve planner problem with b=$25
res = solve_planner_problem(
    x0=x_2017,
    z0=z_2017,
    zbar=zbar,
    gamma=gamma,
    theta=theta,
    time_horizon=T + H,
    price_cattle=pa,
    price_emissions=pee + b,
)

# Compute rolling value of continuing reforestaion
V = []
for t in range(H):
    res_t = PlannerSolution(
        Z=res.Z[t:],
        X=res.X[t:],
        U=res.U[t:],
        V=res.V[t:],
    )
    V.append(
        compute_planner_value(
            pee,
            pa,
            b,
            r=0,
            theta=theta,
            zbar=zbar,
            solution=res_t,
        )["planner_value"]
    )

# Plotting the line plot
plt.plot(V)

# Adding labels and title
plt.xlabel("Time")
plt.ylabel("Value")
plt.title("V(t)")

# Displaying the plot
plt.show()


W = []
for t in range(H):
    res_d = solve_planner_problem(
        x0=res.X[t],
        z0=res.Z[t],
        zbar=zbar,
        gamma=gamma,
        theta=theta,
        time_horizon=T,
        price_cattle=pa,
        price_emissions=pee,
    )

    W.append(
        compute_planner_value(
            pee,
            pa,
            b=0,
            r=0,
            theta=theta,
            zbar=zbar,
            solution=res_d,
        )["planner_value"]
    )


# Plotting the line plots for V and W
plt.plot(V, label="Continuing")
plt.plot(W, label="Defecting")

# Adding legend
plt.legend()

# Adding labels and title
plt.xlabel("Time (years)")
plt.ylabel("$ billion")
# plt.title("Value under b=$25 scheme v.s defecting")

# Displaying the plot
plt.show()

# Year of defection
indices = [i for i, (v, w) in enumerate(zip(V, W)) if v < w]
D = indices[0] if indices else -1

# Compute change in X
X_dot = np.diff(res.X, axis=0)

# Compute net transfers in year 40 dollars
net_transfers = [
    -b * (kappa * np.sum(res.Z[t + 1]) - np.sum(X_dot[t])) * ((1 + delta) ** (D - t))
    for t in range(D)
]
net_transfers = np.sum(net_transfers)

# Compute NPV of TFFF payments over 200 years
discount_factor = 1 / (1 + delta)
TFFF_fund = zbar.sum() * sum(5 * discount_factor**t for t in range(T))

F = zbar.sum() - res.Z.sum(axis=1)
mu = [(W[t] - V[t]) / F[t] if (W[t] - V[t]) > 0 else 0 for t in range(H)]

# Plotting the line plot
plt.plot(mu)

# Adding labels and title
plt.xlabel("Time")
plt.ylabel("$ / ha")
plt.title("mu(t)")

# Displaying the plot
plt.show()
