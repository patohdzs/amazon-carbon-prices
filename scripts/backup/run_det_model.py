import matplotlib.pyplot as plt
from pysrc.services.data_service import load_productivity_params
from pysrc.optimization import solve_planner_problem
from pysrc.sampling import baseline
from pysrc.services.data_service import load_site_data

# Model hyperparameters
solver = "gurobi"
pee = 7.6
pa = 41.11
num_sites = 78
T = 200

# Load site data
(zbar_2017, z_2017, forest_area_2017) = load_site_data(num_sites)

# Set productivity parameters using baseline mean
baseline_fit = baseline.sample(
    num_sites=num_sites,
    iter_sampling=10**3,
    chains=5,
    seed=1,
)

(theta, gamma) = load_productivity_params(num_sites)
# Computing carbon absorbed in start period
x_2017 = gamma * forest_area_2017

# Solve planner problem
results = []
for b in range(0, 5, 5):
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

# Plotting the line plots for V and W
for i, res in enumerate(results):
    plt.plot(res.X.sum(axis=1)[:50], label=f"b=${i*5}")

# Adding legend
plt.legend()

# Adding labels and title
plt.xlabel("Time (years)")
plt.ylabel("Carbon stock (Gigatons)")
