import os
import pickle

from pysrc.analysis import value_decomposition
from pysrc.optimization import solve_planner_problem
from pysrc.sampling import adjusted
from pysrc.services.data_service import load_site_data
from pysrc.services.file_service import get_path

# ## Model scenario
opt = "gurobi"  # need to install gurobi solver
pee = 7.1
pa = 41.11
num_sites = 78
T = 200
b = [0, 10, 15, 20, 25]
pe_values = [pee + bi for bi in b]
xi = 1


# conduct hmc sampling, details see adjusted.py
for pe in pe_values:
    results = adjusted.sample(
        xi=xi,
        pe=pe,
        pa=pa,
        weight=0.25,
        num_sites=num_sites,
        T=200,
        optimizer=opt,
        max_iter=100,
        final_sample_size=5_000,
        iter_sampling=1000,
        iter_warmup=500,
        show_progress=True,
        seed=1,
        inits=0.2,
    )

    output_base_path = os.path.join(
        str(get_path("output")),
        "sampling",
        opt,
        f"{num_sites}sites",
        f"pa_{pa}",
        f"xi_{xi}",
        f"pe_{pe}",
    )
    if not os.path.exists(output_base_path):
        os.makedirs(output_base_path)
    outfile_path = output_base_path + "/results.pcl"
    with open(outfile_path, "wb") as outfile:
        pickle.dump(results, outfile)
        print(f"Results saved to {outfile_path}")

# Load site data
(zbar_2017, z_2017, forest_area_2017) = load_site_data(num_sites)

# Read ambiguity-adjusted parameters
result_folder = os.path.join(
    str(get_path("output")),
    "sampling",
    opt,
    f"{num_sites}sites",
    f"pa_{pa}",
    f"xi_{xi}",
)


# Solving the model with sampled parameters
b = [0, 10, 15, 20, 25]
pe_values = [pee + bi for bi in b]
for pe in pe_values:
    with open(result_folder + f"/pe_{pe}/results.pcl", "rb") as f:
        para = pickle.load(f)

    theta_vals = para["final_sample"][:16000, :78].mean(axis=0)
    gamma_vals = para["final_sample"][:16000, 78:].mean(axis=0)
    x0_vals = gamma_vals * forest_area_2017

    b_val = pe - pee
    results = solve_planner_problem(
        time_horizon=T,
        theta=theta_vals,
        gamma=gamma_vals,
        x0=x0_vals,
        zbar=zbar_2017,
        z0=z_2017,
        price_emissions=pe,
        price_cattle=pa,
    )

    print(
        "result",
        value_decomposition(
            Z=results.Z,
            X=results.X,
            U=results.U,
            V=results.V,
            T=T,
            pee=pee,
            pa=pa,
            b=b_val,
            theta=theta_vals,
        ),
    )
