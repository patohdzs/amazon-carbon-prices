import os
import pickle
import numpy as np
from pysrc.analysis import value_decomposition
from pysrc.optimization import solve_planner_problem
from pysrc.sampling import adjusted
from pysrc.services.data_service import load_site_data
from pysrc.services.file_service import get_path

# ## Model scenario
solver = "gurobi"  # need to install gurobi solver
pee = 4.5
pa = 41.11
num_sites = 78
T = 200
b = [0, 10, 15, 20, 25]
pe_values = [pee + bi for bi in b]
xi = 0.5
pe=pee

# inits = {
#     "log_precision_u_gamma": 4,
#     "log_precision_v_gamma": 4,
#     "log_precision_u_theta": 4,
#     "log_precision_v_theta": 4,
#     "beta_gamma": [0.2] * 6,
#     "beta_theta": [0.2] * 8,
#     "nu_gamma_transform": [0.2] * 78,
#     "nu_theta_transform": [0.2] * 62,
# }


results = adjusted.sample(
    xi=xi,
    pe=pe,
    pa=pa,
    weight=0.20,
    num_sites=num_sites,
    T=200,
    solver=solver,
    max_iter=100,
    final_sample_size=2_000,
    iter_sampling=2000,
    iter_warmup=1000,
    show_progress=True,
    chains=4,
    seed=1,
    inits=0.2,
    adapt_delta=0.99,
    # inits=inits,
    zeta_u=1.66e-4 * 1e9,
    zeta_v=1.00e-4 * 1e9,
    norm_fac=1e9,
)

output_base_path = os.path.join(
    str(get_path("output")),
    "sampling",
    solver,
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