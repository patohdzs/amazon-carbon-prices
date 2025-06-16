import os
import pickle
import pandas as pd
import numpy as np
from pysrc.analysis import value_decomposition
from pysrc.optimization import solve_planner_problem
from pysrc.sampling import adjusted
from pysrc.services.data_service import load_site_data
from pysrc.services.file_service import get_path
import argparse


parser = argparse.ArgumentParser(description="test")
parser.add_argument("--id",type=int,default=1)

args = parser.parse_args()
id=args.id


# ## Model scenario
solver = "gurobi"  # need to install gurobi solver
pee = 2.8
pa = 41.11
num_sites = 1043
T = 200
b = [0, 10, 15, 20, 25]
pe_values = [pee + bi for bi in b]
xi = 0.5
pe=pee


results = adjusted.sample(
    xi=xi,
    pe=pe,
    pa=pa,
    weight=0.25,
    num_sites=num_sites,
    T=200,
    solver=solver,
    max_iter=100,
    final_sample_size=4_000,
    iter_sampling=4000,
    iter_warmup=500,
    show_progress=True,
    seed=id,
    # inits=0.1,
    chains=4,  
    tol=0.005,
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