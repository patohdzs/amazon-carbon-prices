import os
import numpy as np
from pathlib import Path
from pysrc.services.file_service import get_path
from pysrc.services.data_service import load_productivity_params
from pysrc.services.file_service import get_path
from pysrc.mpc.mpc_optimization import mpc_solve_planner_problem
from pysrc.optimization import PlannerSolution
from pysrc.services.data_service import load_site_data_1995

import argparse
parser = argparse.ArgumentParser(description="parameter settings")
parser.add_argument("--pe",type=float,default=20.76)
parser.add_argument("--xi",type=float,default=10000)
parser.add_argument("--type",type=str,default="unconstrained")
args = parser.parse_args()
type = args.type

pe = args.pe
id = 999
xi = args.xi



if type =="constrained":

    price_low = 32.49
    price_high = 42.85
    prob_ll=0.762
    prob_hh=0.959

if type =="unconstrained":

    price_low = 35.76
    price_high = 44.32
    prob_ll=0.707
    prob_hh=0.826




solver="gurobi"
num_sites=78
pa=41.11
model="mpc_shadow_price"

(
    zbar_1995,
    z_1995,
    forest_area_1995,
    z_2008,
    theta,
    gamma,
) = load_site_data_1995(num_sites)

(theta_vals, gamma_vals) = load_productivity_params(num_sites)

x0_vals = gamma_vals * forest_area_1995

results = mpc_solve_planner_problem(
    time_horizon=200,
    theta=theta_vals,
    gamma=gamma_vals,
    x0=x0_vals,
    zbar=zbar_1995,
    z0=z_1995,
    price_emissions=pe,
    price_cattle=pa,
    solver=solver,
    id=id,
    forest_area_2017=forest_area_1995,
    xi=xi,
    mode="sp",
    price_low = price_low,
    price_high = price_high,
    prob_ll=prob_ll,
    prob_hh=prob_hh,
    type=type,
)
print("Results for pe = ", pe)

output_folder = (
    get_path("output")
    / "optimization"
    / model
    / solver
    / f"{num_sites}sites"
    / f"xi_{xi}"
    / f"pa_{pa}"
    / f"pe_{pe}"
    / f"mc_{id}"
    / type
)


def save_planner_solution(results: PlannerSolution, output_dir: Path):
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
                
    np.savetxt(output_dir / "Z.txt", results.Z, delimiter=",")
    np.savetxt(output_dir / "X.txt", results.X, delimiter=",")
    np.savetxt(output_dir / "U.txt", results.U, delimiter=",")
    np.savetxt(output_dir / "V.txt", results.V, delimiter=",")


save_planner_solution(results, output_folder)


print("all done")