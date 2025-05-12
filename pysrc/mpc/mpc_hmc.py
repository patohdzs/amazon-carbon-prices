import os
import time
import numpy as np
import pandas as pd
from pathlib import Path
from gams import GamsWorkspace
import shutil
from pysrc.services.file_service import get_path
from pysrc.services.data_service import load_productivity_params, load_site_data
from pysrc.services.file_service import get_path
from pysrc.mpc.mpc_optimization import mpc_solve_planner_problem
from pysrc.optimization import PlannerSolution

import argparse
parser = argparse.ArgumentParser(description="parameter settings")
parser.add_argument("--id",type=int,default=1)
parser.add_argument("--pe",type=float,default=20.76)
parser.add_argument("--xi",type=float,default=10000)
parser.add_argument("--trig",type=int,default=0)
parser.add_argument("--type",type=str,default="unconstrained")
args = parser.parse_args()

pe = args.pe
id = args.id
xi = args.xi
trig=args.trig
type = args.type

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


if trig==1:
    mode="converge"
    print("using worst case probability distribution")
else:
    mode=None

solver="gurobi"
num_sites=78
pa=41.11
model="mpc"

(
    zbar_2017,
    z_2017,
    forest_area_2017,
) = load_site_data(num_sites)

(theta_vals, gamma_vals) = load_productivity_params(num_sites)

x0_vals = gamma_vals * forest_area_2017

results = mpc_solve_planner_problem(
    time_horizon=200,
    theta=theta_vals,
    gamma=gamma_vals,
    x0=x0_vals,
    zbar=zbar_2017,
    z0=z_2017,
    price_emissions=pe,
    price_cattle=pa,
    solver=solver,
    id=id,
    forest_area_2017=forest_area_2017,
    xi=xi,
    mode=mode,
    price_low = price_low,
    price_high = price_high,
    prob_ll=prob_ll,
    prob_hh=prob_hh,
)
print("Results for pe = ", pe)


if trig==1:
    output_folder = (
        get_path("output")
        / "optimization"
        / "mpc_worstcase"
        / solver
        / f"{num_sites}sites"
        / f"xi_{xi}"
        / f"pa_{pa}"
        / f"pe_{pe}"
        / f"mc_{id}"
        / type
    )
else:
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