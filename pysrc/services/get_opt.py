import pickle
from pathlib import Path
import os
import numpy as np

from pysrc.optimization import PlannerSolution, solve_planner_problem
from pysrc.services.data_service import load_productivity_params, load_site_data
from pysrc.services.file_service import get_path


def first_year_total_xdot_decreasing(results: PlannerSolution, tol: float = 1e-10):
    total_x = np.sum(results.X, axis=1)
    total_x_dot = np.diff(total_x)
    total_x_dot_nonincreasing = np.diff(total_x_dot) <= tol
    suffix_nonincreasing = np.logical_and.accumulate(
        total_x_dot_nonincreasing[::-1]
    )[::-1]
    valid_years = np.flatnonzero(suffix_nonincreasing)
    if valid_years.size == 0:
        return None
    return int(valid_years[0])


def print_terminal_diagnostics(results: PlannerSolution, b, pe):
    max_z_100 = np.max(results.Z[100])
    max_z_100_site = int(np.argmax(results.Z[100]) + 1)
    print(
        "Max z_i at t=100 "
        f"(b={b}, pe={pe}) = {max_z_100:.15g} billion ha "
        f"({max_z_100 * 1e9:.15g} ha), "
        f"site={max_z_100_site}"
    )

    first_decreasing_year = first_year_total_xdot_decreasing(results)
    if first_decreasing_year is None:
        print(
            "First year after which total x_dot is decreasing "
            f"(b={b}, pe={pe}) = NA"
        )
    else:
        print(
            "First year after which total x_dot is decreasing "
            f"(b={b}, pe={pe}) = {first_decreasing_year}"
        )


def get_optimization(
    solver="gurobi",
    num_sites=78,
    pee=5,
    pa=41.11,
    model="det",
    xi=1,
):
    (
        zbar_2017,
        z_2017,
        forest_area_2017,
    ) = load_site_data(num_sites)

    if model == "det":
        # Load baseline productivity params
        (theta_vals, gamma_vals) = load_productivity_params(num_sites)

        # Compute initial carbon stock
        x0_vals = gamma_vals * forest_area_2017

        # Solve for all transfer levels
        b = [0, 10, 15, 20, 25]
        pe_values = [pee + bi for bi in b]
        for bi, pe in zip(b, pe_values):
            results = solve_planner_problem(
                time_horizon=200,
                theta=theta_vals,
                gamma=gamma_vals,
                x0=x0_vals,
                zbar=zbar_2017,
                z0=z_2017,
                price_emissions=pe,
                price_cattle=pa,
                solver=solver,
            )

            print("Results for pe = ", pe)
            print_terminal_diagnostics(results, bi, pe)
            output_folder = (
                get_path("output")
                / "optimization"
                / model
                / solver
                / f"{num_sites}sites"
                / f"pa_{pa}"
                / f"pe_{pe}"
            )

            save_planner_solution(results, output_folder)

    if model == "hmc":
        results_dir = os.path.join(
            get_path("output")
            / "sampling"
            / solver
            / f"{num_sites}sites"
            / f"pa_{pa}"
            / f"xi_{xi}"
        )

        # Solve model for all transfer values
        b = [0, 10, 15, 20, 25]
        pe_values = [pee + bi for bi in b]
        for bi, pe in zip(b, pe_values):
            # Load ambiguity-adjusted params
            with open(results_dir + f"/pe_{pe}/results.pcl", "rb") as f:
                b = pickle.load(f)

            # Take sample
            theta_vals = b["final_sample"][:16000, :num_sites].mean(axis=0)
            gamma_vals = b["final_sample"][:16000, num_sites:].mean(axis=0)

            # Compute carbon stock
            x0_vals = gamma_vals * forest_area_2017

            results = solve_planner_problem(
                time_horizon=200,
                theta=theta_vals,
                gamma=gamma_vals,
                x0=x0_vals,
                zbar=zbar_2017,
                z0=z_2017,
                price_emissions=pe,
                price_cattle=pa,
                solver=solver,
            )
            print("Results for pe = ", pe)
            print_terminal_diagnostics(results, bi, pe)

            output_folder = (
                get_path("output")
                / "optimization"
                / model
                / solver
                / f"{num_sites}sites"
                / f"xi_{xi}"
                / f"pa_{pa}"
                / f"pe_{pe}"
            )
            

            save_planner_solution(results, output_folder)
    return


def save_planner_solution(results: PlannerSolution, output_dir: Path):
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
                
    np.savetxt(output_dir / "Z.txt", results.Z, delimiter=",")
    np.savetxt(output_dir / "X.txt", results.X, delimiter=",")
    np.savetxt(output_dir / "U.txt", results.U, delimiter=",")
    np.savetxt(output_dir / "V.txt", results.V, delimiter=",")
