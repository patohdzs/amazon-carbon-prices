import os
import pickle

import numpy as np
import pandas as pd

from pysrc.services.data_service import load_productivity_params
from pysrc.services.file_service import get_path


def format_float(value):
    return f"{value:.2f}"


def read_theta(num_sites):

    (theta_vals, gamma_vals) = load_productivity_params(num_sites)
    return theta_vals


def read_file(result_directory):
    
    Z = np.loadtxt(os.path.join(result_directory, "Z.txt"), delimiter=",")
    X = np.sum(np.loadtxt(os.path.join(result_directory, "X.txt"), delimiter=","),axis=1)
    Xdot = np.diff(X, axis=0)
    U = np.loadtxt(os.path.join(result_directory, "U.txt"), delimiter=",")
    V = np.loadtxt(os.path.join(result_directory, "V.txt"), delimiter=",")

    return (Z / 1e2, Xdot / 1e2, U / 1e2, V / 1e2)


def value_decom(pee=7.1, num_sites=78, solver="gurobi", pa=41.11, model="det", xi=1):
    b = [0, 10, 15, 20, 25]
    pe = [pee + bi for bi in b]
    kappa = 2.094215255
    zeta_u = 1.66e-4 * 1e11
    zeta_v = 1.00e-4 * 1e11
    output_folder = str(get_path("output")) + "/tables/"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if model == "det":
        dft_np = read_theta(num_sites)

    results = []
    for j in range(5):
        order = j

        result_folder = os.path.join(
            str(get_path("output")),
            "optimization",
            "det",
            solver,
            f"{num_sites}sites",
            f"pa_{pa}",
            f"pe_{pe[order]}",
        )

        if model == "hmc":
            theta_folder = os.path.join(
                str(get_path("output")),
                "sampling",
                solver,
                f"{num_sites}sites",
                f"pa_{pa}",
                f"xi_{xi}",
            )
            with open(theta_folder + f"/pe_{pe[order]}/results.pcl", "rb") as f:
                para_file = pickle.load(f)
            dft_np = para_file["final_sample"][:, :num_sites].mean(axis=0)

        (dfz_np, dfxdot, dfu_np, dfv_np) = read_file(result_folder)

        results_AO = []
        for i in range(200):
            result_AO = pa * np.dot(dfz_np[i + 1], dft_np) / ((1 + 0.02) ** (i))
            results_AO.append(result_AO)
        total_AO = np.sum(results_AO)

        results_NT = []
        for i in range(200):
            result_NT = (
                -b[order]
                * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i])
                / ((1 + 0.02) ** (i))
            )
            results_NT.append(result_NT)
        total_NT = np.sum(results_NT)

        results_CS = []
        for i in range(200):
            result_CS = (
                -pee * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i]) / ((1 + 0.02) ** (i))
            )
            results_CS.append(result_CS)
        total_CS = np.sum(results_CS)

        results_AC = []
        for i in range(200):
            result_AC = (
                ((zeta_u / 2)
                * (np.sum(dfu_np[i]) ) ** 2
                +
                (zeta_v / 2)
                * (np.sum(dfv_np[i]) ) ** 2
                )
                / ((1 + 0.02) ** (i))
            )
            results_AC.append(result_AC)
        total_AC = np.sum(results_AC)

        total_PV = total_AO + total_NT + total_CS - total_AC

        total_AO = total_AO.round(2)
        total_NT = total_NT.round(2)
        total_CS = total_CS.round(2)
        total_AC = total_AC.round(2)
        total_PV = total_PV.round(2)

        iteration_results = {
            "pa": format_float(pa),
            "pe": format_float(pe[order]),
            "b": b[order],
            "agricultural output value": format_float(total_AO),
            "net transfers": format_float(total_NT),
            "forest services": format_float(total_CS),
            "adjustment costs": format_float(total_AC),
            "planner value": format_float(total_PV),
        }

        results.append(iteration_results)

    results_df = pd.DataFrame(results)
    latex_code = results_df.to_latex(index=False)

    with open(
        output_folder + f"present_value_site{num_sites}_pa{pa}_{model}.tex", "w"
    ) as file:
        file.write(latex_code)

    return


def transfer_cost(pee=7.1, num_sites=78, solver="gurobi", pa=41.11, y=30, model="det",xi=5.0):
    b = [0, 10, 15, 20, 25]
    pe = [pee + bi for bi in b]



    if model=="det":
        baseline_folder = os.path.join(
            str(get_path("output")),
            "optimization",
            model,
            solver,
            f"{num_sites}sites",
            f"pa_{pa}",
            f"pe_{pe[0]}",
        )
        
    elif model=="hmc":
        baseline_folder = os.path.join(
            str(get_path("output")),
            "optimization",
            model,
            solver,
            f"{num_sites}sites",
            f"xi_{xi}",
            f"pa_{pa}",
            f"pe_{pe[0]}",
        )  
        
        
    kappa = 2.094215255
    output_folder = str(get_path("output")) + "/tables/"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    
    (dfz_np, dfxdot, dfu_np, dfv_np) = read_file(baseline_folder)

    results_NCE_base = []
    for i in range(y):
        result_NCE_base = -kappa * np.sum(dfz_np[i + 1]) + dfxdot[i]
        results_NCE_base.append(result_NCE_base)
    total_NCE_base = np.sum(results_NCE_base) * 100

    results_table = []
    for j in range(5):
        order = j

        if model=="hmc":
            result_folder = os.path.join(
                str(get_path("output")),
                "optimization",
                model,
                solver,
                f"{num_sites}sites",
                f"xi_{xi}",
                f"pa_{pa}",
                f"pe_{pe[order]}",
            )
        elif model=="det":
                result_folder = os.path.join(
                str(get_path("output")),
                "optimization",
                model,
                solver,
                f"{num_sites}sites",
                f"pa_{pa}",
                f"pe_{pe[order]}",
            )

        (dfz_np, dfxdot, _, _) = read_file(result_folder)

        results_NCE = []
        for i in range(y):
            result_NCE = -kappa * np.sum(dfz_np[i + 1]) + dfxdot[i]
            results_NCE.append(result_NCE)
        total_NCE = np.sum(results_NCE) * 100

        results_NT2 = []
        for i in range(y):
            result_NT2 = (
                -b[order]
                * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i])
                / ((1 + 0.02) ** (i))
            )

            results_NT2.append(result_NT2)
        total_NT2 = np.sum(results_NT2)

        total_EC = b[order]*total_NCE / (total_NCE - total_NCE_base) 

        total_NCE = total_NCE.round(2)
        total_NT2 = total_NT2.round(2)
        total_EC = total_EC.round(2)
        iteration_results = {
            "pe": format_float(pe[order]),
            "b": b[order],
            "net captured emissions": format_float(total_NCE),
            "discounted net transfers": format_float(total_NT2),
            "discounted effective costs": format_float(total_EC),
        }

        results_table.append(iteration_results)

    results_df = pd.DataFrame(results_table)
    latex_code = results_df.to_latex(index=False)


    if model =="det":
        with open(
            output_folder + f"transfer_cost_{num_sites}site_{pa}pa_{y}year_{model}.tex", "w"
        ) as file:
            file.write(latex_code)
    elif model =="hmc":
        with open(
            output_folder + f"transfer_cost_{num_sites}site_{pa}pa_{y}year_{model}_xi_{xi}.tex", "w"
        ) as file:
            file.write(latex_code)            
        
        
    return


def ambiguity_decom(pe_det=7.1, pe_hmc=5.3, num_sites=78, solver="gurobi", pa=41.11, xi=1):
    kappa = 2.094215255
    zeta_u = 1.66e-4 * 1e11
    zeta_v = 1.00e-4 * 1e11
    output_folder = str(get_path("output")) + "/tables/"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    results_det = []
    pee = pe_det
    b = [0, 10, 15, 20, 25]
    pe = [pee + bi for bi in b]
    dft_np = read_theta(num_sites)
    for j in range(5):
        order = j
        result_folder = os.path.join(
            str(get_path("output")),
            "optimization",
            "det",
            solver,
            f"{num_sites}sites",
            f"pa_{pa}",
            f"pe_{pe[order]}",
        )

        (dfz_np, dfxdot, dfu_np, dfv_np) = read_file(result_folder)

        results_AO = []
        for i in range(200):
            result_AO = pa * np.dot(dfz_np[i + 1], dft_np) / ((1 + 0.02) ** (i))
            results_AO.append(result_AO)
        total_AO = np.sum(results_AO)

        results_NT = []
        for i in range(200):
            result_NT = (
                -b[order]
                * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i])
                / ((1 + 0.02) ** (i))
            )
            results_NT.append(result_NT)
        total_NT = np.sum(results_NT)

        results_CS = []
        for i in range(200):
            result_CS = (
                -pee * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i]) / ((1 + 0.02) ** (i))
            )
            results_CS.append(result_CS)
        total_CS = np.sum(results_CS)

        results_AC = []
        for i in range(200):
            result_AC = (
                ((zeta_u / 2)
                * (np.sum(dfu_np[i]) ) ** 2
                +
                (zeta_v / 2)
                * (np.sum(dfv_np[i]) ) ** 2
                )
                / ((1 + 0.02) ** (i))
            )
            results_AC.append(result_AC)
        total_AC = np.sum(results_AC)

        total_PV = total_AO + total_NT + total_CS - total_AC
        
        total_AO = total_AO.round(2)
        total_NT = total_NT.round(2)
        total_CS = total_CS.round(2)
        total_AC = total_AC.round(2)
        total_PV = total_PV.round(2)

        results_det.append(
            {"b": round(b[order]), "total_AO": total_AO, "total_PV": total_PV}
        )

    results_det_df = pd.DataFrame(results_det)

    results_hmc = []
    pee = pe_hmc
    b = [0, 10, 15, 20, 25]
    pe = [pee + bi for bi in b]
    for j in range(5):
        order = j

        result_folder = os.path.join(
            str(get_path("output")),
            "optimization",
            "hmc",
            solver,
            f"{num_sites}sites",
            f"xi_{xi}",
            f"pa_{pa}",
            f"pe_{pe[order]}",
        )

        theta_folder = os.path.join(
            str(get_path("output")),
            "sampling",
            solver,
            f"{num_sites}sites",
            f"pa_{pa}",
            f"xi_{xi}",
        )
        with open(theta_folder + f"/pe_{pe[order]}/results.pcl", "rb") as f:
            para_file = pickle.load(f)
        dft_np = para_file["final_sample"][:16000, :num_sites].mean(axis=0)

        (dfz_np, dfxdot, dfu_np, dfv_np) = read_file(result_folder)

        results_AO = []
        for i in range(200):
            result_AO = pa * np.dot(dfz_np[i + 1], dft_np) / ((1 + 0.02) ** (i))
            results_AO.append(result_AO)
        total_AO = np.sum(results_AO)

        results_NT = []
        for i in range(200):
            result_NT = (
                -b[order]
                * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i])
                / ((1 + 0.02) ** (i))
            )
            results_NT.append(result_NT)
        total_NT = np.sum(results_NT)

        results_CS = []
        for i in range(200):
            result_CS = (
                -pee * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i]) / ((1 + 0.02) ** (i))
            )
            results_CS.append(result_CS)
        total_CS = np.sum(results_CS)

        results_AC = []
        for i in range(200):
            result_AC = (
                ((zeta_u / 2)
                * (np.sum(dfu_np[i]) ) ** 2
                +
                (zeta_v / 2)
                * (np.sum(dfv_np[i]) ) ** 2
                )
                / ((1 + 0.02) ** (i))
            )
            results_AC.append(result_AC)
        total_AC = np.sum(results_AC)

        total_PV = total_AO + total_NT + total_CS - total_AC
        total_AO = total_AO.round(2)
        total_NT = total_NT.round(2)
        total_CS = total_CS.round(2)
        total_AC = total_AC.round(2)
        total_PV = total_PV.round(2)

        results_hmc.append(
            {"b": round(b[order]), "total_AO": total_AO, "total_PV": total_PV}
        )

    results_hmc_df = pd.DataFrame(results_hmc)

    combined_df = pd.DataFrame(
        {
            "b": results_det_df["b"],
            "AO_det": results_det_df["total_AO"],
            "AO_hmc": results_hmc_df["total_AO"],
            "PV_det": results_det_df["total_PV"],
            "PV_hmc": results_hmc_df["total_PV"],
        }
    )

    # Calculate the percentage change for AO and PV
    combined_df["% Change AO"] = ((
        (combined_df["AO_hmc"] - combined_df["AO_det"]) / combined_df["AO_det"]
    ) * 100).round(2)
    combined_df["% Change PV"] =((
        (combined_df["PV_hmc"] - combined_df["PV_det"]) / combined_df["PV_det"]
    ) * 100).round(2)

    combined_df = combined_df[["b","AO_det", "AO_hmc", "% Change AO", "PV_det", "PV_hmc", "% Change PV"]]

    # Display the final DataFrame
    print(combined_df)

    latex_code = combined_df.applymap(lambda x: f"{x:.2f}").to_latex(index=False)

    with open(
        output_folder + f"present_value_site_ambiguity_comparison_xi_{xi}.tex", "w"
    ) as file:
        file.write(latex_code)

    return


def value_decom_mpc(pee=6.9, num_sites=78, solver="gams", model="unconstrained", b=0):
    pe = pee + b
    kappa = 2.094215255
    zeta_u = 1.66e-4 * 1e11
    zeta_v = 1.00e-4 * 1e11
    dft_np = read_theta(num_sites)
    df_ori = pd.read_csv(str(get_path("data")) + "/calibration/hmc/calibration_78_sites.csv")
    x0 = df_ori["x_2017_78Sites"].to_numpy()

    output_folder = str(get_path("output")) + "/tables/"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if solver == "gams":

        def gams_read_dat_file(filename):
            data_dict = {}

            with open(filename, "r") as file:
                lines = file.readlines()

            current_y = None
            values = []
            p_a = None

            for line in lines:
                line = line.strip()

                if line.startswith("Y ="):
                    if current_y is not None:
                        data_dict[current_y] = {"p_a": p_a, "values": values}
                        values = []

                    parts = line.split()
                    current_y = int(parts[2])
                    p_a = float(parts[5])

                else:
                    try:
                        value = float(line)
                        values.append(value)
                    except ValueError:
                        pass

            if current_y is not None:
                data_dict[current_y] = {"p_a": p_a, "values": values}

            return data_dict

        def gams_read_file(result_directory):
            dfz = gams_read_dat_file(result_directory + "/amazon_data_z.dat")
            dfz_np = (
                np.array([np.array(dfz[i]["values"]) for i in range(1, 201)]) / 1e11
            )
            p_a_values = np.array([dfz[i]["p_a"] for i in range(1, 201)])

            dfx = gams_read_dat_file(result_directory + "/amazon_data_x.dat")
            dfx_np = (
                np.array([np.array(dfx[i]["values"]) for i in range(1, 201)]) / 1e11
            )
            dfx_np = np.sum(dfx_np, axis=1)
            dfxdot = np.diff(dfx_np, axis=0)

            dfu = gams_read_dat_file(result_directory + "/amazon_data_u.dat")
            dfu_np = (
                np.array([np.array(dfu[i]["values"]) for i in range(1, 201)]) / 1e11
            )

            dfv = gams_read_dat_file(result_directory + "/amazon_data_v.dat")
            dfv_np = (
                np.array([np.array(dfv[i]["values"]) for i in range(1, 201)]) / 1e11
            )

            dfz_np = np.vstack((dfz_np, dfu_np[-1] - dfv_np[-1]))

            return (dfz_np, p_a_values, dfxdot, dfu_np, dfv_np, dfx_np)

        results = []
        for j in range(99):
            result_directory = (
                str(get_path("output"))
                + f"/optimization/mpc/gams/78sites/model_{model}/"
                + f"pe_{pe}/mc_{j+1}"
            )

            (dfz_np, p_a_values, dfxdot, dfu_np, dfv_np, dfx_np) = gams_read_file(
                result_directory
            )

            x_det = np.insert(dfx_np, 0, np.sum(x0) / 1e11)
            dfxdot = np.diff(x_det, axis=0)

            results_AO = []
            for i in range(200):
                result_AO = (
                    p_a_values[i] * np.dot(dfz_np[i + 1], dft_np) / ((1 + 0.02) ** (i))
                )
                results_AO.append(result_AO)
            total_AO = np.sum(results_AO)

            results_NT = []
            for i in range(200):
                result_NT = (
                    -b
                    * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i])
                    / ((1 + 0.02) ** (i))
                )
                results_NT.append(result_NT)
            total_NT = np.sum(results_NT)

            results_CS = []
            for i in range(200):
                result_CS = (
                    -pee
                    * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i])
                    / ((1 + 0.02) ** (i))
                )
                results_CS.append(result_CS)
            total_CS = np.sum(results_CS)

            results_AC = []
            for i in range(200):
                result_AC = (
                    ((zeta_u / 2)
                    * (np.sum(dfu_np[i]) ) ** 2
                    +
                    (zeta_v / 2)
                    * (np.sum(dfv_np[i]) ) ** 2
                    )
                    / ((1 + 0.02) ** (i))
                )
                results_AC.append(result_AC)
            total_AC = np.sum(results_AC)

            total_PV = total_AO + total_NT + total_CS - total_AC

            iteration_results = {
                "j": j + 1,
                "b": b,
                "total_AO": total_AO,
                "total_NT": total_NT,
                "total_CS": total_CS,
                "total_AC": total_AC,
                "total_PV": total_PV,
            }

            results.append(iteration_results)

        results_df = pd.DataFrame(results)

        mean = results_df.mean()
        sd = results_df.std()
        sd / mean
        quantiles_10 = results_df.quantile(0.10)
        quantiles_50 = results_df.quantile(0.50)
        quantiles_90 = results_df.quantile(0.90)

        # Create a new DataFrame for the summary table
        summary_table_df = pd.DataFrame(
            {
                "": [r"10\%", r"50\%", r"90\%"],
                "agricultural output value": [
                    format_float(quantiles_10["total_AO"]),
                    format_float(quantiles_50["total_AO"]),
                    format_float(quantiles_90["total_AO"]),
                ],
                "net transfers": [
                    format_float(quantiles_10["total_NT"]),
                    format_float(quantiles_50["total_NT"]),
                    format_float(quantiles_90["total_NT"]),
                ],
                "forest services": [
                    format_float(quantiles_10["total_CS"]),
                    format_float(quantiles_50["total_CS"]),
                    format_float(quantiles_90["total_CS"]),
                ],
                "adjustment costs": [
                    format_float(quantiles_10["total_AC"]),
                    format_float(quantiles_50["total_AC"]),
                    format_float(quantiles_90["total_AC"]),
                ],
                "planner value": [
                    format_float(quantiles_10["total_PV"]),
                    format_float(quantiles_50["total_PV"]),
                    format_float(quantiles_90["total_PV"]),
                ],
            }
        )

    with open(output_folder + f"present_value_mpc_b{b}.tex", "w") as file:
        file.write(summary_table_df.to_latex(index=False))

    return print("done")
