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


def value_decom_mpc(pee=5.9, num_sites=78, solver="gurobi", model="unconstrained", b=0,xi=10000,mode=None,price_low=35.76,price_high = 44.25):
    pe = pee + b
    kappa = 2.094215255
    zeta_u = 1.66e-4 * 1e11
    zeta_v = 1.00e-4 * 1e11

    
    dft_np = read_theta(num_sites)

    output_folder = str(get_path("output")) + "/mpc/"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    results = []
    for j in range(49):
        
        
        if mode =="converge":
            if price_low == 32.49:
                pa_file_path = (
                    get_path("output", "simulation", "mpc_path","constrained",f"xi_{xi}",f"pe_{pe}") / f"mc_{j+1}.csv"
                ) 
            else:
                pa_file_path = (
                    get_path("output", "simulation", "mpc_path","unconstrained",f"xi_{xi}",f"pe_{pe}") / f"mc_{j+1}.csv"
                ) 
        else:
            pa_file_path = (
                get_path("output", "simulation", "mpc_path","baseline",f"{model}") / f"mc_{j+1}.csv"
            ) 

        p_a_values = np.array(pd.read_csv(pa_file_path))[:,1]
        p_a_values = np.where(p_a_values == 2, price_high, price_low)
        
        if mode =="converge":
            result_directory = (
                str(get_path("output"))
                + f"/optimization/mpc_worstcase/{solver}/{num_sites}sites/xi_{xi}/pa_41.11/"
                + f"pe_{pe}/mc_{j+1}/{model}"
            )            
            
        else:
            result_directory = (
                str(get_path("output"))
                + f"/optimization/mpc/{solver}/{num_sites}sites/xi_{xi}/pa_41.11/"
                + f"pe_{pe}/mc_{j+1}/{model}"
            )

        (dfz_np, dfxdot, dfu_np, dfv_np) = read_file(
            result_directory
        )


        results_AO = []
        for i in range(200):
            result_AO = p_a_values[i] * np.dot(dfz_np[i + 1], dft_np) / ((1 + 0.02) ** (i))
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

    summary_table_df = pd.DataFrame(
        {
            "  ":  f"b = {b}",
            "agricultural output value": [format_float(mean["total_AO"])],
            "net transfers": [format_float(mean["total_NT"])],
            "forest services": [format_float(mean["total_CS"])],
            "adjustment costs": [format_float(mean["total_AC"])],
            "planner value": [format_float(mean["total_PV"])],
        }
    )


    if mode =="converge":
        with open(output_folder + f"converge_present_value_mpc_b{b}_sites{num_sites}_xi_{xi}_pee_{pee}_{model}.tex", "w") as file:
            file.write(summary_table_df.to_latex(index=False))
    else:
        with open(output_folder + f"present_value_mpc_b{b}_sites{num_sites}_xi_{xi}_pee_{pee}_{model}.tex", "w") as file:
            file.write(summary_table_df.to_latex(index=False))

    return print("done")


# for b in [0,10,15,25]:
#     value_decom_mpc(pee=5.5,num_sites=78,b=b,xi=0.2)
    
# for b in [0,10,15,25]:
#     value_decom_mpc(pee=6.3,num_sites=78,b=b,xi=10000.0)
    
# for b in [0,10,15,25]:
#     value_decom_mpc(pee=5.7,num_sites=78,b=b,xi=1.0,mode="converge")
    
# for b in [0,10,15,25]:
#     value_decom_mpc(pee=5.5,num_sites=78,b=b,xi=0.2,mode="converge")



######### Appendix tables

# for b in [0,10,15,20,25]:
#     value_decom_mpc(pee=5.9,num_sites=78,b=b,xi=10000.0,model="constrained",price_low=32.49,price_high = 42.85)
    
# for b in [0,10,15,20,25]:
#     value_decom_mpc(pee=5.1,num_sites=78,b=b,xi=1.0,mode="converge",model="constrained",price_low=32.49,price_high = 42.85)
    
# for b in [0,10,15,20,25]:
#     value_decom_mpc(pee=5.0,num_sites=78,b=b,xi=0.2,mode="converge",model="constrained",price_low=32.49,price_high = 42.85)
    
    
    
    
    
    
    
    
    
def transfer_cost_mpc(pee=5.9, y=30,num_sites=78, solver="gurobi", model="unconstrained", b=0,xi=10000,mode=None,price_low=35.76,price_high = 44.25):
    kappa = 2.094215255
    pe = pee + b

    output_folder = str(get_path("output")) + "/mpc/"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    results = []
    for j in range(49):
        
        if mode =="converge":
            result_directory = (
                str(get_path("output"))
                + f"/optimization/mpc_worstcase/{solver}/{num_sites}sites/xi_{xi}/pa_41.11/"
                + f"pe_{pe}/mc_{j+1}/{model}"
            )            
            baseline_folder = (
                str(get_path("output"))
                + f"/optimization/mpc_worstcase/{solver}/{num_sites}sites/xi_{xi}/pa_41.11/"
                + f"pe_{pee}/mc_{j+1}/{model}"
            )            
        else:
            result_directory = (
                str(get_path("output"))
                + f"/optimization/mpc/{solver}/{num_sites}sites/xi_{xi}/pa_41.11/"
                + f"pe_{pe}/mc_{j+1}/{model}"
            )
            baseline_folder = (
                str(get_path("output"))
                + f"/optimization/mpc/{solver}/{num_sites}sites/xi_{xi}/pa_41.11/"
                + f"pe_{pee}/mc_{j+1}/{model}"
            )

        (dfz_np, dfxdot, dfu_np, dfv_np) = read_file(baseline_folder)

        results_NCE_base = []
        for i in range(y):
            result_NCE_base = -kappa * np.sum(dfz_np[i + 1]) + dfxdot[i]
            results_NCE_base.append(result_NCE_base)
        total_NCE_base = np.sum(results_NCE_base) * 100

        (dfz_np, dfxdot, dfu_np, dfv_np) = read_file(
            result_directory
        )
    

        results_NCE = []
        for i in range(y):
            result_NCE = -kappa * np.sum(dfz_np[i + 1]) + dfxdot[i]
            results_NCE.append(result_NCE)
        total_NCE = np.sum(results_NCE) * 100

        results_NT2 = []
        for i in range(y):
            result_NT2 = (
                -b
                * (kappa * np.sum(dfz_np[i + 1]) - dfxdot[i])
                / ((1 + 0.02) ** (i))
            )

            results_NT2.append(result_NT2)
        total_NT2 = np.sum(results_NT2)

        total_EC = total_NT2 / (total_NCE - total_NCE_base) * 100

        iteration_results = {
            "j": j + 1,
            "b": b,
            "NCE": total_NCE,
            "NT2": total_NT2,
            "EC": total_EC,
        }

        results.append(iteration_results)

    results_df = pd.DataFrame(results)

    mean = results_df.mean()
    sd = results_df.std()
    sd / mean

    summary_table_df = pd.DataFrame(
        {
            "  ":  f"b = {b}",
            "net captured emissions": [format_float(mean["NCE"])],
            "discounted net transfers": [format_float(mean["NT2"])],
            "discounted effective costs": [format_float(mean["EC"])],
        }
    )


    if mode =="converge":
        with open(output_folder + f"converge_transfer_mpc_b{b}_sites{num_sites}_xi_{xi}_pee_{pee}_{model}.tex", "w") as file:
            file.write(summary_table_df.to_latex(index=False))
    else:
        with open(output_folder + f"transfer_mpc_b{b}_sites{num_sites}_xi_{xi}_pee_{pee}_{model}.tex", "w") as file:
            file.write(summary_table_df.to_latex(index=False))

    return print("done")


for b in [0,10,15,25]:
    transfer_cost_mpc(pee=6.3,num_sites=78,b=b,xi=10000.0)
    
for b in [0,10,15,25]:
    transfer_cost_mpc(pee=5.7,num_sites=78,b=b,xi=1.0,mode="converge")
    
for b in [0,10,15,25]:
    transfer_cost_mpc(pee=5.5,num_sites=78,b=b,xi=0.2,mode="converge")
