import os
import pickle

import numpy as np
import pandas as pd

from pysrc.services.data_service import load_productivity_params
from pysrc.services.file_service import get_path
from pysrc.services.data_service import load_site_data_1995,load_price_data

def format_float(value):
    return f"{value:.3f}"


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


def mpc_sp(pee=5.9, num_sites=78, solver="gurobi", model="unconstrained", b=0,xi=10000):
    pe = pee + b
    mc_id=999
    
    
    (
        zbar_1995,
        z_1995,
        forest_area_1995,
        z_2008,
        theta,
        gamma,
    ) = load_site_data_1995(num_sites)


        
    result_directory = (
        str(get_path("output"))
        + f"/optimization/mpc_shadow_price/{solver}/{num_sites}sites/xi_{xi}/pa_41.11/"
        + f"pe_{pe}/mc_{mc_id}/"+f"{model}"
    )

    (dfz_np, dfxdot, dfu_np, dfv_np) = read_file(
        result_directory
    )


    Z = dfz_np
    z_2008_agg = np.sum(z_2008) / 1e11
    ratio = (np.sum(Z[13]) - z_2008_agg) / z_2008_agg

    return ratio






for xi in (0.5, 1.0, 10000.0):
    for pe in np.linspace(5.0, 6.3, num=14):
        print("xi",xi,"pe",pe,"model unconstrained")
        print("ratio",mpc_sp(pee=pe,num_sites=78,b=0,xi=xi,model="unconstrained"))

for xi in (0.5, 1.0, 10000.0):
    for pe in np.linspace(5.0, 6.3, num=14):
        print("xi",xi,"pe",pe,"model constrained")
        print("ratio",mpc_sp(pee=pe,num_sites=78,b=0,xi=xi,model="constrained"))





#### unconstrained model
#### 6.3 for xi=10000
#### 6.0 for xi=1
#### 5.6 for xi=0.5


#### constrained model
#### 6.0 for xi=10000
#### 5.7 for xi=1
#### 5.2 for xi=0.5