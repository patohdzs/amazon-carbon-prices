import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
from scipy.stats import gaussian_kde
from scipy.stats import entropy
import geopandas as gpd
import pandas as pd
from tqdm import tqdm
from pysrc.analysis.figures import land_allocation,density,trajectory_diff
from matplotlib.backends.backend_pdf import PdfPages
import sys
from pysrc.services.file_service import get_path

import argparse
parser = argparse.ArgumentParser(description="shadow price calculation")
parser.add_argument("--pee",type=float,default=5)
parser.add_argument("--xi",type=float,default=5)
parser.add_argument("--sites",type=int,default=78)

args = parser.parse_args()
pee=args.pee
xi=args.xi
num_sites=args.sites

solver="gurobi"
pa=41.11
pee_det=6.8


result_folder = os.path.join(
    str(get_path("output")),
    "sampling",
    solver,
    f"{num_sites}sites",
    f"pa_{pa}",
    f"xi_{xi}",
)
prior_folder = os.path.join(
    str(get_path("output")),
    "sampling",
    solver,
    f"{num_sites}sites",
    f"pa_{pa}",
    "xi_10000.0",
)

with open(result_folder + f"/pe_{pee}/results.pcl", "rb") as f:
    b0 = pickle.load(f)

with open(result_folder + f"/pe_{pee+15}/results.pcl", "rb") as f:
    b15 = pickle.load(f)

with open(prior_folder + f"/pe_{pee_det+15}/results.pcl", "rb") as f:
    results_unadjusted = pickle.load(f)
    
    
    
theta_unadjusted = results_unadjusted["final_sample"][:16000, :num_sites]
gamma_unadjusted = results_unadjusted["final_sample"][:16000, num_sites:]
theta_adjusted_b0 = b0["final_sample"][:16000, :num_sites]
gamma_adjusted_b0 = b0["final_sample"][:16000, num_sites:]
theta_adjusted_b15 = b15["final_sample"][:16000, :num_sites]
gamma_adjusted_b15 = b15["final_sample"][:16000, num_sites:]




def compute_kl(unadj, adj):
    common_grid = np.linspace(min(unadj.min(), adj.min()), max(unadj.max(), adj.max()), 100)
    p = gaussian_kde(unadj, bw_method='scott')(common_grid) + 1e-20
    q = gaussian_kde(adj, bw_method='scott')(common_grid) + 1e-20
    return entropy(p, q)

kl_data = {"id": np.arange(1, num_sites + 1)}

print("Computing theta KL divergences...")
for label, theta_hmc in zip(['theta_b0', 'theta_b15'], [theta_adjusted_b0, theta_adjusted_b15]):
    kl_data[label] = [
        compute_kl(theta_unadjusted[:, i], theta_hmc[:, i])
        for i in tqdm(range(num_sites), desc=label)
    ]

print("Computing gamma KL divergences...")
for label, gamma_hmc in zip(['gamma_b0', 'gamma_b15'], [gamma_adjusted_b0, gamma_adjusted_b15]):
    kl_data[label] = [
        compute_kl(gamma_unadjusted[:, i], gamma_hmc[:, i])
        for i in tqdm(range(num_sites), desc=label)
    ]

kl_df = pd.DataFrame(kl_data)

output_folder = str(get_path("output")) + f"/figures/entropy/site_{num_sites}/xi{xi}/"
os.makedirs(output_folder, exist_ok=True)
kl_df.to_csv(output_folder + 'kl_divergences_theta_gamma.csv', index=False)

# ------------------------
# Top-2 divergent sites per column
# ------------------------
top_kl_divergences = {}
for col in ['theta_b0', 'theta_b15', 'gamma_b0', 'gamma_b15']:
    top_sites = kl_df.nlargest(2, col)
    top_kl_divergences[col] = top_sites[['id', col]].values.tolist()

print("Top 2 KL divergences per parameter:")
for param, top in top_kl_divergences.items():
    for site, value in top:
        print(f"{param}: id {int(site)} → KL = {value:.4f}")
        
        
        
print("start plot densities")

density(num_sites=num_sites,pee=pee,xi=xi,solver=solver)
