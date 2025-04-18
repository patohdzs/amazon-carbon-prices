import numpy as np
from pysrc.services.get_sample import get_sampling
import argparse
import os
import pickle

from pysrc.sampling import adjusted
from pysrc.services.file_service import get_path

parser = argparse.ArgumentParser(description="shadow price calculation")
parser.add_argument("--pee",type=float,default=5)
parser.add_argument("--xi",type=float,default=5)
parser.add_argument("--id",type=int,default=0)
parser.add_argument("--sites",type=int,default=78)

args = parser.parse_args()
pee=args.pee
xi=args.xi
b=args.id
num_sites=args.sites

if xi==1:
    pee=2.2

pe = pee + b
solver="gams"
pa=41.11

results = adjusted.sample(
    xi=xi,
    pe=pe,
    pa=pa,
    weight=0.25,
    num_sites=num_sites,
    T=200,
    solver=solver,
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