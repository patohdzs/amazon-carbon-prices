import numpy as np

from pysrc.optimization import solve_planner_problem
from pysrc.sampling import adjusted
from pysrc.services.data_service import load_site_data_1995,load_price_data
import argparse

parser = argparse.ArgumentParser(description="shadow price calculation")
parser.add_argument("--id",type=int,default=400)
parser.add_argument("--xi",type=float,default=5)
parser.add_argument("--sites",type=int,default=1043)
args = parser.parse_args()
seed = args.id
seed_i = seed/10
xi=args.xi
num_sites=args.sites


def shadow_price_opt(
    zbar_1995,
    z_1995,
    forest_area_1995,
    z_2008,
    theta,
    gamma,
    sitenum=78,
    solver="gurobi",
    timehzn=200,
    pa=41.11,
    pe=7.1,
    model="det",
):
    
    pa_list = load_price_data()
    price_cattle = np.concatenate((pa_list, np.full(200 - len(pa_list), pa)))
    
    # Computing carbon absorbed in start period
    x0_vals_1995 = gamma * forest_area_1995

    # if model == "mpc":
    #     solve_planner_problem = gams.mpc_shadow_price

    results = solve_planner_problem(
        theta=theta,
        gamma=gamma,
        x0=x0_vals_1995,
        zbar=zbar_1995,
        z0=z_1995,
        price_emissions=pe,
        price_cattle=price_cattle,
        solver=solver,
    )
    Z = results.Z
    z_2008_agg = np.sum(z_2008) / 1e9
    ratio = (np.sum(Z[13]) - z_2008_agg) / z_2008_agg

    return ratio


def shadow_price_cal(sitenum=78, pa=41.11, solver="gurobi", model="det", xi=2, pe_low=5, pe_high=8):
    if model == "det":
        (
            zbar_1995,
            z_1995,
            forest_area_1995,
            z_2008,
            theta,
            gamma,
        ) = load_site_data_1995(sitenum)

        pe_values = np.arange(pe_low, pe_high, 0.1)

        results = np.array(
            [
                shadow_price_opt(
                    zbar_1995,
                    z_1995,
                    forest_area_1995,
                    z_2008,
                    theta,
                    gamma,
                    sitenum=sitenum,
                    solver=solver,
                    timehzn=200,
                    pa=pa,
                    pe=pe,
                )
                for pe in pe_values
            ]
        )

    elif model == "hmc":
        (
            zbar_1995,
            z_1995,
            forest_area_1995,
            z_2008,
            theta,
            gamma,
        ) = load_site_data_1995(sitenum)

        pe_values = np.arange(pe_low, pe_high, 0.1)
        results = []
        for pe in pe_values:
            samples = adjusted.sample(
                xi=xi,
                pe=pe,
                pa=pa,
                weight=0.25,
                num_sites=sitenum,
                T=200,
                solver=solver,
                max_iter=100,
                final_sample_size=4_000,
                iter_sampling=4000,
                iter_warmup=500,
                show_progress=True,
                seed=1,
                # inits=0.1,
            )

            theta = np.mean(samples["final_sample"][:, :sitenum], axis=0)
            gamma = np.mean(samples["final_sample"][:, sitenum:], axis=0)
            result = shadow_price_opt(
                zbar_1995,
                z_1995,
                forest_area_1995,
                z_2008,
                theta,
                gamma,
                sitenum=sitenum,
                solver=solver,
                timehzn=200,
                pa=pa,
                pe=pe,
            )
            results.append(result)
        results = np.array(results)

    min_index = np.argmin(results)
    min_result = results[min_index]
    min_pe = pe_values[min_index]

    return min_result, min_pe

## det 1043 sites
# min_result,det_1043_pe=shadow_price_cal(sitenum=1043,model='det')
# print("min_result",min_result,"min_pe",det_1043_pe)
# min_result,det_78_pe=shadow_price_cal(sitenum=78,model='det')
# print("min_result",min_result,"min_pe",det_78_pe)


if xi==10000:
    min_result, hmc_1043_pe = shadow_price_cal(sitenum=num_sites, model="det", xi=xi,solver='gurobi',pe_low=seed_i, pe_high=(seed_i+0.01))
    print("min_result", min_result, "min_pe", hmc_1043_pe)
else:
    min_result, hmc_1043_pe = shadow_price_cal(sitenum=num_sites, model="hmc", xi=xi,solver='gurobi',pe_low=seed_i, pe_high=(seed_i+0.01))
    print("min_result", min_result, "min_pe", hmc_1043_pe)



