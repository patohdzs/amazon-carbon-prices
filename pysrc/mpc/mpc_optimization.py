import math
import time
from dataclasses import dataclass
import pandas as pd
import numpy as np
import pyomo.environ as pyo
import dill as pickle 
from cmdstanpy import CmdStanModel
import os
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    NonNegativeReals,
    Objective,
    Param,
    RangeSet,
    Var,
    maximize,
    Set,
)
from pyomo.opt import SolverFactory
from ..services.file_service import get_path

@dataclass
class PlannerSolution:
    Z: np.ndarray
    X: np.ndarray
    U: np.ndarray
    V: np.ndarray


def mpc_solve_planner_problem(
    x0,
    z0,
    zbar,
    gamma,
    theta,
    dt=1,
    time_horizon=200,
    price_emissions=20.76,
    price_cattle=44.75,
    price_low = 35.71,
    price_high = 44.25,
    alpha=0.045007414,
    delta=0.02,
    kappa=2.094215255,
    zeta_u=1.66e-4 * 1e9,
    zeta_v=1.00e-4 * 1e9,
    solver="gurobi",
    sto_horizon=8,
    xi=10000,
    prob_ll=0.706,
    prob_hh=0.829,
    id=1,
    max_iter=20000,
    tol=0.005,
    final_sample_size=5_000,
    forest_area_2017=None,
    weight=0.5,
    low_smart_guess_p_ll=None,
    mode=None,
):
    pickle_file = 'stan_model/mpc_compiled_model.pkl'

    if os.path.exists(pickle_file):
        # Load the model from the pickle file
        sampler = pickle.load(open(pickle_file, 'rb'))
        print("Loaded model from pickle.")
    else:
        # Compile the Stan model and save it to the pickle file
        sampler = CmdStanModel(
            stan_file=get_path("stan_model") / "mpc_adjusted.stan",
            cpp_options={"STAN_THREADS": "true"},
            force_compile=True,
        )
        with open(pickle_file, 'wb') as f:
            pickle.dump(sampler, f)
        print("Compiled model and saved to pickle.")
    
    output_base_path = os.path.join(
        str(get_path("output")),
        "sampling",
        solver,
        f"78sites",
        f"mpc",
        f"xi_{xi}",
        f"id_{id}",
        f"pe_{price_emissions}",
    )
    if not os.path.exists(output_base_path):
        os.makedirs(output_base_path)
    
    model = ConcreteModel()

    # Indexing sets for time and sites
    model.T = RangeSet(time_horizon + 1)
    model.S = RangeSet(gamma.size)
    model.J = RangeSet(sto_horizon)

    # Parameters
    model.x0 = Param(model.S, initialize=_np_to_dict(x0))
    model.z0 = Param(model.S, initialize=_np_to_dict(z0))
    model.zbar = Param(model.S, initialize=_np_to_dict(zbar))
    model.gamma = Param(model.S, initialize=_np_to_dict(gamma))
    model.theta = Param(model.S, initialize=_np_to_dict(theta))
    model.delta = Param(initialize=delta)
    model.pe = Param(initialize=price_emissions)


    model.pa_current = Param(initialize=1, mutable=True)
    def initialize_pa(model,  t, j):
        if t == 1:
            return price_low if model.pa_current.value == 1 else price_high
        elif t == 2:
            if j in [1, 2, 3, 4]:
                return price_low
            else:
                return price_high
        elif t == 3:
            if j in [1, 2] or j in [5, 6]:
                return price_low
            else:
                return price_high
        elif t == 4:
            if j in [1, 3, 5, 7]:
                return price_low
            else:
                return price_high
        else:  
            if j in [1, 3, 5, 7]:
                return price_low
            else:
                return price_high
    model.pa = Param(model.T, model.J, mutable=True, initialize=initialize_pa)    
    model.prob = Param(model.J, initialize=0, mutable=True) 
    

    dep = [
        (1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 5, 1), (1, 6, 1), (1, 7, 1), (1, 8, 1),
        (2, 2, 1), (2, 3, 1), (2, 4, 1), (2, 6, 2), (2, 7, 2), (2, 8, 2),
        (3, 2, 1), (3, 4, 3), (3, 6, 5), (3, 8, 7),
    ]

    model.dep = Set(initialize=dep, dimen=3)  
    
    
    
    model.zeta_u = Param(initialize=zeta_u)
    model.zeta_v = Param(initialize=zeta_v)

    model.alpha = Param(initialize=alpha)
    model.kappa = Param(initialize=kappa)
    model.xi = Param(initialize=xi)
    model.dt = Param(initialize=dt)
    
    
    model.p_ll = Param(initialize=prob_ll, mutable=True)
    model.p_hh = Param(initialize=prob_hh, mutable=True)

    # Variables
    model.x = Var(model.T, model.S, model.J)
    model.z = Var(model.T, model.S, model.J, within=NonNegativeReals)
    model.u = Var(model.T, model.S, model.J, within=NonNegativeReals)
    model.v = Var(model.T, model.S, model.J, within=NonNegativeReals)

    # Auxilary variables
    model.w1 = Var(model.T,model.J)
    model.w2 = Var(model.T,model.J)

    # Constraints
    model.zdot_def = Constraint(model.T, model.S,model.J, rule=_zdot_const)
    model.xdot_def = Constraint(model.T, model.S,model.J, rule=_xdot_const)
    model.w1_def = Constraint(model.T, model.J,rule=_w1_const)
    model.w2_def = Constraint(model.T, model.J,rule=_w2_const)


    model.non_x_def = Constraint(
        model.dep,
        model.S,
        rule=non_x_def_rule
    )

    model.non_z_def = Constraint(
        model.dep,
        model.S,
        rule=non_z_def_rule
    )

    model.non_u_def = Constraint(
        model.dep,
        model.S,
        rule=non_u_def_rule
    )

    model.non_v_def = Constraint(
        model.dep,
        model.S,
        rule=non_v_def_rule
    )

    model.non_w_def = Constraint(
        model.dep,
        rule=non_w_def_rule
    )



    # Define the objective
    model.obj = Objective(rule=_planner_obj, sense=maximize)

    Z, X, U, V = [], [], [], []
    Z.append(np.array([model.z0[r] for r in model.S]))
    X.append(np.array([model.x0[r] for r in model.S]))


    if price_low == 32.44:
        file_path = (
            get_path("output", "simulation", "mpc_path","baseline","constrained") / f"mc_{id}.csv"
        ) 
        print("constrained price used")
    else:
        file_path = (
            get_path("output", "simulation", "mpc_path","baseline","unconstrained") / f"mc_{id}.csv"
        ) 
        print("unconstrained price used")


    if mode =="converge":
        if price_low == 32.44:
            file_path = (
                get_path("output", "simulation", "mpc_path","constrained",f"xi_{xi}",f"pe_{price_emissions}") / f"mc_{id}.csv"
            ) 
        else:
            file_path = (
                get_path("output", "simulation", "mpc_path","unconstrained",f"xi_{xi}",f"pe_{price_emissions}") / f"mc_{id}.csv"
            ) 


    if file_path is None:
        raise ValueError(f"Invalid mode '{mode}' or price_low '{price_low}'. Could not determine file_path.")

    
    pa_list = np.array(pd.read_csv(file_path))[:,1]


    collected_ensembles = {}
    uncertain_vals_old = np.array([prob_ll, prob_hh]).copy()
    uncertain_vals = np.array([prob_ll, prob_hh]).copy()
    results = dict(
        tol=tol,
        T=model.T,
        delta_t=delta,
        alpha=alpha,
        kappa=kappa,
        pf=price_emissions,
        xi=xi,
        final_sample_size=final_sample_size,
    )

    iteration_period=time_horizon
    if mode=="sp":
        iteration_period=21
    for Q in range(iteration_period):
    # for Q in range(1):
        
        model.pa_current.set_value(pa_list[Q])

        for t in model.T:
            for j in model.J:
                model.pa[t, j] = initialize_pa(model, t, j)

        print(f"Time {Q+1}: {model.pa_current.value}")

        
        for j in model.J:
            for s in model.S:
                model.z[:, s,j].setub(model.zbar[s])
                model.u[max(model.T), s,j].fix(0)
                model.v[max(model.T), s,j].fix(0)
                model.w1[max(model.T),j].fix(0)
                model.w2[max(model.T),j].fix(0)
                
                if Q ==0:
                    model.x[min(model.T), s,j].fix(model.x0[s])
                    model.z[min(model.T), s,j].fix(model.z0[s])
                if Q >0:
                    model.x[min(model.T), s,j].fix(model.x[2, s,1])
                    model.z[min(model.T), s,j].fix(model.z[2, s,1])
                    
        # Initialize error & iteration counter
        pct_error = np.infty
        cntr = 0

        uncertain_vals_tracker = []
        pct_error_tracker = []
        
        ## initialization
        if Q>0 :
            if model.pa_current.value == 1:
                if low_smart_guess_p_ll == None:
                    model.p_ll = prob_ll
                    model.p_hh = prob_hh
                else:
                    model.p_ll = low_smart_guess_p_ll      
                    model.p_hh = low_smart_guess_p_hh  
            if model.pa_current.value == 2:
                model.p_ll = high_smart_guess_p_ll      
                model.p_hh = high_smart_guess_p_hh  
            uncertain_vals_old = np.array([model.p_ll.value, model.p_hh.value]).copy()
        
        while cntr < max_iter and pct_error > tol:            
            print(f"Optimization Iteration[{cntr+1}/{max_iter}]\n")
            
            
            if cntr>0:
                model.p_ll = uncertain_vals[0]
                model.p_hh = uncertain_vals[1]
            
            # Calculate the probabilities based on `pa_current`
            if model.pa_current.value == 1:
                prob = {
                    1: model.p_ll**3,
                    2: model.p_ll**2 * (1 - model.p_ll),
                    3: model.p_ll * (1 - model.p_ll) * (1 - model.p_hh),
                    4: model.p_ll * (1 - model.p_ll) * model.p_hh,
                    5: (1 - model.p_ll) * (1 - model.p_hh) * model.p_ll,
                    6: (1 - model.p_ll) * (1 - model.p_hh) * (1 - model.p_ll),
                    7: (1 - model.p_ll) * model.p_hh * (1 - model.p_hh),
                    8: (1 - model.p_ll) * model.p_hh * model.p_hh,
                }
            elif model.pa_current.value == 2:
                prob = {
                    1: (1 - model.p_hh) * model.p_ll**2,
                    2: (1 - model.p_hh) * model.p_ll * (1 - model.p_ll),
                    3: (1 - model.p_hh) * (1 - model.p_ll) * (1 - model.p_hh),
                    4: (1 - model.p_hh) * (1 - model.p_ll) * model.p_hh,
                    5: model.p_hh * (1 - model.p_hh) * model.p_ll,
                    6: model.p_hh * (1 - model.p_hh) * (1 - model.p_ll),
                    7: model.p_hh * model.p_hh * (1 - model.p_hh),
                    8: model.p_hh * model.p_hh * model.p_hh,
                }

            # print("prob_Test", prob)
            # print(f"Sum of probabilities = {sum(prob.values())}")

            for j in model.J:
                model.prob[j] = prob[j] 
                        
            

            # Solve the model
            opt = SolverFactory(solver)
            print("Solving the optimization problem...")
            start_time = time.time()
            if solver == "gams":
                opt.solve(model, tee=True, solver="cplex", mtype="qcp")
            else:
                opt.options['Threads'] = 8 
                opt.solve(model, tee=True)

            print(f"Done! Time elapsed: {time.time()-start_time} seconds.")
            
            Z_temp = np.array([[[model.z[t, r, j].value for t in model.T] for r in model.S] for j in model.J])  # Shape: (J, S, T + 1)
            X_temp = np.array([[[model.x[t, r, j].value for t in model.T] for r in model.S] for j in model.J])
            U_temp = np.array([[[model.u[t, r, j].value for t in model.T] for r in model.S] for j in model.J])[:, :, :-1]  # Shape: (J, S, T)
            V_temp = np.array([[[model.v[t, r, j].value for t in model.T] for r in model.S] for j in model.J])[:, :, :-1]  # Shape: (J, S, T)


            # HMC sampling
            print("Starting HMC sampling...\n")
            model_data = dict(
                T=time_horizon,
                S=len(list(model.S)),
                J=len(list(model.J)),
                alpha=alpha,
                zbar_2017=zbar,
                zeta_u=zeta_u,
                zeta_v=zeta_v,
                xi=xi,
                kappa=kappa,
                pa=np.array([[model.pa[t, j].value for j in model.J] for t in model.T]),
                pe=price_emissions,
                pa_current=model.pa_current.value,
                Z=Z_temp,
                U=U_temp,
                V=V_temp,
                X=X_temp,
                gamma=gamma,
                theta=theta,
                **_dynamics_matrices(alpha, delta),
            )

            # Sampling from adjusted distribution
            sampling_time = time.time()
            fit = sampler.sample(
                data=model_data,
                iter_sampling=1000,
                iter_warmup=500,
                show_progress=True,
                seed=123,
                # show_console=True,
            )
            sampling_time = time.time() - sampling_time
            print(f"Finished sampling! Elapsed Time: {sampling_time} seconds\n")
            print(fit.diagnose())

            p_ll_samples = fit.stan_variable("p_ll")[:,0].reshape(-1, 1)
            p_hh_samples = fit.stan_variable("p_hh")[:,0].reshape(-1, 1)


            uncertainty_adj_samples = np.concatenate(
            (p_ll_samples, p_hh_samples), axis=1
            )
            
            collected_ensembles.update({cntr: uncertainty_adj_samples.copy()})

            print(f"Parameters from last iteration: {uncertain_vals_old}\n")
            print(
                f"""Parameters from current iteration:
                {np.mean(uncertainty_adj_samples, axis=0)}\n"""
            )

            # Compute exponentially-smoothened new params
            uncertain_vals = (
                weight * np.mean(uncertainty_adj_samples, axis=0)
                + (1 - weight) * uncertain_vals_old
            )

            uncertain_vals_tracker.append(uncertain_vals.copy())
            print(f"Updated uncertain values: {uncertain_vals}\n")

            pct_error = np.max(
                np.abs(uncertain_vals_old - uncertain_vals) / uncertain_vals_old
            )

            pct_error_tracker.append(pct_error)

            print(
                f"""
                Percentage Error = {pct_error}
                """
            )

            # Exchange parameter values
            uncertain_vals_old = uncertain_vals

            # Increase the counter
            cntr += 1

            # Update results directory
            results.update(
                {   "time_period":Q+1,
                    "cntr": cntr,
                    "pct_error_tracker": np.asarray(pct_error_tracker),
                    "uncertain_vals_tracker": np.asarray(uncertain_vals_tracker),
                    "collected_ensembles": collected_ensembles,
                }
            )
            
            if pct_error <= tol:
                final_samples = uncertain_vals
                results.update({"time_period": Q+1,
                                "final_sample": final_samples})

        # if Q+1 in [1,20,40,60,80,100,120,140,160,180,200]:
        #     df = pd.DataFrame({
        #         "p_ll": p_ll_samples.flatten(),  # Convert to 1D array
        #         "p_hh": p_hh_samples.flatten()
        #     })

        #     # Save to CSV
        #     df.to_csv(os.path.join(output_base_path, f"Year_{Q+1}_p_ll_p_hh_samples.csv"), index=False)




        if model.pa_current.value == 2:
            high_smart_guess_p_ll=uncertain_vals[0]
            high_smart_guess_p_hh=uncertain_vals[1]
        if model.pa_current.value == 1:
            low_smart_guess_p_ll=uncertain_vals[0]
            low_smart_guess_p_hh=uncertain_vals[1]


        Z.append(np.array([model.z[2, r,1].value  for r in model.S]))
        X.append(np.array([model.x[2, r,1].value  for r in model.S]))
        U.append(np.array([model.u[1, r,1].value  for r in model.S]))
        V.append(np.array([model.v[1, r,1].value  for r in model.S]))

        print("year done:",Q+1)
        

    Z = np.array(Z)
    X = np.array(X)
    U = np.array(U)
    V = np.array(V)
    

    
        
    return PlannerSolution(Z, X, U, V)


def vectorize_trajectories(traj: PlannerSolution):
    Z = traj.Z.T
    U = traj.U[:-1, :].T
    V = traj.V[:-1, :].T
    return {
        "Z": Z,
        "U": U,
        "V": V,
    }


def _planner_obj(model):
    
    

    object_value = pyo.quicksum( model.prob[j]*
            pyo.quicksum(
                math.exp(-model.delta * (t * model.dt - model.dt))
                * (
                    -model.pe
                    * pyo.quicksum(
                        model.kappa * model.z[t + 1, s,j]
                        - (model.x[t + 1, s,j] - model.x[t, s,j]) / model.dt
                        for s in model.S
                    )
                    + model.pa[t,j]
                    * pyo.quicksum(model.theta[s] * model.z[t + 1, s,j] for s in model.S)
                    - (model.zeta_u / 2) * (model.w1[t,j] ** 2)
                    - (model.zeta_v / 2) * (model.w2[t,j] ** 2)
                )
                * model.dt
                for t in model.T
                if t < max(model.T)
            )
        for j in model.J       
    )
    
    
    return object_value



def _zdot_const(model, t, s,j):
    if t < max(model.T):
        return (model.z[t + 1, s,j] - model.z[t, s,j]) / model.dt == (
            model.u[t, s,j] - model.v[t, s,j]
        )
    else:
        return Constraint.Skip


def _xdot_const(model, t, s,j):
    if t < max(model.T):
        return (model.x[t + 1, s,j] - model.x[t, s,j]) / model.dt == (
            -model.gamma[s] * model.u[t, s,j]
            - model.alpha * model.x[t, s,j]
            + model.alpha * model.gamma[s] * (model.zbar[s] - model.z[t, s,j])
        )
    else:
        return Constraint.Skip


def _w1_const(model, t,j):
    if t < max(model.T):
        return model.w1[t,j] == pyo.quicksum(model.u[t, s,j] for s in model.S)
    else:
        return Constraint.Skip


def _w2_const(model, t,j):
    if t < max(model.T):
        return model.w2[t,j] == pyo.quicksum(model.v[t, s,j] for s in model.S)
    else:
        return Constraint.Skip


def _np_to_dict(x):
    return dict(enumerate(x.flatten(), 1))




def non_x_def_rule(model, t, j, j1, s):
    return model.x[t+1, s, j] == model.x[t+1, s, j1]

def non_z_def_rule(model, t, j, j1, s):
    return model.z[t+1,s, j] == model.z[ t+1,s, j1]



def non_u_def_rule(model, t, j, j1, s):
    return model.u[t,s, j] == model.u[ t,s, j1]


def non_v_def_rule(model, t, j, j1, s):
    return model.v[t,s, j] == model.v[ t,s, j1]

def non_w_def_rule(model, t, j, j1):
    return model.w1[t,j] == model.w1[ t,j1]

def _dynamics_matrices( alpha, delta, T=200, dt=1):
    # Create dynamics matrices
    arr = np.cumsum(
        np.triu(np.ones((T, T))),
        axis=1,
    ).T
    Bdym = (1 - alpha) ** (arr - 1)
    Bdym[Bdym > 1] = 0.0
    Adym = np.arange(1, T + 1)
    alpha_p_Adym = np.power(1 - alpha, Adym)

    # Other placeholders!
    ds_vect = np.exp(-delta * np.arange(T) * dt)
    ds_vect = np.reshape(ds_vect, (ds_vect.size, 1)).flatten()
    return {"alpha_p_Adym": alpha_p_Adym, "Bdym": Bdym, "ds_vect": ds_vect}
