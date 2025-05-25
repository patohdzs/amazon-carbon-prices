import math
import time
from dataclasses import dataclass
import pandas as pd
import numpy as np
from scipy.special import logsumexp
import pyomo.environ as pyo
from pyomo.environ import value
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
    low_smart_guess=None,
    high_smart_guess=None,
    mode=None,
    type=None,
):
    # pickle_file = 'stan_model/mpc_compiled_model.pkl'

    # if os.path.exists(pickle_file):
    #     # Load the model from the pickle file
    #     sampler = pickle.load(open(pickle_file, 'rb'))
    #     print("Loaded model from pickle.")
    # else:
    #     # Compile the Stan model and save it to the pickle file
    #     sampler = CmdStanModel(
    #         stan_file=get_path("stan_model") / "mpc_adjusted.stan",
    #         cpp_options={"STAN_THREADS": "true"},
    #         force_compile=True,
    #     )
    #     with open(pickle_file, 'wb') as f:
    #         pickle.dump(sampler, f)
    #     print("Compiled model and saved to pickle.")
    
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
    # model.prob = Param(model.J, initialize=0, mutable=True) 
    

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
    
    
    # model.p_ll = Param(initialize=prob_ll, mutable=True)
    # model.p_hh = Param(initialize=prob_hh, mutable=True)
    
    model.p_C0_1 = Param(initialize=prob_ll, mutable=True)
    model.p_C1_1 = Param(initialize=prob_ll, mutable=True)
    model.p_C1_2 = Param(initialize=prob_ll, mutable=True)
    model.p_C2_1 = Param(initialize=prob_ll, mutable=True)
    model.p_C2_2 = Param(initialize=prob_ll, mutable=True)
    model.p_C2_3 = Param(initialize=prob_ll, mutable=True)
    model.p_C2_4 = Param(initialize=prob_ll, mutable=True)
    
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


    if type == "constrained":
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
        if type == "constrained":
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



    results = dict(
        tol=tol,
        T=model.T,
        delta_t=delta,
        alpha=alpha,
        kappa=kappa,
        pf=price_emissions,
        xi=xi,
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
        
        # initialization
        
        
        model.p_C1_1.value = prob_ll
        model.p_C1_2.value = 1-prob_hh
        model.p_C2_1.value = prob_ll
        model.p_C2_2.value = 1-prob_hh
        model.p_C2_3.value = prob_ll
        model.p_C2_4.value = 1-prob_hh
        if model.pa_current.value == 1:
            model.p_C0_1.value = prob_ll
        if model.pa_current.value == 2:
            model.p_C0_1.value = 1-prob_hh
        
        
        model.p_C0_1_ori = model.p_C0_1.value
        model.p_C1_1_ori = model.p_C1_1.value
        model.p_C1_2_ori = model.p_C1_2.value
        model.p_C2_1_ori = model.p_C2_1.value
        model.p_C2_2_ori = model.p_C2_2.value
        model.p_C2_3_ori = model.p_C2_3.value
        model.p_C2_4_ori = model.p_C2_4.value
        
        
        if Q>0 :
            if model.pa_current.value == 1 and low_smart_guess is not None:

                (model.p_C0_1.value,model.p_C1_1.value,model.p_C1_2.value,model.p_C2_1.value,model.p_C2_2.value,model.p_C2_3.value,model.p_C2_4.value) = low_smart_guess

            if model.pa_current.value == 2 and high_smart_guess is not None:

                (model.p_C0_1.value,model.p_C1_1.value,model.p_C1_2.value,model.p_C2_1.value,model.p_C2_2.value,model.p_C2_3.value,model.p_C2_4.value) = high_smart_guess




            
        uncertain_vals_old = np.array([ model.p_C0_1.value, model.p_C1_1.value,model.p_C1_2.value,model.p_C2_1.value,model.p_C2_2.value,model.p_C2_3.value,model.p_C2_4.value]).copy()
        uncertain_vals =  np.array([ model.p_C0_1.value, model.p_C1_1.value,model.p_C1_2.value,model.p_C2_1.value,model.p_C2_2.value,model.p_C2_3.value,model.p_C2_4.value]).copy()

        
        
        while cntr < max_iter and pct_error > tol:            
            print(f"Optimization Iteration[{cntr+1}/{max_iter}]\n")
            
            
            # if cntr>0:
            #     model.p_ll = uncertain_vals[0]
            #     model.p_hh = uncertain_vals[1]
            

            
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
            
            g_0_1, g_1_1, g_1_2, g_2_1, g_2_2, g_2_3, g_2_4=adjust(model,T=time_horizon,S=gamma.size,J=sto_horizon)
            # solved_obj_val = value(model.obj)
            # print(f"Objective value: {solved_obj_val}")
            # print("adjustment value:",object_value_C_t0)
            print("g_0_1:",g_0_1)
            print("g_1_1:",g_1_1)
            print("g_1_2:",g_1_2)
            print("g_2_1:",g_2_1)
            print("g_2_2:",g_2_2)
            print("g_2_3:",g_2_3)
            print("g_2_4:",g_2_4)

            model.p_C0_1.value = min(g_0_1*model.p_C0_1_ori,0.999)
            model.p_C1_1.value = min(g_1_1*model.p_C1_1_ori,0.999)
            model.p_C1_2.value = min(g_1_2*model.p_C1_2_ori,0.999)
            model.p_C2_1.value = min(g_2_1*model.p_C2_1_ori,0.999)
            model.p_C2_2.value = min(g_2_2*model.p_C2_2_ori,0.999)
            model.p_C2_3.value = min(g_2_3*model.p_C2_3_ori,0.999)
            model.p_C2_4.value = min(g_2_4*model.p_C2_4_ori,0.999)

            uncertain_vals = np.array([model.p_C0_1.value, model.p_C1_1.value,model.p_C1_2.value,model.p_C2_1.value,model.p_C2_2.value,model.p_C2_3.value,model.p_C2_4.value]).copy()

            print(f"Parameters from last iteration: {uncertain_vals_old}\n")
            print(
                f"""Parameters from current iteration:
                {uncertain_vals}\n"""
            )


            uncertain_vals_tracker.append(uncertain_vals.copy())
            
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
                }
            )
            

        # if Q+1 in [1,20,40,60,80,100,120,140,160,180,200]:
        #     df = pd.DataFrame({
        #         "p_ll": p_ll_samples.flatten(),  # Convert to 1D array
        #         "p_hh": p_hh_samples.flatten()
        #     })

        #     # Save to CSV
        #     df.to_csv(os.path.join(output_base_path, f"Year_{Q+1}_p_ll_p_hh_samples.csv"), index=False)




        if model.pa_current.value == 2:
            high_smart_guess=uncertain_vals
        if model.pa_current.value == 1:
            low_smart_guess=uncertain_vals



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













def _planner_obj(model):
    
    def compute_flow(model, t, j):
        return (
            -model.pe
            * pyo.quicksum(
                model.kappa * model.z[t + 1, s, j]
                - (model.x[t + 1, s, j] - model.x[t, s, j]) / model.dt
                for s in model.S
            )
            + model.pa[t, j]
            * pyo.quicksum(model.theta[s] * model.z[t + 1, s, j] for s in model.S)
            - (model.zeta_u / 2) * (model.w1[t, j] ** 2)
            - (model.zeta_v / 2) * (model.w2[t, j] ** 2)
        )

    
    object_value_C_t3 = {
        j: (1 - math.exp(-model.delta)) * pyo.quicksum(
            math.exp(-model.delta * (t * model.dt - 4 * model.dt)) *
            compute_flow(model, t, j) * model.dt
            for t in model.T
            if t < max(model.T) and t > 3
        )
        for j in model.J
    }

    
    # Group states into state primes: {1, 2}, {3, 4}, {5, 6}, {7, 8}
    state_C2_C3_mapping = {
        1: [1, 2],
        2: [1, 2],
        3: [3, 4],
        4: [3, 4],
        5: [5, 6],
        6: [5, 6],
        7: [7, 8],
        8: [7, 8],
    }

    state_C2_C3_prob = {
        1: model.p_C2_1, 
        2: 1-model.p_C2_1, 
        3: model.p_C2_2,  
        4: 1-model.p_C2_2,  
        5: model.p_C2_3, 
        6: 1-model.p_C2_3, 
        7: model.p_C2_4,  
        8: 1-model.p_C2_4,  
    }


    object_value_C_t2 = {
        j: (1 - pyo.exp(-model.delta)) *  compute_flow(model, 3, j)
        + pyo.exp(-model.delta) * pyo.quicksum(
            state_C2_C3_prob[mapping] * object_value_C_t3[mapping]
            for mapping in state_C2_C3_mapping[j]
        )
        for j in model.J
    }


    state_C1_C2_prob = {
        1: model.p_C1_1,   
        2: model.p_C1_1,   
        3: 1-model.p_C1_1,  
        4: 1-model.p_C1_1, 
        5: model.p_C1_2,
        6: model.p_C1_2,  
        7: 1-model.p_C1_2,   
        8: 1-model.p_C1_2,
    }

    object_value_C_t1_unique = {}

    object_value_C_t1_unique[0] = (
        (1 - pyo.exp(-model.delta)) *  compute_flow(model, 2, 1) 
        + pyo.exp(-model.delta) *  pyo.quicksum(
            state_C1_C2_prob[j] * object_value_C_t2[j]
            for j in [1, 3]
        )
    )

    object_value_C_t1_unique[1] = (
        (1 - pyo.exp(-model.delta)) * compute_flow(model, 2, 5) 
        + pyo.exp(-model.delta) *  pyo.quicksum(
            state_C1_C2_prob[j] * object_value_C_t2[j]
            for j in [5, 7]
        )
    )


    object_value_C_t0 = (
        (1 - pyo.exp(-model.delta)) *  compute_flow(model, 1, 1) 
        + pyo.exp(-model.delta) * (
            model.p_C0_1 * object_value_C_t1_unique[0]
            + (1-model.p_C0_1) * object_value_C_t1_unique[1]
        )
    )



    return object_value_C_t0




def adjust(model,T,S,J):
    

    dt = model.dt.value
    delta = model.delta.value
    pe = model.pe.value
    kappa = model.kappa.value
    zeta_u = model.zeta_u.value
    zeta_v = model.zeta_v.value
    xi = model.xi.value
    
    z = np.array([[[model.z[t, r, j].value for t in model.T] for r in model.S] for j in model.J])  # Shape: (J, S, T + 1)
    x = np.array([[[model.x[t, r, j].value for t in model.T] for r in model.S] for j in model.J])
    pa = np.array([[model.pa[t, j].value for t in model.T] for j in model.J])
    theta = np.array([model.theta[s] for s in model.S])
    w1 = np.array([[model.w1[t, j].value for t in model.T] for j in model.J])
    w2 = np.array([[model.w2[t, j].value for t in model.T] for j in model.J])


    def compute_flow2(t, j):
        dz = z[j, :, t+1]
        dx = (x[j, :, t+1] - x[j, :, t]) / dt
        flow = (
            -pe * np.sum(kappa * dz - dx)
            + pa[j, t] * np.sum(theta * dz)
            - (zeta_u / 2) * (w1[j, t] ** 2)
            - (zeta_v / 2) * (w2[j, t] ** 2)
        )
        return flow
    
    
    object_value_C_t3 = {}
    for j in range(J):
        object_value_C_t3[j+1] = (1 - np.exp(-delta)) * sum(
            np.exp(-delta * (t * dt - 3 * dt)) * compute_flow2(t, j)
            for t in range(3, T-1)
        ) * dt
        

    g_2_1= np.exp(-1/xi*object_value_C_t3[1])/(model.p_C2_1_ori*np.exp(-1/xi*object_value_C_t3[1])+(1-model.p_C2_1_ori)*np.exp(-1/xi*object_value_C_t3[2]))
    g_2_2= np.exp(-1/xi*object_value_C_t3[3])/(model.p_C2_2_ori*np.exp(-1/xi*object_value_C_t3[3])+(1-model.p_C2_2_ori)*np.exp(-1/xi*object_value_C_t3[4]))
    g_2_3= np.exp(-1/xi*object_value_C_t3[5])/(model.p_C2_3_ori*np.exp(-1/xi*object_value_C_t3[5])+(1-model.p_C2_3_ori)*np.exp(-1/xi*object_value_C_t3[6]))
    g_2_4= np.exp(-1/xi*object_value_C_t3[7])/(model.p_C2_4_ori*np.exp(-1/xi*object_value_C_t3[7])+(1-model.p_C2_4_ori)*np.exp(-1/xi*object_value_C_t3[8]))
    

    # Group states into state primes: {1, 2}, {3, 4}, {5, 6}, {7, 8}
    state_C2_C3_mapping = {
        1: [1, 2],
        2: [1, 2],
        3: [3, 4],
        4: [3, 4],
        5: [5, 6],
        6: [5, 6],
        7: [7, 8],
        8: [7, 8],
    }

    state_C2_C3_prob = {
        1: model.p_C2_1_ori, 
        2: 1-model.p_C2_1_ori,  
        3: model.p_C2_2_ori,  
        4: 1-model.p_C2_2_ori,   
        5: model.p_C2_3_ori, 
        6: 1-model.p_C2_3_ori, 
        7: model.p_C2_4_ori, 
        8: 1-model.p_C2_4_ori,  
    }



    object_value_C_t2 = {
        j+1: (1 - np.exp(-delta)) * compute_flow2(2, j)
        - np.exp(-delta) * xi * np.log(sum(
            state_C2_C3_prob[m] * np.exp(-1/xi*object_value_C_t3[m]) for m in state_C2_C3_mapping[j+1]
        ))
        for j in range(J) 
    }


    g_1_1 = np.exp(-1/xi*object_value_C_t2[1])/(model.p_C1_1_ori*np.exp(-1/xi*object_value_C_t2[1])+(1-model.p_C1_1_ori)*np.exp(-1/xi*object_value_C_t2[3]))
    g_1_2 = np.exp(-1/xi*object_value_C_t2[5])/(model.p_C1_2_ori*np.exp(-1/xi*object_value_C_t2[5])+(1-model.p_C1_2_ori)*np.exp(-1/xi*object_value_C_t2[7]))


    state_C1_C2_prob = {
        1: model.p_C1_1_ori,    
        2: model.p_C1_1_ori,    
        3: 1-model.p_C1_1_ori,   
        4: 1-model.p_C1_1_ori,  
        5: model.p_C1_2_ori, 
        6: model.p_C1_2_ori,   
        7: 1-model.p_C1_2_ori,    
        8: 1-model.p_C1_2_ori, 
    }

    obj_C_t1 = {}
    obj_C_t1[0] = (1 - np.exp(-delta)) * compute_flow2(1, 0) - np.exp(-delta) * xi * np.log( sum(
        state_C1_C2_prob[j] * np.exp(-1/xi*object_value_C_t2[j]) for j in [1,3]
    ))
    obj_C_t1[1] = (1 - np.exp(-delta)) * compute_flow2(1, 4) - np.exp(-delta) * xi * np.log( sum(
        state_C1_C2_prob[j] * np.exp(-1/xi*object_value_C_t2[j]) for j in [5,7]
    ))


    g_0_1 = np.exp(-1/xi*obj_C_t1[0])/(model.p_C0_1_ori*np.exp(-1/xi*obj_C_t1[0])+(1-model.p_C0_1_ori)*np.exp(-1/xi*obj_C_t1[1]))



    return g_0_1, g_1_1, g_1_2, g_2_1, g_2_2, g_2_3, g_2_4




