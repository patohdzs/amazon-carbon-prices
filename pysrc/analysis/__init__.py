import numpy as np

from pysrc.optimization import PlannerSolution


def compute_planner_value(
    pee: float,
    pa: float,
    b: float,
    r: float,
    theta: np.ndarray,
    zbar: np.ndarray,
    solution: PlannerSolution,
    T: int = 200,
    delta: float = 0.02,
    kappa: float = 2.094215255,
    zeta_u: float = 1.66e-4 * 1e9,
    zeta_v: float = 1.00e-4 * 1e9,
):
    # Compute change in X
    X_dot = np.diff(solution.X, axis=0)

    # Compute agricultural output
    agr_output = [
        pa * np.dot(solution.Z[t + 1], theta) / ((1 + delta) ** t) for t in range(T)
    ]
    agr_output = np.sum(agr_output)

    # Compute net transfers
    net_transfers = [
        -b * (kappa * np.sum(solution.Z[t + 1]) - np.sum(X_dot[t])) / ((1 + delta) ** t)
        for t in range(T)
    ]
    net_transfers = np.sum(net_transfers)

    # Compute forest services
    forest_services = [
        -pee
        * (kappa * np.sum(solution.Z[t + 1]) - np.sum(X_dot[t]))
        / ((1 + delta) ** t)
        for t in range(T)
    ]
    forest_services = np.sum(forest_services)

    # Compute adjustment costs
    adj_costs = [
        (
            (zeta_u / 2) * (np.sum(solution.U[t])) ** 2
            + (zeta_v / 2) * (np.sum(solution.V[t])) ** 2
        )
        / ((1 + delta) ** t)
        for t in range(T)
    ]
    adj_costs = np.sum(adj_costs)

    tfff = [
        r * max(np.sum(zbar - 101 * solution.Z[t + 1]), 0) / ((1 + delta) ** t)
        for t in range(T)
    ]

    tfff = np.sum(tfff)

    # Compute total net present value
    planner_value = agr_output + net_transfers + forest_services - adj_costs + tfff

    return {
        "pe": pee + b,
        "b": b,
        "agr_output": agr_output,
        "net_transfers": net_transfers,
        "forest_services": forest_services,
        "adj_costs": adj_costs,
        "tfff": tfff,
        "planner_value": planner_value,
    }


def compute_transfers(
    res,
    res_base,
    pee: float,
    b,
    num_years,
    delta=0.02,
    kappa=2.094215255,
):
    # Compute change in X
    X_dot = np.diff(res.X, axis=0)
    X_dot_base = np.diff(res_base.X, axis=0)

    # Compute net captured emissions for base case
    net_emissions_base = [
        -kappa * res_base.Z[t + 1] + X_dot_base[t] for t in range(num_years)
    ]
    net_emissions_base = np.sum(net_emissions_base)

    # Compute net captured emissions
    net_emissions = [-kappa * res.Z[t + 1] + X_dot[t] for t in range(num_years)]
    net_emissions = np.sum(net_emissions)

    net_transfers = [
        -b * (kappa * res.Z[t + 1] - X_dot[t]) / ((1 + delta) ** t)
        for t in range(num_years)
    ]
    net_transfers = np.sum(net_transfers)

    effective_cost = net_transfers / (net_emissions - net_emissions_base)

    return {
        "pe": pee + b,
        "b": b,
        "net captured emissions": net_emissions,
        "discounted net transfers": net_transfers,
        "discounted effective costs": effective_cost,
    }
