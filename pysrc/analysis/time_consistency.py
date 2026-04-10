"""
Carrot Policy Functions for Amazon Carbon Prices Project

This module contains functions for computing fund balances and finding optimal
tau values (when fund stops accruing interest) for carrot policy mechanisms.


"""

import warnings
import numpy as np
from typing import Optional, List
from datetime import datetime
from pathlib import Path
import time
import argparse

from pysrc.analysis import value_decomposition
from pysrc.optimization import PlannerSolution, solve_planner_problem
from pysrc.services.data_service import load_productivity_params, load_site_data


def compute_fund_balance(
    X: np.ndarray,
    Z: np.ndarray,
    bf: float,
    tau_f: int,
    delta: float = 0.02,
    kappa: float = 2.094215255
) -> List[float]:
    """
    Compute fund balance over time according to carrot policy dynamics.
    
    Indexing Convention:
    --------------------
    - X, Z are arrays of shape (T+1, n_sites) representing states at times 0, dt, 2dt, ..., T*dt
    - Period t spans from time t to time t+1
    - During period t, flows use END-OF-PERIOD state Z[t+1]
    - X_dot[t] = X[t+1] - X[t] is the carbon change during period t
    
    Fund Dynamics:
    ------------------------------------------
    - Net capture during period t = X_dot[t] - kappa * Z[t+1]
    - Deposit during period t = bf * net_capture[t]
    - Interest ALWAYS accrues at rate delta on fund balance
    - Before tau_f: B[t] = exp(delta) * B[t-1] + deposit[t] (interest retained in fund, continuous compounding)
    - After tau_f:  B[t] = B[t-1] + deposit[t] (interest paid to planner, not added to fund)
    
    The tau_f parameter controls when interest payments to planner begin:
    - Smaller tau_f: Planner receives interest earlier (fund balance grows slower)
    - Larger tau_f: Interest retained in fund longer (fund balance grows larger)
    
    Parameters
    ----------
    X : np.ndarray of shape (T+1, n_sites)
        Carbon stock trajectories
    Z : np.ndarray of shape (T+1, n_sites)
        Agricultural land trajectories
    bf : float
        Fund contribution rate (b_f in the equations)
    tau_f : int
        Period when interest payments to planner begin (0-indexed).
        Before tau_f: interest retained in fund (compounds with balance).
        After tau_f: interest paid to planner (not added to fund balance).
    delta : float, optional
        Discount/interest rate per period (default: 0.02)
    kappa : float, optional
        Agricultural emissions factor in tons CO2 per hectare (default: 2.094215255)
    
    Returns
    -------
    List[float]
        Fund balance at dates t=0,...,T (length T+1).
        balance[0] = B_0 = 0 (initial state, no fund at start)
        balance[t] = B_t, the fund value at date t (end of period t-1, start of period t)
        balance[t+1] = B_{t+1} = exp(delta) * B_t + deposit_t (if t < tau_f)
    
    Raises
    ------
    ValueError
        If tau_f is negative
    
    """
    # Validate inputs
    n_periods = len(X) - 1  # Number of periods (time steps)
    
    if tau_f < 0:
        raise ValueError(f"tau_f must be non-negative, got {tau_f}")
    
    # Note: tau_f > n_periods is allowed - it means interest accrues throughout all periods
    
    # Compute carbon stock changes: X_dot[t] = X[t+1] - X[t] (change during period t)
    X_dot = np.diff(X, axis=0)
    
    # Compute fund deposits for each period t
    # Net capture during period t = X_dot[t] - kappa * Z[t+1]
    # Deposit = bf * net_capture
    fund_deposits = [
        bf * (np.sum(X_dot[t]) - kappa * np.sum(Z[t+1])) 
        for t in range(len(X)-1)
    ]
    
    # Initialize balance: B[0] = 0 (no fund at start)
    # B_t at dates t=0,...,T (length T+1)
    balance = [0.0]  # B[0] = 0
    
    # Before tau_f: interest retained in fund (compounds continuously)
    # B[t+1] = exp(delta) * B[t] + deposits[t]
    # Interest accrues on balance and stays in fund
    growth = float(np.exp(delta))  # Continuous compounding
    for t in range(n_periods):
        if t < tau_f:
            balance.append(growth * balance[t] + fund_deposits[t])
        else:
            balance.append(balance[t] + fund_deposits[t])
    
    return balance


def earliest_tau_search(
    X: np.ndarray,
    Z: np.ndarray,
    V: List[float],
    W: List[float],
    bf: float,
    terminal_rule: str = "lt_zero",
) -> Optional[int]:
    """
    Find the earliest tau (year when fund starts distributing interest) that ensures compliance.
    
    This function searches for the earliest year (tau) that satisfies BOTH:
    1) W[t] - fund_balance[t] < V[t] for ALL time periods t (no defection incentive), AND
    2) Terminal criterion, one of:
       - "lt_zero": W[last] - fund_balance[last] < 0 (strictly negative terminal defection),
       - "lt_v": W[last] - fund_balance[last] < V[last].

    Condition (2) is a stronger terminal requirement than (1) alone when V[last] > 0.
    The search automatically covers the entire horizon from 0 to len(W) - 1.

    The search first validates that compliance is achievable with the given bf by checking
    the case where tau→∞ (fund balance is largest because interest is retained forever in 
    the fund, never paid to planner). If defection occurs even with this maximum fund balance
    to lose, then no finite tau would ensure compliance and an exception is raised.
    
    Parameters
    ----------
    X : np.ndarray
        Carbon stock trajectories
    Z : np.ndarray
        Agricultural land trajectories
    V : List[float]
        Continuation values for each period (length h)
    W : List[float]
        Defection values for each period (length h)
    bf : float
        Fund contribution rate
    terminal_rule : str, optional
        Terminal criterion to apply ("lt_zero" or "lt_v"). Defaults to "lt_zero".
    
    Returns
    -------
    int or None
        The earliest tau value that ensures compliance, or None if no solution found
    
    Raises
    ------
    ValueError
        If even with tau→∞ (interest retained in fund forever, creating maximum fund balance
        and strongest defection deterrent), the planner would still defect at some point.
        This indicates bf is insufficient and must be increased.
    ValueError
        If even with tau→∞, the terminal defection value is non-negative:
        W[last] - fund_balance[last] >= 0. This indicates bf is insufficient to make
        terminal defection strictly negative.
    
    Warnings
    --------
    UserWarning
        If tau=0 already ensures compliance (planner receives interest immediately each period,
        minimal fund balance is sufficient deterrent)
    UserWarning
        If no tau in the entire horizon ensures compliance (shouldn't happen after validation,
        suggests numerical issues)
    
    Notes
    -----
    
    This implementation assumes that the continuation value converges to zero as all forests reach a
    steady state and no further captures remain to generate transfers. If the continuation value for
    some reason converges to a positive constant, e.g. due to an optimal trajectory that involves 
    maintaining a positive level of agriculture (Z), this implementation will overestimate tau_f.
    
    The tau parameter controls when interest payments to planner begin:
    - tau = 0: Interest paid to planner immediately each period (fund balance = deposits only)
    - tau = t: Interest retained in fund until year t, then paid to planner from year t onward
    - tau → ∞: Interest always retained in fund (maximum fund balance, planner never receives interest)
    
    Higher tau means:
    - Planner waits longer to receive interest payments
    - Fund balance grows larger (more compounding before payouts begin)
    - Stronger deterrent against defection (larger fund to lose)
    
    Note: Uses continuous compounding (exp(delta)) before tau_f, matching the continuous-time
    formulation of the model.
    """

    h = len(W)  # Horizon length
    last_year = h - 1

    if terminal_rule not in {"lt_zero", "lt_v"}:
        raise ValueError(
            f"Unknown terminal_rule='{terminal_rule}'. Use 'lt_zero' or 'lt_v'."
        )

    def terminal_condition(terminal_value: float) -> bool:
        if terminal_rule == "lt_zero":
            return terminal_value < 0
        return terminal_value < V[last_year]

    # First check if compliance is achievable with the given bf.
    # With tau→∞, the fund balance is maximized because interest is retained in the fund forever
    # (never paid to planner), compounding indefinitely. This creates the largest possible
    # fund balance to lose upon defection - the strongest deterrent.
    # If planner still defects in this case, no finite tau will help.
    fund_balance_max = compute_fund_balance(X, Z, bf, int(1e6))  # tau→∞ (effectively infinite)
    defection_value_max = [W[t] - fund_balance_max[t] for t in range(h)]
    
    # Check if any defection value exceeds continuation value
    if any(d > c for d, c in zip(defection_value_max, V)):
        first_defection = next(t for t, (d, c) in enumerate(zip(defection_value_max, V)) if d > c)
        raise ValueError(
            f'Even with tau→∞ (interest retained in fund forever, maximum deterrent), '
            f'planner would defect at year {first_defection}. '
            f'Current bf={bf} is insufficient. Increase bf or b to ensure compliance.'
        )

    # Check terminal requirement under maximum deterrent
    if not terminal_condition(defection_value_max[last_year]):
        if terminal_rule == "lt_zero":
            raise ValueError(
                f'Even with tau→∞ (maximum fund deterrent), terminal defection is not negative: '
                f'W[last]-fund_balance[last]={defection_value_max[last_year]:.4f} >= 0. '
                f'Current bf={bf} is insufficient to make terminal defection strictly negative.'
            )
        raise ValueError(
            f'Even with tau→∞ (maximum fund deterrent), terminal condition fails: '
            f'W[last]-fund_balance[last]={defection_value_max[last_year]:.4f} >= V[last]={V[last_year]:.4f}. '
            f'Current bf={bf} is insufficient for W[last]-fund_balance[last] < V[last].'
        )

    # Conducting search    
    for tau in range(h):
        fund_balance = compute_fund_balance(X, Z, bf, tau)
        
        # Check compliance at ALL time periods: W[t] - fund_balance[t] < V[t] for all t
        defection_values = [W[t] - fund_balance[t] for t in range(h)]
        compliance_check = all(d < v for d, v in zip(defection_values, V))
        terminal_negative_check = terminal_condition(defection_values[last_year])
        
        if compliance_check and terminal_negative_check:
            if tau == 0:
                warnings.warn(
                    'tau=0 already ensures compliance. '
                    'Planner receives interest payments immediately (minimal fund balance needed).',
                    UserWarning
                )
            return tau
    
    # Handle case where no solution is found (shouldn't happen after validation)
    warnings.warn(
        f'No tau in range [0, {h-1}] ensures compliance. '
        f'This is unexpected after validation - possible numerical issues.',
        UserWarning
    )
    return None


def tau0_feasibility_check(
    X: np.ndarray,
    Z: np.ndarray,
    V: List[float],
    W: List[float],
    bf: float,
    check_all_periods: bool,
    terminal_rule: str,
    terminal_year: int,
    delta: float = 0.02,
    kappa: float = 2.094215255,
) -> bool:
    """
    Check whether tau_f=0 is feasible for a given bf.

    Feasibility conditions:
    - If check_all_periods=True: W[t] - fund[t] < V[t] for all t.
    - Always enforce terminal condition according to terminal_rule at terminal_year.
    """
    if terminal_rule not in {"lt_zero", "lt_v"}:
        raise ValueError(
            f"Unknown terminal_rule='{terminal_rule}'. Use 'lt_zero' or 'lt_v'."
        )

    h = len(W)
    if terminal_year < 0 or terminal_year >= h:
        raise ValueError(f"terminal_year={terminal_year} is outside [0, {h-1}]")

    fund_balance = compute_fund_balance(X, Z, bf, tau_f=0, delta=delta, kappa=kappa)

    if check_all_periods:
        if any(np.isnan(W[t]) for t in range(h)):
            raise ValueError("W contains NaN values, cannot run all-period feasibility check.")
        for t in range(h):
            if W[t] - fund_balance[t] >= V[t]:
                return False

    terminal_value = W[terminal_year] - fund_balance[terminal_year]
    if terminal_rule == "lt_zero":
        return terminal_value < 0
    return terminal_value < V[terminal_year]


def smallest_bf_search_tau0(
    X: np.ndarray,
    Z: np.ndarray,
    V: List[float],
    W: List[float],
    check_all_periods: bool,
    terminal_rule: str,
    terminal_year: int,
    delta: float = 0.02,
    kappa: float = 2.094215255,
    tol: float = 1e-8,
) -> Optional[float]:
    """
    Find the smallest non-negative bf that makes tau_f=0 feasible.

    This search uses the linear structure of tau_f=0:
      fund_balance_t(bf) = bf * fund_balance_t(1),
    so constraints in bf are linear inequalities.

    Returns
    -------
    float or None
        Smallest bf (up to numerical tolerance) if feasible, else None.
    """
    if terminal_rule not in {"lt_zero", "lt_v"}:
        raise ValueError(
            f"Unknown terminal_rule='{terminal_rule}'. Use 'lt_zero' or 'lt_v'."
        )

    h = len(W)
    if terminal_year < 0 or terminal_year >= h:
        raise ValueError(f"terminal_year={terminal_year} is outside [0, {h-1}]")

    if check_all_periods and any(np.isnan(W[t]) for t in range(h)):
        raise ValueError("W contains NaN values, cannot run all-period bf search.")
    if np.isnan(W[terminal_year]):
        raise ValueError("W[terminal_year] is NaN, cannot run terminal bf search.")

    # At tau=0, fund_balance is linear in bf, so evaluate once at bf=1.
    fund_unit = compute_fund_balance(X, Z, bf=1.0, tau_f=0, delta=delta, kappa=kappa)

    lower_bound = -np.inf
    upper_bound = np.inf

    def apply_constraint(t: int, threshold: float) -> bool:
        """
        Apply strict inequality:
            W[t] - bf * fund_unit[t] < threshold
        which is equivalent to:
            bf * fund_unit[t] > W[t] - threshold
        """
        nonlocal lower_bound, upper_bound

        s = fund_unit[t]
        rhs = W[t] - threshold
        eps = 1e-14

        if abs(s) <= eps:
            # No bf leverage at this t: must already satisfy strict inequality.
            return 0.0 > rhs

        bound = rhs / s
        if s > 0:
            lower_bound = max(lower_bound, bound)
        else:
            upper_bound = min(upper_bound, bound)
        return True

    if check_all_periods:
        for t in range(h):
            if not apply_constraint(t=t, threshold=V[t]):
                return None

    terminal_threshold = 0.0 if terminal_rule == "lt_zero" else V[terminal_year]
    if not apply_constraint(t=terminal_year, threshold=terminal_threshold):
        return None

    # Enforce non-negative bf.
    lower_bound = max(lower_bound, 0.0)

    # Because inequalities are strict, the feasible interval is (lower_bound, upper_bound).
    if upper_bound <= lower_bound + tol:
        return None

    # Prefer bf=0 exactly when feasible and admissible.
    if lower_bound == 0.0 and tau0_feasibility_check(
        X=X,
        Z=Z,
        V=V,
        W=W,
        bf=0.0,
        check_all_periods=check_all_periods,
        terminal_rule=terminal_rule,
        terminal_year=terminal_year,
        delta=delta,
        kappa=kappa,
    ):
        return 0.0

    candidate = lower_bound + tol
    if candidate >= upper_bound:
        return None

    # Guard against floating-point edge cases near strict boundaries.
    for _ in range(10):
        if tau0_feasibility_check(
            X=X,
            Z=Z,
            V=V,
            W=W,
            bf=candidate,
            check_all_periods=check_all_periods,
            terminal_rule=terminal_rule,
            terminal_year=terminal_year,
            delta=delta,
            kappa=kappa,
        ):
            return candidate
        candidate += tol
        if candidate >= upper_bound:
            return None

    return None


def save_curves_to_h5(V, W, new_W, fund_value, results, metadata, output_dir=None):
    import h5py

    if output_dir is None:
        project_root = Path(__file__).parent.parent.parent
        output_dir = project_root / "output" / "data"
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    b_str = str(metadata["b"]).replace(".", "p")
    bf_str = str(metadata.get("bf", 1)).replace(".", "p")
    filename = output_dir / f"curves_b{b_str}_bf{bf_str}_{timestamp}.h5"

    with h5py.File(filename, "w") as file:
        file.create_dataset("continuation_value", data=np.array(V))
        file.create_dataset("defecting_value", data=np.array(W))
        file.create_dataset("defecting_minus_fund", data=np.array(new_W))
        file.create_dataset("fund_value", data=np.array(fund_value))

        solution_group = file.create_group("planner_solution")
        solution_group.create_dataset("Z", data=results.Z, compression="gzip")
        solution_group.create_dataset("X", data=results.X, compression="gzip")
        solution_group.create_dataset("U", data=results.U, compression="gzip")
        solution_group.create_dataset("V", data=results.V, compression="gzip")

        for key, value in metadata.items():
            if isinstance(value, (int, float, str, bool)):
                file.attrs[key] = value
            elif isinstance(value, np.ndarray):
                file.create_dataset(f"metadata_{key}", data=value)
            else:
                file.attrs[key] = str(value)

        file.attrs["timestamp"] = timestamp
        file.attrs["created_at"] = datetime.now().isoformat()

    print(f"Saved curves and full solution to: {filename}")
    return filename


def _run_workflow(
    b: float,
    bf: float,
    pee: float,
    pa: float,
    num_sites: int,
    T: int,
    h: int,
    kappa: float,
    delta: float,
    solver: str,
    save_to_h5: bool,
    generate_plots: bool,
    output_dir: str = None,
    check_all_periods: bool = True,
):
    """
    Main workflow function.

    Parameters
    ----------
    b : float
        Base transfer payment ($/ton CO2)
    bf : float
        Fund contribution ($/ton CO2)
    pee : float
        Emissions price ($/ton CO2)
    pa : float
        Cattle price ($/hectare)
    num_sites : int
        Number of sites (78 or 1043)
    T : int
        Full optimization horizon
    h : int
        Analysis horizon (years)
    kappa : float
        Emissions factor (tons CO2 per hectare)
    delta : float
        Discount rate
    solver : str
        Solver name (e.g., 'gurobi')
    save_to_h5 : bool
        Set to True to save results to H5 file
    generate_plots : bool
        Set to True to generate and display plots
    output_dir : str, optional
        Output directory for H5 file (None = default)
    check_all_periods : bool, optional
        If True, enforce both no-defection-at-all-periods and terminal negativity.
        If False, enforce only terminal negativity at the final period.
    """
    import matplotlib.pyplot as plt

    start_time = time.time()
    print("=" * 80)
    print("FUND WORKFLOW ANALYSIS")
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)
    print("\nScenario Parameters:")
    print(f"  b (base transfer): {b} $/ton CO2")
    print(f"  bf (fund contribution): {bf} $/ton CO2")
    print(f"  pee (emissions price): {pee} $/ton CO2")
    print(f"  pa (cattle price): {pa} $/hectare")
    print(f"  Total price (pee + b + bf): {pee + b + bf} $/ton CO2")
    print("\nModel Parameters:")
    print(f"  Number of sites: {num_sites}")
    print(f"  Optimization horizon (T): {T}")
    print(f"  Analysis horizon (h): {h}")
    print(f"  kappa (emissions factor): {kappa} tons CO2/hectare")
    print(f"  delta (discount rate): {delta}")
    print(f"  Solver: {solver}")
    print("=" * 80)

    print("\n[Step 1] Loading data...")
    step_start = time.time()
    print(f"  Loading site data for {num_sites} sites...")
    zbar_2017, z_2017, forest_area_2017 = load_site_data(num_sites)
    print("  ✓ Site data loaded: zbar, z_2017, forest_area_2017")

    print("  Loading productivity parameters...")
    theta, gamma = load_productivity_params(num_sites)
    print("  ✓ Productivity parameters loaded: theta, gamma")

    print("  Computing initial carbon stocks (x0 = gamma * forest_area_2017)...")
    x0_vals = gamma * forest_area_2017

    elapsed = time.time() - step_start
    print(f"\n✓ Step 1 complete in {elapsed:.2f} seconds")
    print(f"  Loaded {num_sites} sites")
    print(f"  Total zbar: {np.sum(zbar_2017):.6f} billion hectares")
    print(f"  Total z_2017: {np.sum(z_2017):.6f} billion hectares")
    print(f"  Total forest area: {np.sum(forest_area_2017):.6f} billion hectares")
    print(f"  Total initial carbon stock (x0): {np.sum(x0_vals):.6f} billion tons")
    print(f"  Theta range: [{np.min(theta):.4f}, {np.max(theta):.4f}]")
    print(f"  Gamma range: [{np.min(gamma):.4f}, {np.max(gamma):.4f}]")

    print("\n[Step 2] Solving planner problem...")
    step_start = time.time()
    print("  Parameters:")
    print(f"    Time horizon: {T + h} periods (T={T} + h={h})")
    print(f"    Number of sites: {num_sites}")
    print(f"    Price emissions: ${pee + b + bf:.2f}/ton CO2 (pee={pee} + b={b} + bf={bf})")
    print(f"    Price cattle: ${pa:.2f}/hectare")
    print(f"    Solver: {solver}")
    print("    Note: Using T+h periods to ensure enough data for continuation values")
    print("  Solving optimization problem (this may take several minutes)...")

    results = solve_planner_problem(
        time_horizon=T + h,
        theta=theta,
        gamma=gamma,
        x0=x0_vals,
        zbar=zbar_2017,
        z0=z_2017,
        price_emissions=pee + b + bf,
        price_cattle=pa,
        solver=solver,
    )

    elapsed = time.time() - step_start
    print(f"\n✓ Step 2 complete in {elapsed:.2f} seconds ({elapsed / 60:.2f} minutes)")
    print("  Solution arrays shape:")
    print(f"    Z (agricultural area): {results.Z.shape}")
    print(f"    X (carbon stock): {results.X.shape}")
    print(f"    U (deforestation): {results.U.shape}")
    print(f"    V (value): {results.V.shape}")
    print("  Initial values:")
    print(f"    Z[0] sum: {np.sum(results.Z[0]):.6f} billion hectares")
    print(f"    X[0] sum: {np.sum(results.X[0]):.6f} billion tons")
    print("  Terminal values:")
    print(f"    Z[{T + h - 1}] sum: {np.sum(results.Z[T + h - 1]):.6f} billion hectares")
    print(f"    X[{T + h - 1}] sum: {np.sum(results.X[T + h - 1]):.6f} billion tons")

    print(f"\n[Step 3] Computing continuation values V(t) for {h} periods...")
    step_start = time.time()
    print(f"  Computing continuation value from each starting period t in [0, {h - 1}]...")
    print(f"  Transfer payment: b + bf = ${b + bf:.2f}/ton CO2")
    print("  Progress: ", end="", flush=True)

    V = []
    for t in range(h):
        if (t + 1) % 20 == 0 or t == 0 or t == h - 1:
            print(f"t={t}... ", end="", flush=True)

        value = value_decomposition(
            pee=pee,
            pa=pa,
            b=b + bf,
            theta=theta,
            solution=PlannerSolution(
                results.Z[t:],
                results.X[t:],
                results.U[t:],
                results.V[t:],
            ),
            T=T,
        )["total_PV"]

        V.append(value)

    elapsed = time.time() - step_start
    print(f"\n✓ Step 3 complete in {elapsed:.2f} seconds ({elapsed / 60:.2f} minutes)")
    print(f"  Computed {len(V)} continuation values")
    print(f"  V[0] = ${V[0]:.2f} billion (initial period)")
    print(f"  V[{h // 4}] = ${V[h // 4]:.2f} billion (25% through horizon)")
    print(f"  V[{h // 2}] = ${V[h // 2]:.2f} billion (50% through horizon)")
    print(f"  V[{3 * h // 4}] = ${V[3 * h // 4]:.2f} billion (75% through horizon)")
    print(f"  V[{h - 1}] = ${V[h - 1]:.2f} billion (final period)")

    terminal_year = h - 1 if check_all_periods else 100
    if terminal_year >= h:
        raise ValueError(
            f"terminal_year={terminal_year} is outside h={h}. "
            "Increase h or adjust terminal-year rule."
        )

    print(f"\n[Step 4] Computing defection values W(t) for {h} periods...")
    step_start = time.time()
    if check_all_periods:
        print(f"⚠️  WARNING: This step can take many hours (solving {h} optimization problems)")
        print(f"  Each optimization problem solves for {T} periods with {num_sites} sites")
        print("  Estimated time per problem: 1-5 minutes (varies by solver and problem size)")
        print(f"  Total estimated time: {h * 2:.0f}-{h * 5:.0f} minutes")
        print("\n  Starting computation...")
        print("  Progress: ", end="", flush=True)

        W = []
        for t in range(h):
            show_progress = (t + 1) % 10 == 0 or t == 0 or t == h - 1
            if show_progress:
                elapsed_so_far = time.time() - step_start
                avg_time_per_problem = elapsed_so_far / (t + 1) if t > 0 else 0.0
                remaining_problems = h - (t + 1)
                eta_seconds = avg_time_per_problem * remaining_problems if avg_time_per_problem > 0 else 0.0
                eta_minutes = eta_seconds / 60
                print(f"\n    t={t} ({t + 1}/{h}) | Elapsed: {elapsed_so_far / 60:.1f} min", end="")
                if eta_minutes > 0:
                    print(f" | ETA: {eta_minutes:.1f} min", end="")
                print(" | Solving defection problem... ", end="", flush=True)
            elif (t + 1) % 5 == 0:
                print(".", end="", flush=True)

            iter_start = time.time()

            defection_results = solve_planner_problem(
                time_horizon=T,
                theta=theta,
                gamma=gamma,
                x0=results.X[t],
                zbar=zbar_2017,
                z0=results.Z[t],
                price_emissions=pee,
                price_cattle=pa,
                solver=solver,
            )

            if show_progress:
                iter_elapsed = time.time() - iter_start
                print(f"({iter_elapsed:.1f}s) | Computing value... ", end="", flush=True)

            defection_value = value_decomposition(
                pee=pee,
                pa=pa,
                b=0,
                theta=theta,
                solution=PlannerSolution(
                    defection_results.Z,
                    defection_results.X,
                    defection_results.U,
                    defection_results.V,
                ),
                T=T,
            )["total_PV"]

            W.append(defection_value)

            if show_progress:
                print(f"W[{t}] = ${defection_value:.2f}B")

        elapsed = time.time() - step_start
        print(f"\n✓ Step 4 complete in {elapsed:.2f} seconds ({elapsed / 60:.2f} minutes)")
        print(f"  Computed {len(W)} defection values")
        print(f"  W[0] = ${W[0]:.2f} billion (initial period)")
        print(f"  W[{h // 4}] = ${W[h // 4]:.2f} billion (25% through horizon)")
        print(f"  W[{h // 2}] = ${W[h // 2]:.2f} billion (50% through horizon)")
        print(f"  W[{3 * h // 4}] = ${W[3 * h // 4]:.2f} billion (75% through horizon)")
        print(f"  W[{h - 1}] = ${W[h - 1]:.2f} billion (final period)")
        print(f"  Average time per defection problem: {elapsed / h:.2f} seconds")
    else:
        print("  Terminal-only mode enabled for `check_all_periods=False`.")
        print(f"  Solving only one defection problem at t={terminal_year} (using X[{terminal_year}], Z[{terminal_year}]).")

        defection_results = solve_planner_problem(
            time_horizon=T,
            theta=theta,
            gamma=gamma,
            x0=results.X[terminal_year],
            zbar=zbar_2017,
            z0=results.Z[terminal_year],
            price_emissions=pee,
            price_cattle=pa,
            solver=solver,
        )

        defection_value = value_decomposition(
            pee=pee,
            pa=pa,
            b=0,
            theta=theta,
            solution=PlannerSolution(
                defection_results.Z,
                defection_results.X,
                defection_results.U,
                defection_results.V,
            ),
            T=T,
        )["total_PV"]

        W = [np.nan] * h
        W[terminal_year] = defection_value

        elapsed = time.time() - step_start
        print(f"\n✓ Step 4 complete in {elapsed:.2f} seconds ({elapsed / 60:.2f} minutes)")
        print(f"  Computed terminal-only defection value W[{terminal_year}] = ${W[terminal_year]:.2f} billion")

    print("\n[Step 5] Finding optimal tau...")
    step_start = time.time()
    bf_candidates = [0.15 * float(b), 0.20 * float(b)]
    bf_candidate_label = ", ".join(f"{candidate:.6f}" for candidate in bf_candidates)
    tau_results = []
    if check_all_periods:
        print(f"  Running earliest_tau_search for bf in {{{bf_candidate_label}}}:")
        print("    1) W[t] - fund_balance[t] < V[t] for all t")
        print(f"    2a) W({terminal_year}) - fund_balance({terminal_year}) < 0")
        print(f"    2b) W({terminal_year}) - fund_balance({terminal_year}) < V({terminal_year})")
        for bf_candidate in bf_candidates:
            tau_zero = None
            tau_lt_v = None
            try:
                tau_zero = earliest_tau_search(
                    X=results.X,
                    Z=results.Z,
                    V=V,
                    W=W,
                    bf=bf_candidate,
                    terminal_rule="lt_zero",
                )
                if tau_zero is None:
                    tau_zero = h - 1
            except ValueError as error:
                print(f"  ✗ bf={bf_candidate:.2f} terminal<0 -> {error}")

            try:
                tau_lt_v = earliest_tau_search(
                    X=results.X,
                    Z=results.Z,
                    V=V,
                    W=W,
                    bf=bf_candidate,
                    terminal_rule="lt_v",
                )
                if tau_lt_v is None:
                    tau_lt_v = h - 1
            except ValueError as error:
                print(f"  ✗ bf={bf_candidate:.2f} terminal<V[last] -> {error}")

            tau_results.append(
                {
                    "bf": bf_candidate,
                    "tau_terminal_lt_zero": tau_zero,
                    "tau_terminal_lt_v": tau_lt_v,
                }
            )
            tau_zero_text = "NA" if tau_zero is None else str(tau_zero)
            tau_lt_v_text = "NA" if tau_lt_v is None else str(tau_lt_v)
            print(
                f"  ✓ bf={bf_candidate:.2f} -> "
                f"tau_f[terminal<0]={tau_zero_text}, "
                f"tau_f[terminal<V[last]]={tau_lt_v_text}"
            )
    else:
        print(f"  Running terminal-only tau search for bf in {{{bf_candidate_label}}}:")
        print(f"    2a) W({terminal_year}) - fund_balance({terminal_year}) < 0")
        print(f"    2b) W({terminal_year}) - fund_balance({terminal_year}) < V({terminal_year})")
        last_year = terminal_year

        def terminal_only_tau(bf_candidate: float, use_v_threshold: bool) -> Optional[int]:
            max_fund_balance = compute_fund_balance(
                X=results.X,
                Z=results.Z,
                bf=bf_candidate,
                tau_f=int(1e6),
                delta=delta,
                kappa=kappa,
            )
            threshold = V[last_year] if use_v_threshold else 0.0
            max_terminal_value = W[last_year] - max_fund_balance[last_year]
            if max_terminal_value >= threshold:
                return None

            tau_candidate = None
            for tau in range(h):
                fund_balance_tau = compute_fund_balance(
                    X=results.X,
                    Z=results.Z,
                    bf=bf_candidate,
                    tau_f=tau,
                    delta=delta,
                    kappa=kappa,
                )
                if W[last_year] - fund_balance_tau[last_year] < threshold:
                    tau_candidate = tau
                    break

            if tau_candidate is None:
                tau_candidate = h - 1
            return tau_candidate

        for bf_candidate in bf_candidates:
            tau_zero = terminal_only_tau(
                bf_candidate=bf_candidate,
                use_v_threshold=False,
            )
            tau_lt_v = terminal_only_tau(
                bf_candidate=bf_candidate,
                use_v_threshold=True,
            )
            tau_results.append(
                {
                    "bf": bf_candidate,
                    "tau_terminal_lt_zero": tau_zero,
                    "tau_terminal_lt_v": tau_lt_v,
                }
            )
            tau_zero_text = "NA" if tau_zero is None else str(tau_zero)
            tau_lt_v_text = "NA" if tau_lt_v is None else str(tau_lt_v)
            print(
                f"  ✓ bf={bf_candidate:.2f} -> "
                f"tau_f[terminal<0]={tau_zero_text}, "
                f"tau_f[terminal<V[last]]={tau_lt_v_text}"
            )

    print("\n  Running smallest bf search for fixed tau_f=0:")
    if check_all_periods:
        print("    1) W[t] - fund_balance[t] < V[t] for all t")
    else:
        print(f"    1) Terminal-only mode at t={terminal_year} (no all-period checks)")
    print(f"    2a) W({terminal_year}) - fund_balance({terminal_year}) < 0")
    print(f"    2b) W({terminal_year}) - fund_balance({terminal_year}) < V({terminal_year})")

    bf_tau0_lt_zero = smallest_bf_search_tau0(
        X=results.X,
        Z=results.Z,
        V=V,
        W=W,
        check_all_periods=check_all_periods,
        terminal_rule="lt_zero",
        terminal_year=terminal_year,
        delta=delta,
        kappa=kappa,
    )
    bf_tau0_lt_v = smallest_bf_search_tau0(
        X=results.X,
        Z=results.Z,
        V=V,
        W=W,
        check_all_periods=check_all_periods,
        terminal_rule="lt_v",
        terminal_year=terminal_year,
        delta=delta,
        kappa=kappa,
    )

    bf_tau0_lt_zero_text = "NA" if bf_tau0_lt_zero is None else f"{bf_tau0_lt_zero:.8f}"
    bf_tau0_lt_v_text = "NA" if bf_tau0_lt_v is None else f"{bf_tau0_lt_v:.8f}"
    print(
        "  ✓ smallest bf with tau_f=0 -> "
        f"bf_min[terminal<0]={bf_tau0_lt_zero_text}, "
        f"bf_min[terminal<V[last]]={bf_tau0_lt_v_text}"
    )

    elapsed = time.time() - step_start
    print(f"\n✓ Step 5 complete in {elapsed:.2f} seconds")

    selected_result = None
    for result in tau_results:
        if abs(float(result["bf"]) - float(bf)) <= 1e-12:
            selected_result = result
            break
    if selected_result is None:
        if not tau_results:
            print("✗ No tau result produced for any bf candidate.")
            return
        selected_result = tau_results[0]
        print(f"  Input bf={bf:.2f} was not solved; selecting bf={selected_result['bf']:.2f} for Step 6.")

    tau_optimal = selected_result["tau_terminal_lt_zero"]
    if tau_optimal is None:
        tau_optimal = selected_result["tau_terminal_lt_v"]
        if tau_optimal is None:
            print(f"✗ bf={selected_result['bf']:.2f} has no feasible tau for either terminal criterion.")
            return
        print(
            f"  bf={selected_result['bf']:.2f} has no tau for terminal<0; "
            f"using tau_f={tau_optimal} from terminal<V[last] for Step 6."
        )

    selected_bf = float(selected_result["bf"])

    print(f"\n[Step 6] Computing fund balance with optimal tau={tau_optimal}...")
    step_start = time.time()
    print("  Parameters:")
    print(f"    Fund contribution rate (bf): ${selected_bf:.2f}/ton CO2")
    print(f"    Interest start year (tau): {tau_optimal}")
    print(f"    Discount rate (delta): {delta}")
    print(f"    Emissions factor (kappa): {kappa} tons CO2/hectare")
    print("  Computing fund balance trajectory...")

    fund_balance_optimal = compute_fund_balance(
        X=results.X,
        Z=results.Z,
        bf=selected_bf,
        tau_f=tau_optimal,
        delta=delta,
        kappa=kappa,
    )

    print("  Computing defection value minus fund balance (W - fund)...")
    new_W = [W[t] - fund_balance_optimal[t] for t in range(h)]

    elapsed = time.time() - step_start
    print(f"\n✓ Step 6 complete in {elapsed:.2f} seconds")
    print(f"  Fund balance trajectory computed for {len(fund_balance_optimal)} periods")
    print("\n  Key statistics:")
    print(f"    Fund balance at t=0: ${fund_balance_optimal[0]:.2f} billion")
    print(f"    Fund balance at t={tau_optimal}: ${fund_balance_optimal[tau_optimal]:.2f} billion")
    print(f"    Fund balance at t={h - 1}: ${fund_balance_optimal[h - 1]:.2f} billion")
    print(
        "    Maximum fund balance: "
        f"${np.max(fund_balance_optimal):.2f} billion (at t={np.argmax(fund_balance_optimal)})"
    )
    print(f"\n  Compliance check at terminal period (t={terminal_year}):")
    print(f"    W (defection value): ${W[terminal_year]:.2f} billion")
    print(f"    Fund balance: ${fund_balance_optimal[terminal_year]:.2f} billion")
    print(f"    W - Fund: ${new_W[terminal_year]:.2f} billion")
    print(
        "    Condition 2 (W-Fund < 0): "
        f"{new_W[terminal_year] < 0} ({'✓ PASS' if new_W[terminal_year] < 0 else '✗ FAIL'})"
    )
    terminal_negative = new_W[terminal_year] < 0
    print(
        f"  Terminal criterion result: {terminal_negative} "
        f"({'✓ PASS' if terminal_negative else '✗ FAIL'})"
    )

    if generate_plots:
        print("\n[Step 7] Generating plots...")
        step_start = time.time()
        print("  Creating value paths plot...")
        plt.figure(figsize=(12, 6))
        plt.plot(V, label="Continuation Value (V)", linestyle="-", linewidth=2.5)
        plt.plot(W, label="Defection Value (W)", linestyle="-", linewidth=2)
        plt.plot(
            new_W,
            label=f"Defection Value - Fund (tau={tau_optimal})",
            linestyle="--",
            linewidth=2,
        )
        plt.axhline(y=0, color="black", linestyle=":", linewidth=1, alpha=0.5)
        plt.axvline(
            x=tau_optimal,
            color="red",
            linestyle=":",
            linewidth=1,
            alpha=0.5,
            label=f"tau={tau_optimal}",
        )
        plt.xlabel("Time (years)")
        plt.ylabel("$ billion")
        plt.title(f"Value Paths for b={b}, bf={selected_bf} (Optimal tau={tau_optimal})")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plot_dir = Path("output") / "figures"
        plot_dir.mkdir(parents=True, exist_ok=True)
        plot_filename = plot_dir / f"curves_b{b}_bf{bf}_tau{tau_optimal}.png"
        plt.savefig(plot_filename, dpi=150, bbox_inches="tight")
        elapsed = time.time() - step_start
        print(f"✓ Step 7 complete in {elapsed:.2f} seconds")
        print(f"  ✓ Plot saved to: {plot_filename}")
        plt.close()

    if save_to_h5:
        print("\n[Step 8] Saving results to H5 file...")
        step_start = time.time()
        print("  Preparing metadata...")
        metadata = {
            "n_sites": num_sites,
            "pee": pee,
            "pa": pa,
            "b": b,
            "bf_input": bf,
            "bf": selected_bf,
            "tau_f": tau_optimal,
            "T": T,
            "h": h,
            "kappa": kappa,
            "delta": delta,
            "solver": solver,
            "initial_zbar_sum": float(np.sum(zbar_2017)),
            "initial_z_sum": float(np.sum(z_2017)),
            "initial_forest_area_sum": float(np.sum(forest_area_2017)),
            "initial_x0_sum": float(np.sum(x0_vals)),
        }
        print(f"  Metadata keys: {list(metadata.keys())}")
        print("  Saving curves (V, W, new_W, fund_value) and full solution to H5...")
        output_file = save_curves_to_h5(
            V=V,
            W=W,
            new_W=new_W,
            fund_value=fund_balance_optimal,
            results=results,
            metadata=metadata,
            output_dir=output_dir,
        )
        elapsed = time.time() - step_start
        print(f"\n✓ Step 8 complete in {elapsed:.2f} seconds")
        print(f"  ✓ Saved results to: {output_file}")
        print("  File contains:")
        print(f"    - Continuation values (V): {len(V)} values")
        print(f"    - Defection values (W): {len(W)} values")
        print(f"    - Adjusted defection values (new_W): {len(new_W)} values")
        print(f"    - Fund balance trajectory: {len(fund_balance_optimal)} values")
        print(f"    - Full planner solution (Z, X, U, V) for {T} periods")
        print("    - Metadata dictionary with all parameters")

    total_elapsed = time.time() - start_time
    print("\n" + "=" * 80)
    print("WORKFLOW COMPLETE")
    print(f"Total execution time: {total_elapsed:.2f} seconds ({total_elapsed / 60:.2f} minutes)")
    print(f"Completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)


def main(
    b: float,
    bf: float,
    pee: float,
    pa: float,
    num_sites: int,
    T: int,
    h: int,
    kappa: float,
    delta: float,
    solver: str,
    save_to_h5: bool,
    generate_plots: bool,
    output_dir: str = None,
    check_all_periods: bool = True,
):
    return _run_workflow(
        b=b,
        bf=bf,
        pee=pee,
        pa=pa,
        num_sites=num_sites,
        T=T,
        h=h,
        kappa=kappa,
        delta=delta,
        solver=solver,
        save_to_h5=save_to_h5,
        generate_plots=generate_plots,
        output_dir=output_dir,
        check_all_periods=check_all_periods,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run time-consistency workflow."
    )
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument(
        "--check-all-periods",
        dest="check_all_periods",
        action="store_true",
        help="Solve defection values for all periods and enforce no-defection at all t (default).",
    )
    mode_group.add_argument(
        "--terminal-only",
        dest="check_all_periods",
        action="store_false",
        help="Solve only terminal defection value and enforce terminal criterion only.",
    )
    parser.set_defaults(check_all_periods=True)
    cli_args = parser.parse_args()

    print("Starting time-consistency workflow...\n", flush=True)
    if cli_args.check_all_periods:
        print("Mode: check_all_periods=True (full no-defection path checks)\n", flush=True)
    else:
        print("Mode: check_all_periods=False (terminal-only checks)\n", flush=True)

    try:
        main(
            b=25,
            bf=3.75,
            pee=6.8,
            pa=41.11,
            num_sites=1043,
            T=200,
            h=101,
            kappa=2.094215255,
            delta=0.02,
            solver="gurobi",
            save_to_h5=False,
            generate_plots=False,
            output_dir=None,
            check_all_periods=cli_args.check_all_periods,
        )
    except KeyboardInterrupt:
        print("\n\n" + "=" * 80)
        print("WORKFLOW INTERRUPTED BY USER")
        print("=" * 80)
    except Exception as error:
        print("\n\n" + "=" * 80)
        print("ERROR: Workflow failed")
        print("=" * 80)
        print(f"Error type: {type(error).__name__}")
        print(f"Error message: {str(error)}")
        import traceback

        print("\nTraceback:")
        traceback.print_exc()
        raise
