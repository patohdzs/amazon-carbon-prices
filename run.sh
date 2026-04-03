#!/bin/bash

# Full replication pipeline: ./run.sh -spcdtmha
# Or run individual stages: ./run.sh -s, ./run.sh -p, etc.
#
# Flags:
#   -s  Setup (venv, Python/R packages, CmdStan)
#   -p  Processing (data/raw -> data/clean -> data/processed)
#   -c  Calibration (data/processed -> data/calibration, baseline sampling)
#   -d  Deterministic model (shadow prices, optimization, maps)
#   -t  Time-consistency tau_f checks
#   -h  HMC ambiguity model (adjusted sampling, shadow prices, relative entropy)
#   -m  MPC price-risk model (MC samples, shadow prices, optimization)
#   -a  Analysis (remaining figures, tables, maps)

set -e

usage() {
    echo "Usage: $0 [-spcdtmha]" 1>&2
    exit 1
}

setup_flag='false'
processing_flag='false'
calibration_flag='false'
det_model_flag='false'
time_consistency_flag='false'
mpc_model_flag='false'
hmc_model_flag='false'
analysis_flag='false'

while getopts 'spcdtmha' flag; do
    case "${flag}" in
    s) setup_flag='true' ;;
    p) processing_flag='true' ;;
    c) calibration_flag='true' ;;
    d) det_model_flag='true' ;;
    t) time_consistency_flag='true' ;;
    m) mpc_model_flag='true' ;;
    h) hmc_model_flag='true' ;;
    a) analysis_flag='true' ;;
    *) usage ;;
    esac
done

# ============================================================================
# -s  SETUP
# ============================================================================
if [ "$setup_flag" = "true" ]; then
    # Check prerequisites
    if ! command -v python3 &>/dev/null; then
        echo "Python3 is not installed. Please install Python first."
        exit 1
    fi

    if ! command -v Rscript &>/dev/null; then
        echo "R is not installed. Please install R first."
        exit 1
    fi

    # Create and activate virtual environment
    echo "Creating Python virtual environment..."
    python3 -m venv .venv
    source .venv/bin/activate
    echo "Done!"

    # Install Python dependencies (pinned in pyproject.toml)
    echo "Installing Python dependencies..."
    python -m pip install -e '.[all]'
    echo "Done!"

    # Install CmdStan (pinned version matching cmdstanpy 1.2.0)
    echo "Installing CmdStan..."
    install_cmdstan --version 2.33.1 --overwrite
    echo "Done!"

    # Install R dependencies (pinned in renv.lock)
    echo "Restoring R packages..."
    Rscript -e "renv::restore()"
    echo "Done!"
fi

# Activate venv for all subsequent steps
if [ -d ".venv" ]; then
    source .venv/bin/activate
fi

# ============================================================================
# -p  PROCESSING: data/raw -> data/clean -> data/processed
# ============================================================================
if [ "$processing_flag" = "true" ]; then
    echo "Cleaning raw data..."
    Rscript rsrc/cleaning/_masterfile.R
    echo "Done!"

    echo "Processing clean data..."
    Rscript rsrc/processing/_masterfile.R
    echo "Done!"
fi

# ============================================================================
# -c  CALIBRATION: data/processed -> data/calibration + baseline sampling
# ============================================================================
if [ "$calibration_flag" = "true" ]; then
    echo "Calibrating model parameters..."
    Rscript rsrc/calibration/_masterfile.R
    echo "Done!"

    echo "Running baseline sampling (1043 sites)..."
    python3 pysrc/sampling/baseline.py --sites 1043
    echo "Done!"

    echo "Running baseline sampling (78 sites)..."
    python3 pysrc/sampling/baseline.py --sites 78
    echo "Done!"
fi

# ============================================================================
# -d  DETERMINISTIC MODEL
# ============================================================================
if [ "$det_model_flag" = "true" ]; then
    echo "Running deterministic shadow price calibration..."
    for id in $(seq 60 70); do
        python3 pysrc/bash/shadow_price.py --xi 10000 --sites 1043 --id "$id"
    done
    echo "Done!"

    echo "Running deterministic model..."
    python3 scripts/conduction_det.py
    echo "Done!"

    echo "Generating deterministic maps..."
    Rscript rsrc/analysis/map_1043_det.R
    echo "Done!"
fi

# ============================================================================
# -t  TIME CONSISTENCY TAU_F CHECKS
# ============================================================================
if [ "$time_consistency_flag" = "true" ]; then
    echo "Running time-consistency tau_f checks..."
    python3 pysrc/analysis/time_consistency.py
    echo "Done!"
fi

# ============================================================================
# -h  HMC AMBIGUITY MODEL
# ============================================================================
if [ "$hmc_model_flag" = "true" ]; then
    echo "Running HMC adjusted sampling..."
    for xi in 0.5 1 2 10000; do
        if [[ "$xi" == "0.5" ]]; then
            peearray=(2.9)
        elif [[ "$xi" == "1" ]]; then
            peearray=(6.6 4.7)
        elif [[ "$xi" == "2" ]]; then
            peearray=(5.5)
        else
            peearray=(6.6 4.7)
        fi
        for pee in "${peearray[@]}"; do
            for id in 0 10 15 20 25; do
                python3 pysrc/bash/hmc_sampling.py --id "$id" --xi "$xi" --sites 1043 --pee "$pee"
            done
        done
    done
    echo "Done!"

    echo "Running HMC shadow price calibration..."
    for xi in 0.5 1 2; do
        if [[ "$xi" == "0.5" ]]; then
            idarray=($(seq 20 39))
        elif [[ "$xi" == "1" ]]; then
            idarray=($(seq 40 50))
        elif [[ "$xi" == "2" ]]; then
            idarray=($(seq 50 60))
        fi
        for id in "${idarray[@]}"; do
            python3 pysrc/bash/shadow_price.py --id "$id" --xi "$xi" --sites 1043
        done
    done
    echo "Done!"

    echo "Computing relative entropy..."
    python3 pysrc/bash/relative_entropy.py --xi 1.0 --sites 1043 --pee 4.7
    echo "Done!"

    echo "Running HMC conduction..."
    python3 scripts/conduction_hmc.py
    echo "Done!"
fi

# ============================================================================
# -m  MPC PRICE-RISK MODEL
# ============================================================================
if [ "$mpc_model_flag" = "true" ]; then
    echo "Preparing MPC Monte Carlo samples..."
    for type in baseline constrained shadow_price converge_uncon converge_con; do
        python3 pysrc/mpc/mpc_simulating.py --type "$type"
    done
    echo "Done!"

    echo "Running MPC shadow price search..."
    for xi in 0.5 1 10000; do
        for type in unconstrained constrained; do
            for pe in $(seq 5.0 0.1 6.3); do
                python3 pysrc/mpc/mpc_hmc_sp.py --pe "$pe" --xi "$xi" --type "$type"
            done
        done
    done
    echo "Done!"

    echo "Computing MPC shadow prices..."
    python3 pysrc/mpc/mpc_compute_sp.py
    echo "Done!"

    echo "Running MPC optimization..."
    # Config block 1: det baseline (unconstrained)
    pee=6.3; xi=10000; trig=0; type="unconstrained"
    pearray=($pee $(echo "$pee + 10" | bc) $(echo "$pee + 15" | bc) $(echo "$pee + 20" | bc) $(echo "$pee + 25" | bc))
    for id in $(seq 1 50); do
        for pe in "${pearray[@]}"; do
            python3 pysrc/mpc/mpc_hmc.py --id "$id" --pe "$pe" --xi "$xi" --trig "$trig" --type "$type"
        done
    done

    # Config block 2: hmc xi=1 (unconstrained)
    pee=6.0; xi=1; trig=1; type="unconstrained"
    pearray=($pee $(echo "$pee + 10" | bc) $(echo "$pee + 15" | bc) $(echo "$pee + 20" | bc) $(echo "$pee + 25" | bc))
    for id in $(seq 1 50); do
        for pe in "${pearray[@]}"; do
            python3 pysrc/mpc/mpc_hmc.py --id "$id" --pe "$pe" --xi "$xi" --trig "$trig" --type "$type"
        done
    done

    # Config block 3: hmc xi=0.5 (unconstrained)
    pee=5.6; xi=0.5; trig=1; type="unconstrained"
    pearray=($pee $(echo "$pee + 10" | bc) $(echo "$pee + 15" | bc) $(echo "$pee + 20" | bc) $(echo "$pee + 25" | bc))
    for id in $(seq 1 50); do
        for pe in "${pearray[@]}"; do
            python3 pysrc/mpc/mpc_hmc.py --id "$id" --pe "$pe" --xi "$xi" --trig "$trig" --type "$type"
        done
    done

    # Config block 4: det baseline (constrained)
    pee=6.0; xi=10000; trig=0; type="constrained"
    pearray=($pee $(echo "$pee + 10" | bc) $(echo "$pee + 15" | bc) $(echo "$pee + 20" | bc) $(echo "$pee + 25" | bc))
    for id in $(seq 1 50); do
        for pe in "${pearray[@]}"; do
            python3 pysrc/mpc/mpc_hmc.py --id "$id" --pe "$pe" --xi "$xi" --trig "$trig" --type "$type"
        done
    done

    # Config block 5: hmc xi=1 (constrained)
    pee=5.7; xi=1; trig=1; type="constrained"
    pearray=($pee $(echo "$pee + 10" | bc) $(echo "$pee + 15" | bc) $(echo "$pee + 20" | bc) $(echo "$pee + 25" | bc))
    for id in $(seq 1 50); do
        for pe in "${pearray[@]}"; do
            python3 pysrc/mpc/mpc_hmc.py --id "$id" --pe "$pe" --xi "$xi" --trig "$trig" --type "$type"
        done
    done

    # Config block 6: hmc xi=0.5 (constrained)
    pee=5.2; xi=0.5; trig=1; type="constrained"
    pearray=($pee $(echo "$pee + 10" | bc) $(echo "$pee + 15" | bc) $(echo "$pee + 20" | bc) $(echo "$pee + 25" | bc))
    for id in $(seq 1 50); do
        for pe in "${pearray[@]}"; do
            python3 pysrc/mpc/mpc_hmc.py --id "$id" --pe "$pe" --xi "$xi" --trig "$trig" --type "$type"
        done
    done
    echo "Done!"

    echo "Computing MPC results..."
    python3 pysrc/mpc/mpc_compute.py
    echo "Done!"

    echo "Plotting MPC trajectories..."
    python3 scripts/mpc_trajectory.py
    echo "Done!"
fi

# ============================================================================
# -a  ANALYSIS: remaining figures, tables, maps
# ============================================================================
if [ "$analysis_flag" = "true" ]; then
    echo "Generating calibration maps..."
    Rscript rsrc/analysis/calibration_maps_1043_sites.R
    Rscript rsrc/analysis/calibration_maps_78_sites.R
    echo "Done!"

    echo "Running price estimation..."
    python3 scripts/price_estimation.py
    echo "Done!"

    echo "Computing Bayesian R2..."
    python3 scripts/bayesian_R2.py
    echo "Done!"

    echo "Generating HMC maps..."
    Rscript rsrc/analysis/map_1043_hmc_xi05.R
    Rscript rsrc/analysis/map_1043_hmc_xi1.R
    echo "Done!"
fi
