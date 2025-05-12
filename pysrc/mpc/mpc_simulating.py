import os
import shutil

import numpy as np
import pandas as pd

from pysrc.sampling import mpc_estimation
from pysrc.services.file_service import get_path
import argparse
parser = argparse.ArgumentParser(description="mpc simulation")
parser.add_argument("--type",type=str,default="baseline")
args = parser.parse_args()
type = args.type






def mc_samples(location,p_low,p_high,price_low=35.76,price_high=44.32):
    np.random.seed(123)


    num_simulations = 200

    states = {"low": price_low, "high":  price_high }
    initial_state = "high"

    probability_matrix = {
        "low": {"low": p_low, "high": 1-p_low},
        "high": {"low": 1-p_high, "high": p_high},
    }

    # Number of observations
    num_observations = 200

    for i in range(1, num_simulations + 1):
        # Generating the Markov chain for each simulation
        current_state = initial_state
        markov_chain = [states[current_state]]

        for _ in range(num_observations - 1):
            next_state = np.random.choice(
                list(probability_matrix[current_state].keys()),
                p=list(probability_matrix[current_state].values()),
            )
            markov_chain.append(states[next_state])
            current_state = next_state

        transformed_markov_chain = [
            1 if price == 35.76 else 2 for price in markov_chain
        ]

        # Creating an index from 1 to 200
        index = range(1, num_observations + 1)

        # Combine index and transformed Markov chain into a DataFrame
        markov_chain_df = pd.DataFrame(
            {"Index": index, "scenario": transformed_markov_chain}
        )

        # Specify the filename (e.g., mc_1.csv,..., mc_100.csv)
        csv_filename = f"/mc_{i}.csv"
        output = location
        # Output to CSV file
        markov_chain_df.to_csv(output + csv_filename, index=False)

    return "mc sampling is done"






if type =="baseline":
    
    location = os.path.join(
        str(get_path("output")),
        "simulation",
        "mpc_path",
        "baseline",
        "unconstrained",
    )

    # Create the directory if it doesn't exist
    os.makedirs(location, exist_ok=True)

    # Run the simulation
    mc_samples(location=location, p_low=0.707, p_high=0.826)

    print("Mpc simulation done for baseline")




if type =="constrained":
    
    location = os.path.join(
        str(get_path("output")),
        "simulation",
        "mpc_path",
        "baseline",
        "constrained",
    )

    # Create the directory if it doesn't exist
    os.makedirs(location, exist_ok=True)

    # Run the simulation
    mc_samples(location=location, p_low=0.762, p_high=0.959,price_low=32.49,price_high=42.85)

    print("Mpc simulation done for constrained,baseline")




if type == "shadow_price":
    location1 = os.path.join(
        str(get_path("output")),
        "simulation",
        "mpc_path",
        "baseline",
        "unconstrained",
    )
    location2 = os.path.join(
        str(get_path("output")),
        "simulation",
        "mpc_path",
        "baseline",
        "constrained",
    )


    num_observations = 200
    index = range(1, num_observations + 1)
    
    predefined_states = [
        2, 2, 2, 2, 2, 2, 2, 2, 1, 1,
        1, 1, 1, 2, 2, 2, 2, 1, 1, 2
    ]

    transformed_markov_chain = predefined_states + [1] * (num_observations - 20)
    
    markov_chain_df = pd.DataFrame(
        {"Index": index, "scenario": transformed_markov_chain}
    )

    markov_chain_df.to_csv(location1 +  f"/mc_999.csv", index=False)
    markov_chain_df.to_csv(location2 +  f"/mc_999.csv", index=False)



    transformed_markov_chain = [1] * num_observations  # All 200 years set to state 1
    markov_chain_df = pd.DataFrame( {"Index": index, "scenario": transformed_markov_chain})
    markov_chain_df.to_csv(location1 + f"/mc_997.csv", index=False)
    markov_chain_df.to_csv(location2 + f"/mc_997.csv", index=False)
    
    transformed_markov_chain = [2] * num_observations  # All 200 years set to state 1
    markov_chain_df = pd.DataFrame( {"Index": index, "scenario": transformed_markov_chain})
    markov_chain_df.to_csv(location1 + f"/mc_998.csv", index=False)
    markov_chain_df.to_csv(location2 + f"/mc_998.csv", index=False)
    
    print("MPC simulation done for shadow_price")










if type =="converge_uncon":


    xi_values = [0.2,1.0]
    b_values = [0, 10, 15, 20, 25]


    prob_dict = {
        1: {
            "b0": {"prob_ll": 0.96, "prob_hh": 0.45},
            "b10": {"prob_ll": 0.82, "prob_hh": 0.73},
            "b15": {"prob_ll": 0.76, "prob_hh": 0.79},
            "b20": {"prob_ll": 0.74, "prob_hh": 0.81},
            "b25": {"prob_ll": 0.74, "prob_hh": 0.81},
        },
        0.2: {
            "b0": {"prob_ll": 0.99, "prob_hh": 0.25},
            "b10": {"prob_ll": 0.95, "prob_hh": 0.47},
            "b15": {"prob_ll": 0.91, "prob_hh": 0.59},
            "b20": {"prob_ll": 0.85, "prob_hh": 0.68},
            "b25": {"prob_ll": 0.83, "prob_hh": 0.72},
        }
    }


    for xi in xi_values:
        for b in b_values:
            if xi==1.0:
                pee=5.7
            elif xi==0.2:
                pee=5.5
            pe = pee +b
            # Get the probability values for the current b value
            p_low = prob_dict[xi][f"b{b}"]["prob_ll"]
            p_high = prob_dict[xi][f"b{b}"]["prob_hh"]

            print("xi",xi,"b",b,"p_low",p_low,"p_high",p_high)

            # Define the output location
            location = os.path.join(
                str(get_path("output")),
                "simulation",
                "mpc_path",
                "unconstrained",
                f"xi_{xi}",
                f"pe_{pe}",
            )

            # Create the directory if it doesn't exist
            os.makedirs(location, exist_ok=True)

            # Run the simulation
            mc_samples(location=location, p_low=p_low, p_high=p_high)

    print("Monte Carlo simulations completed for all xi and b values. unconstrained")
    
    
    
    
    
    
    
if type =="converge_con":


    xi_values = [0.2,1.0]
    b_values = [0, 10, 15, 20, 25]


    prob_dict = {
        1: {
            "b0": {"prob_ll": 0.97, "prob_hh": 0.41},
            "b10": {"prob_ll": 0.83, "prob_hh": 0.72},
            "b15": {"prob_ll": 0.77, "prob_hh": 0.78},
            "b20": {"prob_ll": 0.75, "prob_hh": 0.80},
            "b25": {"prob_ll": 0.74, "prob_hh": 0.81},
        },
        0.2: {
            "b0": {"prob_ll": 0.99, "prob_hh": 0.18},
            "b10": {"prob_ll": 0.96, "prob_hh": 0.45},
            "b15": {"prob_ll": 0.91, "prob_hh": 0.57},
            "b20": {"prob_ll": 0.87, "prob_hh": 0.66},
            "b25": {"prob_ll": 0.85, "prob_hh": 0.69},
        }
    }


    for xi in xi_values:
        for b in b_values:
            if xi==1:
                pee=5.1
            elif xi==0.2:
                pee=5.0
            pe = pee +b
            # Get the probability values for the current b value
            p_low = prob_dict[xi][f"b{b}"]["prob_ll"]
            p_high = prob_dict[xi][f"b{b}"]["prob_hh"]

            print("xi",xi,"b",b,"p_low",p_low,"p_high",p_high)

            # Define the output location
            location = os.path.join(
                str(get_path("output")),
                "simulation",
                "mpc_path",
                "constrained",
                f"xi_{xi}",
                f"pe_{pe}",
            )

            # Create the directory if it doesn't exist
            os.makedirs(location, exist_ok=True)

            # Run the simulation
            mc_samples(location=location, p_low=p_low, p_high=p_high)

    print("Monte Carlo simulations completed for all xi and b values. constrained")
