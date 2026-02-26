from pathlib import Path
from scipy.linalg import expm
from scipy.linalg import fractional_matrix_power

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import click

@click.command()
@click.argument('popsize', type=int)
@click.argument('selcoef', type=float)
@click.argument('mutrate', type=float)

@click.option('--input', '-i', default='results/escsim', help='Input folder for simulation results (default: results/escsim)')
@click.option('--output', '-o', default='results/coalescent_densities', help='Output folder for coalescent densities (default: results/coalescent_densities)')

def write_file(popsize, selcoef, mutrate, input, output):
    """
    Reads the results from the forward-simulation,
    calculates the coalescent densities and effective population size,
    and writes the results to an output file.
    """
    file_path = Path(input)
    if not file_path.exists():
        click.echo(f"Error: File {file_path} does not exist.")
        return
    
    results_file = pd.read_csv(file_path, sep="\t")
    mut_burden_profile = np.array(results_file["profile"].iloc[0].split(","), dtype=float)
    v = read_params(results_file)[3]
    max_t = 5000
    delta_t = 0.001
    t = np.arange(0, max_t, delta_t)
    coalescent_densities = calc_effective_coalescent_rates(mut_burden_profile, popsize, mutrate, selcoef, v, delta_t=delta_t, max_t=max_t)
    effective_pop_size = calc_effective_population_size(coalescent_densities)
    
     # Ensure output folder exists
    if os.path.exists(output):
        click.echo("[INFO] Output folder already exists.")
    else:
        os.makedirs(output, exist_ok=True)
        click.echo("[INFO] Created output folder.")

    with open(f"{output}/N{popsize}_U{mutrate}_s{selcoef}.out", "w") as f:
        header = ["popsize", "selcoef", "mutrate", "velocity", "density", "effective_pop_size", "time"]
        f.write("\t".join(header) + "\n")

        str_coalescent_densities = ','.join(map(str, coalescent_densities))
        str_effective_pop_size = ','.join(map(str, effective_pop_size))
        str_t = ','.join(map(str, t))

        line = [str(popsize), str(selcoef), str(mutrate), str(v), str_coalescent_densities, str_effective_pop_size, str_t]
        f.write("\t".join(line) + "\n")


def read_params(df):
    """
    Reads the parameters from the results file and converts them to variables.
    
    df: DataFrame containing the results from the forward-simulation.

    Returns:
    N: Population size.
    U: Mutation rate.
    sel_coef: Selection coefficient.
    v: Velocity of the wave.
    """
    N = int(df["popsize"].iloc[0])
    U = df["mutrate"].iloc[0]
    sel_coef = df["selcoef"].iloc[0]
    v = df["velocity"].iloc[0]
    return N, U, sel_coef, v

def create_transition_matrix(lineage_dist, N, U, sel_coef, v, delta_t):
    """ Creates the transition matrix for the backward simulation. 
    mut_burden_profile: List of mutation burden profiles for each time point. 
    N: Population size. U: Mutation rate. sel_coef: Selection coefficient. 
    
    Returns: transition_matrix: A 2D numpy array representing the transition probabilities between states. """ 
    num_states = len(lineage_dist)

    eps = 1e-15

    transition_rate_matrix = np.zeros((num_states, num_states), dtype=np.longdouble)

    # calculate off-diagonal rates
    for k in range(num_states):

        if lineage_dist[k] == 0:
           continue  # skip states with zero lineage distribution to avoid division by zero
            
        if k > 0:  # ensure we don't divide by zero
            transition_rate_matrix[k, k-1] = U * lineage_dist[k-1] / lineage_dist[k]

        if k < num_states - 1:  # ensure we don't divide by zero
            transition_rate_matrix[k, k+1] = U  * v * lineage_dist[k+1] / lineage_dist[k]

    # fill diagonals by enforcing row sums to zero
    row_sums = transition_rate_matrix.sum(axis=1)
    transition_rate_matrix[np.arange(num_states), np.arange(num_states)] -= row_sums
    scaling_factor = max(np.max(-np.diag(transition_rate_matrix)), eps)
    transition_rate_matrix_scaled = transition_rate_matrix / scaling_factor

    transition_prob_matrix = np.eye(num_states) + transition_rate_matrix_scaled * delta_t
    return transition_prob_matrix

def calc_next_lineage_dist(lineage_dist, N, U, sel_coef, v, delta_t):
    """
    Retrieves the mutational burden profile at a specific time point.
    
    mut_burden_profile: List of mutation burden profiles for each time point.
    t: The time point for which to retrieve the mutational burden profile.

    Returns:
    A numpy array representing the mutation burden profile at the specified time point.
    """
    p = create_transition_matrix(lineage_dist, N, U, sel_coef, v, delta_t)
    return np.dot(lineage_dist, p)

def calc_coalescent_rates_per_class(mut_burden_profile, lineage_dist, N):
    """
    Calculates the coalescent rates for a given mutational burden profile and parameters.
    
    mut_burden_profile: Initial mutational burden profiles.
    N: Population size.
    U: Mutation rate.
    sel_coef: Selection coefficient.
    v: Velocity of the wave.
    t: The time point for which to calculate the coalescent rates.

    Returns:
    A numpy array representing the coalescent rates at the specified time point.
    """

    return np.square(lineage_dist) / (N * mut_burden_profile)


def calc_effective_coalescent_rates(mut_burden_profile, N, U, sel_coef, v, delta_t, max_t=5000):
    """
    Calculates the effective coalescent rate by summing over all classes.
        
    mut_burden_profile: Initial mutational burden profiles.
    N: Population size.
    U: Mutation rate.
    sel_coef: Selection coefficient.
    v: Velocity of the wave.
    t: vector of time points for which to calculate the effective coalescent rates.
    
    Returns:
    A numpy array representing the overall coalescent rate up to the specified time point.
    """
    threshold = 1e-40

    num_steps = int(max_t / delta_t)

    overall_coalescent_rates = []
    survival_prob = np.zeros_like(mut_burden_profile, dtype=np.longdouble)  # initialize survival probabilities to 0
    lineage_dist = mut_burden_profile.copy()
    for _ in range(num_steps): 
        #lineage_dist[lineage_dist < threshold] = 0.0                                                 
        # Calculate the lineage distribution at the current time point,
        # let mutation and fixation happen.
        lineage_dist = calc_next_lineage_dist(lineage_dist, N, U, sel_coef, v, delta_t)   

        # Calculate the coalescent rates for each class at the current time point.
        coalescent_rates_per_class = calc_coalescent_rates_per_class(mut_burden_profile,
                                                                    lineage_dist, N)

        # Update the survival probabilities for each class based on the coalescent rates.
        survival_prob += np.log1p(-coalescent_rates_per_class*delta_t)

        # Update the lineage distribution by removing the lineages that coalesce at this time point.
        lineage_dist *= (1-np.exp(survival_prob)*coalescent_rates_per_class*delta_t)
        
        # Calculate the overall coalescent rate by summing over all classes
        overall_coalescent_rate = np.sum(coalescent_rates_per_class)
        overall_coalescent_rates.append(overall_coalescent_rate)

    return np.array(overall_coalescent_rates, dtype=np.longdouble)

def calc_effective_population_size(coalescent_rates: np.array):
    """
    Calculates the effective population size from the coalescent rates.
    
    coalescent_rates: A list of coalescent rates over time.

    Returns:
    A numpy array representing the effective population size over time.
    """
  
    return 1 / coalescent_rates

if __name__ == "__main__":
    write_file()