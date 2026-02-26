from pathlib import Path
from scipy.linalg import expm
from scipy import interpolate, integrate

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import click
import warnings

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
    delta_t = 1
    time_points = np.arange(0, max_t, delta_t)
    coalescent_rates = calc_effective_coalescent_rates(mut_burden_profile, popsize, mutrate, selcoef, v, time_points=time_points)
    effective_pop_size = calc_effective_population_size(coalescent_rates)
    coalescent_densities = calc_coalescent_density(effective_pop_size, time_points)
    
     # Ensure output folder exists
    if os.path.exists(output):
        click.echo("[INFO] Output folder already exists.")
    else:
        os.makedirs(output, exist_ok=True)
        click.echo("[INFO] Created output folder.")

    with open(f"{output}/N{popsize}_U{mutrate}_s{selcoef}.out", "w") as f:
        header = ["popsize", "selcoef", "mutrate", "velocity", "rates", "effective_pop_size", "density", "time"]
        f.write("\t".join(header) + "\n")

        str_coalescent_rates = ','.join(map(str, coalescent_rates))
        str_effective_pop_size = ','.join(map(str, effective_pop_size))
        str_coalescent_densities = ','.join(map(str, coalescent_densities))
        str_t = ','.join(map(str, time_points))

        line = [str(popsize), str(selcoef),
                str(mutrate), str(v),
                str_coalescent_rates, str_effective_pop_size,
                str_coalescent_densities, str_t]
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

def create_transition_matrix(mut_burden_profile, N, U, sel_coef, v):
    """ Creates the transition matrix for the backward simulation. 
    mut_burden_profile: List of mutation burden profiles for each time point. 
    N: Population size. U: Mutation rate. sel_coef: Selection coefficient. 
    
    Returns: transition_matrix: A 2D numpy array representing the transition probabilities between states. """ 
    num_states = len(mut_burden_profile)

    transition_rate_matrix = np.zeros((num_states, num_states), dtype=np.float64)

    # calculate off-diagonal rates
    for k in range(num_states):

        if mut_burden_profile[k] == 0:
           continue  # skip states with zero mutational burden to avoid division by zero
            
        if k > 0:  # ensure we don't divide by zero
            transition_rate_matrix[k, k-1] = U * mut_burden_profile[k-1] / mut_burden_profile[k]

        if k < num_states - 1:  # ensure we don't divide by zero
            transition_rate_matrix[k, k+1] = U  * v * mut_burden_profile[k+1] / mut_burden_profile[k]

    # fill diagonals by enforcing row sums to zero
    row_sums = transition_rate_matrix.sum(axis=1)
    transition_rate_matrix[np.arange(num_states), np.arange(num_states)] -= row_sums

    return transition_rate_matrix

def calc_next_lineage_dist(mut_burden_profile, t, Q):
    """
    Retrieves the mutational burden profile at a specific time point.
    
    mut_burden_profile: List of mutation burden profiles for each time point.
    t: The time point for which to retrieve the mutational burden profile.

    Returns:
    A numpy array representing the mutation burden profile at the specified time point.
    """

    return np.dot(mut_burden_profile, expm(Q*t))

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


def calc_effective_coalescent_rates(mut_burden_profile, N, U, sel_coef, v, time_points):
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

    overall_coalescent_rates = []
    lineage_dist = mut_burden_profile.copy()
    transition_rate_matrix = create_transition_matrix(mut_burden_profile, N, U, sel_coef, v)
    for t in time_points:                                                
        # Calculate the lineage distribution at the current time point,
        # let mutation and fixation happen.
        lineage_dist = calc_next_lineage_dist(mut_burden_profile, t, transition_rate_matrix)   

        # Calculate the coalescent rates for each class at the current time point.
        coalescent_rates_per_class = calc_coalescent_rates_per_class(mut_burden_profile,
                                                                    lineage_dist, N)
        
        # Calculate the overall coalescent rate by summing over all classes
        overall_coalescent_rate = np.sum(coalescent_rates_per_class)
        overall_coalescent_rates.append(overall_coalescent_rate)

    return np.array(overall_coalescent_rates, dtype=np.float64)

def calc_effective_population_size(coalescent_rates: np.array):
    """
    Calculates the effective population size from the coalescent rates.
    
    coalescent_rates: A list of coalescent rates over time.

    Returns:
    A numpy array representing the effective population size over time.
    """
  
    return 1 / coalescent_rates

def calc_coalescent_density(pop_size, time_points, n_samples=2):
    """
    Calculates the coalescent density for a given mutational burden profile and parameters.
    
    mut_burden_profile: Initial mutational burden profiles.
    N: Population size.
    U: Mutation rate.
    sel_coef: Selection coefficient.
    v: Velocity of the wave.
    t: vector of time points for which to calculate the coalescent density.

    Returns:
    A numpy array representing the coalescent density at the specified time points.
    """

    n_combinations = int(n_samples * (n_samples - 1) / 2)
    
    pop_size_t = interpolate.interp1d(
        x=time_points,
        y=pop_size,
        kind='previous',
        bounds_error=False,
        fill_value=(pop_size[0], pop_size[-1]),
        assume_sorted=True
    )

    # Calculate integral for each time point
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=integrate.IntegrationWarning)
        integral_values = np.array([
            integrate.quad(lambda s: n_combinations / pop_size_t(s), 0, t)[0]
            for t in time_points
        ])
    
    # Calculate coalescent density: psi(t) = (n choose 2) * (1/N(t)) * exp(-integral)
    coalescent_densities = np.where(
        integral_values == 0,
        n_combinations / pop_size,  # At t=0, just return the rate
        (n_combinations / pop_size_t(time_points)) * np.exp(-integral_values)
    )

    return coalescent_densities

if __name__ == "__main__":
    write_file()