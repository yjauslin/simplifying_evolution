from pathlib import Path

import warnings
import os

from scipy.linalg import expm
from scipy import interpolate, integrate
from math import comb

import numpy as np
import pandas as pd
import click

def format_fs_flat(value):
    v = float(value)

    if v == 0:
        return "0.0"

    # keep 15 decimals but avoid scientific drift, then strip zeros
    s = f"{v:.15f}".rstrip('0').rstrip('.')

    # ensure at least one decimal for consistency
    return s if '.' in s else s + ".0"


@click.command()
@click.argument('pop_size', type=int)
@click.argument('sel_coef', type=float)
@click.argument('mut_rate', type=float)
@click.argument('sigma', type=float, required=False, default=0.0)

@click.option('--input_folder', '-i', default='results/escsim',
            help='Input folder for simulation results (default: results/escsim)')
@click.option('--output', '-o', default='results/coalescent_densities',
              help='Output folder for coalescent densities '
              '(default: results/coalescent_densities)')
@click.option('--mode', '-m', default='f', help='Determines file names of result files, if run with fixed selection coefficient (f) files have the prefix fixed,' \
'if run with normally distributed selection coefficient (n) files have the prefix normal.')

def write_file(pop_size, sel_coef, mut_rate, sigma, input_folder, output, mode):
    """
    Reads the results from the forward-simulation,
    calculates the coalescent densities and effective population size,
    and writes the results to an output file.
    """
    file_path = Path(input_folder)
    if not file_path.exists():
        click.echo(f"Error: File {file_path} does not exist.")
        return
    
    # read in the results file and extract the mutational burden profile and velocity
    results_file = pd.read_csv(file_path, sep="\t")
    mut_burden_profile = get_mean_profile(results_file)
    velocity = get_mean_velocity(results_file)

    # set time points for which to calculate coalescent densities
    time_points = create_time_psi(tmax=13000, ndens=300)

    # calculate coalescent rates, effective population size, and coalescent densities
    time_points, coalescent_rates = calc_effective_coalescent_rates(mut_burden_profile, pop_size,
                                                       mut_rate, velocity,
                                                       time_points=time_points)
    effective_pop_size = calc_effective_population_size(coalescent_rates)
    time_points, coalescent_densities = calc_coalescent_density(effective_pop_size, time_points)
    
     # Ensure output folder exists
    if os.path.exists(output):
        click.echo("[INFO] Output folder already exists.")
    else:
        os.makedirs(output, exist_ok=True)
        click.echo("[INFO] Created output folder.")
    
    if mode == 'n':
        file_name = os.path.join(output, f"N{pop_size}_U{format_fs_flat(mut_rate)}_s{format_fs_flat(sel_coef)}_sd{format_fs_flat(sigma)}.out")
    else:
        file_name = os.path.join(output, f"N{pop_size}_U{format_fs_flat(mut_rate)}_s{format_fs_flat(sel_coef)}.out")

    # Write results to output file
    with open(file_name, "w") as f:
        if mode == 'n':
            header = ["popsize", "selcoef",
                      "mutrate", "sigma",
                      "velocity",
                      "effective_pop_size",
                      "density", "time"]
        else:
            header = ["popsize", "selcoef",
                      "mutrate", "velocity",
                      "effective_pop_size",
                      "density", "time"]

        f.write("\t".join(header) + "\n")

        str_effective_pop_size = ','.join(map(str, effective_pop_size))
        str_coalescent_densities = ','.join(map(str, coalescent_densities))
        str_t = ','.join(map(str, time_points))
        if mode == 'n':
             line = [str(pop_size), str(sel_coef),
                    str(mut_rate), str(sigma),
                    str(velocity),
                    str_effective_pop_size,
                    str_coalescent_densities, str_t]
        else:
            line = [str(pop_size), str(sel_coef),
                    str(mut_rate), str(velocity),
                    str_effective_pop_size,
                    str_coalescent_densities, str_t]
        f.write("\t".join(line) + "\n")
        click.echo(f"[INFO] Results saved to {file_name}")


def read_params(df):
    """
    Reads the parameters from the results file and converts them to variables.
    
    df: DataFrame containing the results from the forward-simulation.

    Returns:
    pop_size: Population size.
    mut_rate: Mutation rate.
    sel_coef: Selection coefficient.
    v: Velocity of the wave.
    """
    pop_size = int(df["popsize"].iloc[0])
    mut_rate = df["mutrate"].iloc[0]
    sel_coef = df["selcoef"].iloc[0]
    return pop_size, mut_rate, sel_coef

def create_time_psi(tmin=0, tmax=13000, ndens=300):
    """
    Create logarithmically spaced time points.
    Includes t=0 explicitly.
    """

    assert tmin == 0, "Present time must be zero."

    return np.array(
        [tmin]
        + list(
            np.logspace(
                np.log10(1),
                np.log10(tmax),
                num=ndens,
            )
        )
    )

def get_mean_velocity(df):
    """
    Calculates the mean velocity from a
    DataFrame containing simulation results.
    """
    if "velocity" not in df.columns:
        click.echo("[WARNING] 'velocity' column not found in results.")
        return 0.0
    
    # Calculate the mean of all rows in the velocity column
    mean_v = df["velocity"].mean()
    
    # Apply your existing bounds checks to the final mean
    if mean_v > 1:
        mean_v = 1.0
    elif mean_v < 0:
        mean_v = 0.0
        
    return mean_v

def get_mean_profile(df):
    """
    Calculates the mean profile from 
    a DataFrame containing simulation results.
    """

    all_profiles = []
    max_len = 0
    
    # First pass: collect and find the longest profile
    for profile_str in df["profile"]:
        try:
            p = np.fromstring(profile_str, sep=",", dtype=float)
        except TypeError:
            p = np.array([profile_str], dtype=float)
        all_profiles.append(p)
        if len(p) > max_len:
            max_len = len(p)
            
    if not all_profiles:
        return None
        
    # Second pass: pad ONLY the right side (the end) to match max_len
    padded_profiles = [
        np.pad(p, (0, max_len - len(p)), mode='constant', constant_values=0.0) 
        for p in all_profiles
    ]
    
    # Compute mean across the rows
    return np.mean(padded_profiles, axis=0)

def create_transition_matrix(mut_burden_profile, mut_rate, velocity):
    """ 
    Creates the transition matrix for the backward simulation. 
    mut_burden_profile: List of mutation burden profiles for each time point. 
    mut_rate: Mutation rate. 
    velocity: Velocity of the wave.
    
    Returns: transition_matrix: A 2D numpy array representing the transition probabilities
                                between states.
    """

    # Get the number of states from the length of the mutational burden profile
    num_states = len(mut_burden_profile)

    # Initialize the transition rate matrix with zeros
    transition_rate_matrix = np.zeros((num_states, num_states), dtype=np.float64)

    # calculate off-diagonal rates
    for k in range(num_states):

        # skip states with zero mutational burden to avoid division by zero
        if mut_burden_profile[k] == 0:
            continue

        # fill in the matrix with backward transition rates
        if k > 0:
            transition_rate_matrix[k, k-1] = (mut_rate * mut_burden_profile[k-1] /
                                              mut_burden_profile[k])

        # fill in the matrix with forward transition rates
        if k < num_states - 1:
            transition_rate_matrix[k, k+1] = (mut_rate * velocity *
                                              mut_burden_profile[k+1] /
                                              mut_burden_profile[k])

    # fill diagonals by enforcing row sums to zero
    row_sums = transition_rate_matrix.sum(axis=1)
    transition_rate_matrix[np.arange(num_states), np.arange(num_states)] -= row_sums

    return transition_rate_matrix

def calc_next_lineage_dist(mut_burden_profile, t, transition_rate_matrix):
    """
    Retrieves the mutational burden profile at a specific time point.
    
    mut_burden_profile: List of mutation burden profiles for each time point.
    t: The time point for which to retrieve the mutational burden profile.
    transition_rate_matrix: The transition rate matrix defining transition rates
                            between states.

    Returns:
    The dot product of the initial mutational burden profile
    and the matrix exponential of the transition rate matrix multiplied by time t,
    which gives the lineage distribution at time t.
    """

    # return the lineage distribution at time t
    # as the dot product of the initial mutational burden profile
    # and the matrix exponential of the transition rate matrix multiplied by time t
    return np.dot(mut_burden_profile, expm(transition_rate_matrix*t))

def calc_coalescent_rates_per_class(mut_burden_profile, lineage_dist, pop_size):
    """
    Calculates the coalescent rates for a given mutational burden profile and parameters.
    
    mut_burden_profile: Initial mutational burden profiles.
    lineage_dist: The distribution of lineages
                across fitness classes at a given time point.
    pop_size: Population size.


    Returns:
    A numpy array representing the coalescent rates at the specified time point.
    """
    # Create an array of zeros with the same shape/dtype as lineage_dist
    coalescent_rates = np.zeros_like(lineage_dist, dtype=np.float64)

    # Only calculate rates for classes where mut_burden_profile is strictly positive
    mask = mut_burden_profile > 0

    # calculate coalescent rates per class as
    # (lineage distribution)^2 / (population size * mutational burden profile)
    coalescent_rates[mask] = np.square(lineage_dist[mask]) / (pop_size * mut_burden_profile[mask])
    
    return coalescent_rates


def calc_effective_coalescent_rates(mut_burden_profile, pop_size, mut_rate, velocity, time_points):
    """
    Calculates the effective coalescent rate by summing over all classes.
        
    mut_burden_profile: Initial mutational burden profiles.
    pop_size: Population size.
    mut_rate: Mutation rate.
    velocity: Velocity of the wave.
    time_points: vector of time points for which to calculate the effective coalescent rates.
    
    Returns:
    A numpy array representing the overall coalescent rate up to the specified time point.
    """

    # initialize list to store overall coalescent rates at each time point
    overall_coalescent_rates = []
    # initialize lineage distribution with the initial mutational burden profile
    lineage_dist = mut_burden_profile.copy()
    # Create the transition matrix once, since it stays constant over time
    transition_rate_matrix = create_transition_matrix(mut_burden_profile, mut_rate, velocity)
    for t in time_points:

        # Calculate the lineage distribution at the current time point,
        lineage_dist = calc_next_lineage_dist(mut_burden_profile, t, transition_rate_matrix) 

        # Calculate the coalescent rates for each class at the current time point.
        coalescent_rates_per_class = calc_coalescent_rates_per_class(mut_burden_profile,
                                                                    lineage_dist, pop_size)
        
        # Calculate the overall coalescent rate by summing over all classes
        overall_coalescent_rate = np.sum(coalescent_rates_per_class)
        overall_coalescent_rates.append(overall_coalescent_rate)

    return time_points, np.array(overall_coalescent_rates, dtype=np.float64)

def calc_effective_population_size(coalescent_rates: np.array):
    """
    Calculates the effective population size from the coalescent rates.
    
    coalescent_rates: A numpy array of coalescent rates over time.

    Returns:
    A numpy array representing the effective population size over time.
    """
    # return effective population size as the inverse of the coalescent rates
    return 1 / coalescent_rates

def calc_coalescent_density(pop_size, time_points, n_samples=2):
    """
    Calculates the coalescent density for a given mutational burden
    profile and parameters. 
    pop_size: Population size. 
    time_points: vector of time points for which to calculate the coalescent density. 
    n_samples: number of samples to consider. 
    Returns: A numpy array representing the coalescent density at the specified time points. 
    """ 
    times = np.array(time_points)
    popsizes = np.array(pop_size)
    n_combinations = comb(n_samples, 2)

    if len(times) != len(popsizes):
        raise ValueError("times and popsizes must have the same length")

    # Sort by time if not already sorted
    sorted_indices = np.argsort(times)
    times = times[sorted_indices]
    popsizes = popsizes[sorted_indices]

    lambda_t = n_combinations / popsizes

    cumulative_hazard = integrate.cumulative_trapezoid(
        lambda_t,
        times,
        initial=0,
    )
    
    # Create interpolation function for population size
    # Use 'previous' to keep popsize constant backward in time until next value
    # pop_size_t = interpolate.interp1d(
    #    x=times,
    #    y=popsizes,
    #    kind='previous',
    #    bounds_error=False,
    #    fill_value=(popsizes[0], popsizes[-1]),
    #    assume_sorted=True
    # )
    
    # Calculate integral for each time point
    # with warnings.catch_warnings():
    #    warnings.filterwarnings('ignore', category=integrate.IntegrationWarning)
    #    integral_values = np.array([
    #        integrate.quad(lambda s: n_combinations / pop_size_t(s), 0, t)[0]
    #        for t in times
    #    ])

    # Calculate coalescent density: psi(t) = (n choose 2) * (1/N(t)) * exp(-integral)
    # densities = np.where(
    #    integral_values == 0,
    #    n_combinations / popsizes,  # At t=0, just return the rate
    #    (n_combinations / pop_size_t(times)) * np.exp(-integral_values)
    # )

    densities = lambda_t * np.exp(-cumulative_hazard)
    
    return times, densities
    

    

    

    
    return times, densities
if __name__ == "__main__":
    write_file()
