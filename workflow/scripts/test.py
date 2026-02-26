from pathlib import Path
from scipy import interpolate, integrate

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings

def calc_popsize_esc(times, popsize, mutrate, velocity, profile):
    """Calculate effective population size over time using the ESC (Evolutionary Stepping Stone Coalescent) model.
    
    This function implements the stepping stone approach for calculating effective population size
    under the influence of deleterious mutations and selective sweeps.
    
    Parameters:
        times (array-like): Time steps at which to calculate effective population size (must start at 0).
        popsize (float): Fixed population size parameter used in simulations.
        mutrate (float): Deleterious mutation rate.
        velocity (float): Wave velocity scaled to mutation rate (should be between 0 and 1).
        profile (array-like): Fitness class frequency distribution (also called simhk).
            Should sum to 1.
    
    Returns:
        numpy.ndarray: Effective population size at each time step.
    
    The function:
    1. Creates a backward migration matrix based on the mutation rate, velocity, and profile
    2. Calculates lineage weight distribution over time using matrix operations
    3. Computes effective population size from the lineage weights and profile
    """
    times = np.array(times)
    profile = np.array(profile)
    
    if times[0] != 0:
        raise ValueError("Time vector must start at 0")
    
    if not np.isclose(profile.sum(), 1.0, rtol=1e-5):
        raise ValueError(f"Profile must sum to 1, but sums to {profile.sum()}")
    

    # if velocity is out of bounds, if it close to 0 make it zero
    if np.isclose(velocity, 0, atol=1e-3):
        velocity = 0.0
    elif np.isclose(velocity, 1, atol=1e-3):
        velocity = 1.0
    else:
        if velocity < 0:
            velocity = 0.0
            print(f"[WARNING] Velocity is negative ({velocity}), setting to 0.", file=sys.stderr)
        elif velocity > 1:
            velocity = 1.0
            print(f"[WARNING] Velocity is greater than 1 ({velocity}), setting to 1.", file=sys.stderr)

    if not (0 <= velocity <= 1):
        raise ValueError(f"Velocity must be between 0 and 1, but is {velocity}")
    
    # Create backward migration matrix
    number_fitness_classes = len(profile)
    migmat = np.zeros((number_fitness_classes, number_fitness_classes))
    
    # Migration rates backward in time to fitter classes (due to deleterious mutations)
    migrates_bwd_to_fitter = [
        profile[k - 1] * mutrate / profile[k]
        for k in range(1, number_fitness_classes)
    ]
    np.fill_diagonal(migmat[1:], migrates_bwd_to_fitter)
    
    # Migration rates backward in time to lower fitness classes (due to wave movement)
    migrates_bwd_to_lower = [
        profile[k + 1] * velocity * mutrate / profile[k]
        for k in range(number_fitness_classes - 1)
    ]
    np.fill_diagonal(migmat[:, 1:], migrates_bwd_to_lower)
    
    # Diagonal: probability of staying in same class (must sum to 1)
    np.fill_diagonal(migmat, 1 - migmat.sum(axis=1))
    
    # Calculate lineage weight distribution over time
    # pkt[t] gives the probability distribution of lineage being in each fitness class at time t
    time_ints = times.astype(int)
    pkt = np.array([
        profile @ np.linalg.matrix_power(migmat, t)
        for t in time_ints
    ])
    
    # Calculate effective population size at each time point
    effective_popsize = np.zeros(len(times))
    for t_idx, t in enumerate(time_ints):
        lineage_weights = pkt[t_idx]
        
        # Coalescence rate in each fitness class weighted by lineage probability
        coalescence_rates = [
            (p**2) / (popsize * fhk)
            for p, fhk in zip(lineage_weights, profile)
        ]
        summed_coal_rate = sum(coalescence_rates)
        
        effective_popsize[t_idx] = 1 / summed_coal_rate
    
    return times, effective_popsize

def estimate_coaldens_from_popsize(times, popsizes, nsam=2):
    """Using the coaldens integral (psi) to estimate coalescent densities from popsizes over time.
    
    Parameters:
        times (array-like): Time steps for which population sizes are provided.
        popsizes (array-like): Population sizes at each time step. Population size is assumed 
            constant backward in time until the next value appears.
        nsam (int): Number of samples (default: 2).
    
    Returns:
        tuple: (times, densities) where densities is the coalescent density at each time point.
    
    The coalescent density is calculated as:
        psi(t) = (n choose 2) * (1/N(t)) * exp(-integral from 0 to t of (n choose 2)/N(s) ds)
    
    where N(t) is the population size at time t.
    """
    times = np.array(times)
    popsizes = np.array(popsizes)
    
    if len(times) != len(popsizes):
        raise ValueError("times and popsizes must have the same length")
    
    # Sort by time if not already sorted
    sorted_indices = np.argsort(times)
    times = times[sorted_indices]
    popsizes = popsizes[sorted_indices]
    
    # Number of pairs
    ncomb = nsam * (nsam - 1) // 2  # equivalent to math.comb(nsam, 2)
    
    # Create interpolation function for population size
    # Use 'previous' to keep popsize constant backward in time until next value
    pop_size_t = interpolate.interp1d(
        x=times,
        y=popsizes,
        kind='previous',
        bounds_error=False,
        fill_value=(popsizes[0], popsizes[-1]),
        assume_sorted=True
    )
    
    # Calculate integral for each time point
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=integrate.IntegrationWarning)
        integral_values = np.array([
            integrate.quad(lambda s: ncomb / pop_size_t(s), 0, t)[0]
            for t in times
        ])
    
    # Calculate coalescent density: psi(t) = (n choose 2) * (1/N(t)) * exp(-integral)
    densities = np.where(
        integral_values == 0,
        ncomb / popsizes,  # At t=0, just return the rate
        (ncomb / pop_size_t(times)) * np.exp(-integral_values)
    )
    
    return times, densities

results_folder = Path("results/coalescent_densities")

df = pd.read_csv(f"{results_folder}/N5000_U0.006_s0.0.out", sep="\t")
N = df["popsize"].iloc[0]
U = df["mutrate"].iloc[0]
sel_coef = df["selcoef"].iloc[0]
v = df["velocity"].iloc[0]
coalescent_rates = np.array(df["rates"].iloc[0].split(","), dtype=float)
density = np.array(df["density"].iloc[0].split(","), dtype=float)
effective_pop_size = np.array(df["effective_pop_size"].iloc[0].split(","), dtype=float)
time = np.array(df["time"].iloc[0].split(","), dtype=float)

profile_df =  pd.read_csv("results/escsim/escsim_N5000_U0.006_s0.0.out", sep="\t")
profile = np.array(profile_df["profile"].iloc[0].split(","), dtype=float)
plt.bar(x=np.arange(len(profile)), height=profile)
plt.xlabel("Mutation Burden Class")
plt.ylabel("Frequency")
plt.title(f"N={N}, U={U}, s={sel_coef}")
plt.show()

sns.lineplot(x=time, y=coalescent_rates)
plt.title(f"N={N}, U={U}, s={sel_coef}")
plt.xlabel("Time (Generations)")
plt.ylabel("Coalescent Rate")

plt.show()

sns.lineplot(x=time, y=effective_pop_size)
plt.title(f"N={N}, U={U}, s={sel_coef}")
plt.xlabel("Time (Generations)")
plt.ylabel("Effective Population Size")

plt.show()

esc_pop_size = calc_popsize_esc(time, N, U, v, profile)[1]
esc_density = estimate_coaldens_from_popsize(time, esc_pop_size)[1]

sns.lineplot(x=time, y=density, label="Yannick")
sns.lineplot(x=time, y=esc_density, label="Stefan")
plt.title(f"N={N}, U={U}, s={sel_coef}")
plt.xlabel("Time (Generations)")
plt.ylabel("Coalescent Density")
plt.legend()

plt.show()
