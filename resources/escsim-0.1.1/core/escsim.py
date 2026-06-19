import os
import sys
import hashlib
import subprocess
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from scipy import interpolate, integrate
import sklearn.linear_model as lm
from tqdm import tqdm


def run_single_sim(cmd, sim_id, folder):
    """Run a single simulation with the given command and simulation ID."""
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # get args from result
    args = extract_args(result.args)
    mutrate = float(args.get("mutrate"))

    # Check if tree mode was ON
    is_tree_mode = args.get("WRITE_TREE") == "T"

    if result.returncode != 0:
        print(f"[ERROR] Simulation {sim_id} failed with return code {result.returncode}.", file=sys.stderr)
        print(f"[ERROR] stderr: {result.stderr}", file=sys.stderr)
        
    time, wave = parse_result(result.stderr)

    wave_raw = wave.copy()

    # Fit linear model to calculate velocity
    time_flat = np.repeat(time, wave.shape[1])
    model_full = lm.LinearRegression().fit(time_flat.reshape(-1, 1), wave.flatten("C"))
    velocity = model_full.coef_[0] / mutrate
    
    # Fit linear model to calculate the profile
    mean_burden_pred = model_full.predict(time_flat.reshape(-1, 1)).reshape(wave.shape)
    wave = wave - mean_burden_pred
    wave = wave - wave.min()
    del mean_burden_pred, time_flat, model_full

    # Make burden histogram
    bin_edges = np.arange(-0.05, np.ceil(wave.max()) + 1.5, 1)
    profile, _ = np.histogram(wave, bins=bin_edges)
    profile = np.trim_zeros(profile, trim="b")  # remove trailing zeros from back
    profile = profile/profile.sum()  # get frequency

    # Write wave data to temp file, to keep main process memory-light
    tmp_path = os.path.join(folder, f"sim_{sim_id}.tmp")
    with open(tmp_path, 'w') as f_tmp:
        for t, wave_row in zip(time, wave_raw):
            wave_str = ",".join(map(str, wave_row))
            f_tmp.write(f"{t}\t{wave_str}\n") 

    del bin_edges, _, wave, wave_raw

    return {
        "sim_id": sim_id,
        "seed": args.get("seed"),
        "popsize": args.get("popsize"),
        "selcoef": args.get("selcoef"),
        "sigma": args.get("sigma"),
        "mutrate": mutrate,
        "velocity": velocity,
        "profile": profile,
        "tmp_path": tmp_path
    }


def extract_args(args):
    """Extract simulation parameters from command arguments."""
    params = {}
    arg_iter = iter(args)
    for arg in arg_iter:
        if arg == "-seed":
            seed_value = next(arg_iter)
            params["seed"] = int(seed_value)
        elif arg == "-d":
            key_value = next(arg_iter)
            key, value = key_value.split("=")
            if key in ["popsize", "selcoef", "sigma", "mutrate", "seqlen", "burnin", "ending"]:
                if key in ["popsize", "selcoef", "sigma", "mutrate"]:
                    params[key] = float(value)
                else:
                    params[key] = int(value)
    return params


def parse_result(stderr: str):
    """Parse the stderr output from the simulation to extract relevant metrics."""
    lines = [l for l in stderr.splitlines() if l and not l.startswith("#")]

    time = np.empty(len(lines), dtype=int)
    burden0 = np.empty(len(lines), dtype=int)
    max_class = np.empty(len(lines), dtype=int)
    counts = []

    for i, line in enumerate(lines):
        _time, _burden0, _counts_str = line.split(":")
    
        burden0[i] = int(_burden0)
        time[i] = int(_time)

        _counts = np.array(_counts_str.split(","), dtype=int)
        max_class[i] = burden0[i] + len(_counts)

        counts.append(_counts)
    
    max_class = np.max(max_class)

    wave = np.zeros((len(time), max_class+1), dtype=int)
    for i, _counts in enumerate(counts):
        start_idx = burden0[i]
        wave[i, start_idx:start_idx+len(_counts)] = _counts


    # Expand wave to full population size
    flat_wave = np.array([
        np.repeat(np.arange(wave.shape[1]), row)
        for row in wave
    ], dtype=int)

    return time, flat_wave


def run_external(seeds, folder, **kwargs):
    """Run external command in parallele with different seeds and return result."""
    # Check for required parameters
    if not seeds or not isinstance(seeds, list):
        raise ValueError("A list of seeds must be provided.")
    
    for seed in seeds:
        if not isinstance(seed, int):
            raise ValueError("Each seed must be an integer.")
        
    popssize = float(kwargs.get("popsize"))
    selcoef = float(kwargs.get("selcoef"))
    mutrate = float(kwargs.get("mutrate"))
    sigma = float(kwargs.get("sigma")) if kwargs.get("sigma") is not None else 0.0
    chrmlen = int(kwargs.get("chrmlen")) if kwargs.get("chrmlen") is not None else 15000
    burnin = int(kwargs.get("burnin")) if kwargs.get("burnin") is not None else int(popssize)
    gens = int(kwargs.get("gens")) if kwargs.get("gens") is not None else int(popssize) * 2
    jobs = int(kwargs.get("jobs"))
    mode = kwargs.get("mode")


    # Double the gens if popsize is less or equal to 1000
    if popssize <= 1000:
        gens = gens * 2
        print(f"[INFO] Population size is {popssize}, doubling generations to {gens}.", file=sys.stderr)

    pars = [popssize, selcoef, mutrate, chrmlen, burnin, gens, jobs]
    for p in pars:
        if p is None:
            raise ValueError("All simulation parameters must be provided.")
    del pars, p
        
    
    # Get absolute path to the slim script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    if mode == "n":
        slim_script = os.path.join(script_dir, f"escsim_normal.slim")
    else:
        slim_script = os.path.join(script_dir, f"escsim.slim")


    sim_ids = []
    cmd_list = []
    for sid, s in enumerate(seeds):
        # Prepare the exact command as a list
        cmd = [
            "slim",
            "-seed", str(s),
            "-d", f"popsize={int(popssize)}",
            "-d", f"selcoef={float(selcoef)}",
            "-d", f"sigma={sigma:.15f}",   # Forces clean standard float text
            "-d", f"mutrate={mutrate:.15f}", # Forces clean standard float text
            "-d", f"seqlen={chrmlen}",
            "-d", f"burnin={burnin}",
            "-d", f"ending={gens}",
            "-d", f"OUTPUT_FOLDER='{folder}'",
            slim_script
            ]
        cmd_list.append(cmd)
        sim_ids.append(sid)


    with ProcessPoolExecutor(max_workers=jobs) as executor:
        futures = {executor.submit(run_single_sim, cmd, sim_id, folder): sim_id 
                   for cmd, sim_id in zip(cmd_list, sim_ids)}
        
        for future in tqdm(as_completed(futures), total=len(futures), desc="Running simulations"):
            result = future.result()
            yield result
            del result


def create_seeds(n: int, base_value: str, max_value: int = 2**32 - 1):
    """Create a list of unique seeds based on a base value."""
    # Use a string to seed hash logic for reproducibility
    seeds = set()
    counter = 0
    while len(seeds) < n:
        seed_input = f"{base_value}_{counter}".encode("utf-8")
        seed_hash = hashlib.md5(seed_input).hexdigest()
        seed_int = int(seed_hash, 16) % max_value
        seeds.add(seed_int)
        counter += 1

        if counter > n * 10:
            raise ValueError("Unable to generate enough unique seeds.")
        
    return list(seeds)


def calc_phi(N, s, U):
    if s == 0:
        return 0

    s = np.abs(s)
    return N * s * np.exp(-U/s)


def calc_popsize_sc(popsize, selcoef, mutrate, tmin, tmax, ntimes):
    """Calculate the coalescent density for the strong selection case."""
    selcoef = np.abs(selcoef)

    times = np.linspace(tmin, tmax, ntimes)
    if selcoef > 0:
        effective_popsize = popsize * np.exp( (-mutrate/selcoef) * (1-np.exp( -selcoef * times))**2)
    else:
        effective_popsize = np.full_like(times, popsize, dtype=float)

    return times, effective_popsize


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


if __name__ == "__main__":
    print("[ERROR] This module is intended to be imported, not run directly.", file=sys.stderr)
   
