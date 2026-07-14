from pathlib import Path
import re
from time import time

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd
import click
import math

from visualize_mut_burden_dist import parse_multiple_floats


def parse_filename(path):
    """
    Function to find all result files and return their parameters.
    """
    match = re.search(
        r"N(\d+)_U([\deE.+-]+)_s([\deE.+-]+)(?:_sd([\deE.+-]+))?",
        path.stem
    )

    if match is None:
        raise ValueError(f"Filename {path.name} does not match expected pattern.")

    pop_size = int(match.group(1))
    mut_rate = float(match.group(2))
    sel_coef = float(match.group(3))

    sigma = match.group(4)
    sigma = float(sigma) if sigma is not None else None

    return pop_size, mut_rate, sel_coef, sigma

def read_density_file(path, prefix, pop_size=None, mut_rate=None, sel_coef=None, sigma=None):
    """
    Function to find all result files with a given prefix and return their paths sorted by parameters.
    """
    file_info = []

    # Select files based on whether they contain 'sd'
    if prefix == "fixed":
        for s in sel_coef:
            for u in mut_rate:
                file = path.glob(f"N{pop_size}_U{u}_s{s}.out")
                file_info.append((file, pop_size, u, s, None))
    else:
        for s in sel_coef:
            for u in mut_rate:
                for sd in sigma:
                    file = path.glob(f"N{pop_size}_U{u}_s{s}_sd{sd}.out")
                    file_info.append((file, pop_size, u, s, sd))

    # Sort files
    if prefix == "fixed":
        # fixed files have no sigma dimension
        file_info.sort(key=lambda x: (x[3], x[1], x[2]))
    else:
        # normal files include sigma
        file_info.sort(key=lambda x: (x[3], x[4], x[1], x[2]))

    return [x[0] for x in file_info]

def read_tree_file(path, prefix, pop_size=None, mut_rate=None, sel_coef=None, sigma=None):
    """
    Function to find all tree files with a given prefix and return their paths sorted by parameters.
    """
    file_info = []
    if prefix == "fixed":
        for s in sel_coef:
            for u in mut_rate:
                file = path.glob(f"N{pop_size}_U{u}_s{s}.txt")
                file_info.append((file, pop_size, u, s, None))
    else:
        for s in sel_coef:
            for u in mut_rate:
                for sd in sigma:
                    file = path.glob(f"N{pop_size}_U{u}_s{s}_sd{sd}.txt")
                    file_info.append((file, pop_size, u, s, sd))

    # sort the files first by selection coefficient then sigma (if normal), population size and finally
    # mutation rate
    if prefix == "fixed":
        file_info.sort(key=lambda x: (x[3], x[1], x[2]))
    else:
        file_info.sort(key=lambda x: (x[3], x[4], x[1], x[2]))

    return [x[0] for x in file_info]

def build_fixed_lookup(files):
    """
    Build dictionary:
    key = (N, U, s)
    value = file path
    """
    lookup = {}

    for f in files:
        pop_size, mut_rate, sel_coef, sigma = parse_filename(f)

        key = (pop_size, mut_rate, sel_coef)

        # keep first occurrence
        if key not in lookup:
            lookup[key] = f

    return lookup

def get_effective_selection_coefficient(pop_size, sel_coef, mut_rate, sigma):
    """
    Returns the effective selection coefficient for a given population size, selection coefficient,
    mutation rate, and standard deviation (sigma) by reading the corresponding results file.
    """
    df = pd.read_csv(f'results/min_values/s_eff/s_N{pop_size}_U{mut_rate}_s{sel_coef}.txt', sep="\t")
    
    mask = np.isclose(df['sd'], sigma)
    if not mask.any():
        raise ValueError(f"No matching sigma found for value: {sigma} in s_eff file.")
        
    return df.loc[mask, 's'].values[0]

def get_effective_mutation_rate(pop_size, sel_coef, mut_rate, sigma):
    """
    Returns the effective mutation rate for a given population size, selection coefficient,
    mutation rate, and standard deviation (sigma) by reading the corresponding results file.
    """
    df = pd.read_csv(f'results/min_values/U_eff/U_N{pop_size}_U{mut_rate}_s{sel_coef}.txt', sep="\t")

    mask = np.isclose(df['sd'], sigma)
    if not mask.any():
        raise ValueError(f"No matching sigma found for value: {sigma} in U_eff file.")

    return df.loc[mask, 'U'].values[0]

def time_intervals(nsam=40, max_tmrca=2, popsize=1, include_0=True, include_inf=True):
    """
    Generates log-spaced time intervals. 
    Scaled by popsize to convert from coalescent units to generations.
    """
    # Create the base log-spaced breaks
    my_breaks = np.zeros((nsam,))
    for my_index in range(nsam):
        my_breaks[my_index] = (
            0.1 * np.exp((my_index + 1) / nsam * np.log(1 + 10 * max_tmrca)) - 0.1
        )

    if include_0:
        my_breaks = np.array([0] + list(my_breaks))

    if include_inf:
        my_breaks = np.array(list(my_breaks) + [np.inf])

    return my_breaks * popsize

def calc_density(tmrca_values, pop_size, nbin=20, tmax=3):
    """
    Calculates the coalescent density from TMRCA values using histogram binning.
    """
    
    # define bins
    bin_edges = time_intervals(nsam=nbin, max_tmrca=tmax, popsize=pop_size)
    bin_widths = np.diff(bin_edges)
    
    # get coaldens
    time, counts = np.unique(np.array(tmrca_values), return_counts=True)

    # convert to numpy arrays
    simcoal, simweights = np.array(time), np.array(counts)

    # calculate histogram counts with weights
    hist_counts, _ = np.histogram(simcoal, bins=bin_edges, weights=simweights)
    
    # normalize to get density
    simprobs = hist_counts / (hist_counts.sum() * bin_widths)

    return simprobs, bin_widths, bin_edges

def to_latex_sci(val):
    fmt = f"{val:.1e}"
    base, exponent = fmt.split("e")
    return rf"{base} \times 10^{{{int(exponent)}}}"

@click.command()
@click.option('--pop-size', '-N', type=int, default=5000, required=False, help='Population size to match.')

@click.option('--mut-rate', '-u', type=str, required=True, callback=parse_multiple_floats,
              help='Comma-separated list of mutation rates.')

@click.option('--sel-coef', '-s', type=str, required=True, callback=parse_multiple_floats,
              help='Comma-separated list of selection coefficients.')
@click.option('--sigma', '-sd', type=str, required=True, callback=parse_multiple_floats,
              help='Comma-separated list of standard deviation / sigma values.')
@click.option('--input_folder', '-i', default='results/coalescent_densities',
            help='Input folder for wave results from forward simulations'
            '(default: results/coalescent_densities)')
@click.option('--output', '-o', default='results/escsim_figures',
              help='Output folder for wave plots summary file '
              '(default: results/escsim_figures)')
def visualizing_densities(pop_size, mut_rate, sel_coef, sigma, input_folder, output):
    """
    Visualizes coalescent densities from fixed forward simulation results and the corresponding effective selection coefficients.
    """
    # Get path to results folder
    results_folder = Path(input_folder)

    sel_coef = list(sel_coef)
    mut_rate = list(mut_rate)

    # enforce mutual exclusivity, meaning: you cannot provide multiple selection coefficients and multiple mutation rates at the same time
    if len(sel_coef) > 1 and len(mut_rate) > 1:
        raise click.UsageError(
            "You cannot provide multiple --sel_coef AND multiple --mut_rate at the same time."
        )

    if len(sel_coef) == 0 and len(mut_rate) == 0:
        raise click.UsageError(
            "You must provide at least one of --sel_coef or --mut_rate."
        )

    # The number of rows is determined by whichever list has more entries
    num_rows = max(len(sel_coef), len(mut_rate))
    num_cols = len(sigma)

    sns.set_context("paper")
    sns.set_style("ticks")

    plt.rcParams.update({
        "text.usetex": False,
        "mathtext.fontset": "cm",
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
    })

    cb_palette = sns.color_palette("colorblind")
    color_fixed = cb_palette[0]
    color_normal = cb_palette[1]     
    color_wf_fixed = '#519e8a'   

    # Fig width corresponds to column width in LaTeX (72.27 points per inch)
    fig_width = 426.79134 / 72.27  
    fig_height = fig_width * (num_rows / num_cols) * 1.2

    fig, axes = plt.subplots(figsize=(fig_width, fig_height), sharey=True, sharex=True, 
    nrows=num_rows, ncols=num_cols)
    
    # Ensure axes is always a 2D array even if num_rows or num_cols == 1
    axes = np.atleast_2d(axes)
    
    fig.subplots_adjust(hspace=0.35, wspace=0.18, top=0.88, bottom=0.15, left=0.25)

    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-3, 3))

    # Determine execution flow mode
    s_eff_mode = len(sel_coef) > len(mut_rate) or (len(sel_coef) == len(mut_rate) == 1)

    # if more selection coefficients than mutation rates are provided, we will verify effective selection coefficient,  otherwise we will verify effective mutation rate
    if s_eff_mode:
        click.echo(f"[INFO] Visualizing coalescent densities for effective selection coefficients with population size {pop_size} and mutation rate {mut_rate[0]}...")
        for row in range(num_rows):
            s = sel_coef[row]
            u = mut_rate[0]
            sd = [s * sig for sig in sigma]
            click.echo(sd)
        
            for col in range(num_cols):
                ax = axes[row, col]

                # Fetching the effective selection coefficient
                s_eff = get_effective_selection_coefficient(pop_size, s, u, sd[col])

                # Gather all available fixed files for this N and U
                fixed_dir = results_folder / "fixed"
                matching_files = list(fixed_dir.glob(f"N{pop_size}_U{u}_s*.out"))

                if not matching_files:
                    raise FileNotFoundError(f"No fixed files found matching N{pop_size}_U{u} in {fixed_dir}")

                # Parse out the 's' values from the filenames to find the closest one
                available_s = []
                for f in matching_files:
                    # Extract the string between '_s' and '.out'
                    s_str = f.stem.split("_s")[-1]
                    available_s.append(float(s_str))

                # Find the index of the closest selection coefficient
                closest_idx = np.argmin(np.abs(np.array(available_s) - s_eff))
                closest_file = matching_files[closest_idx]

                # Read the closest file
                df_fixed = pd.read_csv(closest_file, sep="\t")

                normal_dir = results_folder / "normal"
                matching_normal_files = list(normal_dir.glob(f"N{pop_size}_U*_s{s}_sd*.out"))

                if not matching_normal_files:
                    raise FileNotFoundError(f"No normal files found matching N{pop_size} and s{s} in {normal_dir}")

                target_normal_file = None

                for f in matching_normal_files:
                    # Extracts both the mutation rate (U) and standard deviation (sd) securely
                    match = re.search(r"N\d+_U([\deE.+-]+)_s[\deE.+-]+_sd([\deE.+-]+)\.out", f.name)
                    if match:
                        file_u = float(match.group(1))
                        file_sd = float(match.group(2))
                        
                        # Match exact mutation rate (u) and floating-tolerant standard deviation (sd[col])
                        if file_u == u and np.isclose(file_sd, sd[col]):
                            target_normal_file = f
                            break

                if target_normal_file is None:
                    raise FileNotFoundError(
                        f"No normal file found matching U={u} and sd={sd[col]} (checked with floating-point tolerance) in {normal_dir}"
                    )

                # Set up file paths with multiple suffixes cleanly using .with_suffix()
                out_file_path = target_normal_file
                txt_file_path = target_normal_file.with_suffix(".txt")

                if not txt_file_path.exists():
                    raise FileNotFoundError(f"Expected normal sister file missing: {txt_file_path}")

                # Read the corresponding normal files
                df_normal = pd.read_csv(out_file_path, sep="\t")
                tree_normal = pd.read_csv(txt_file_path, sep="\t")

                normal_tmrca_values = []
                for entry in tree_normal['tmrca_list']:
                    values = [float(x) for x in str(entry).split(',')]
                    normal_tmrca_values.extend(values)

                simprobs, bin_widths, bin_edges = calc_density(normal_tmrca_values, pop_size=2*pop_size, tmax=3)

                density_fixed = np.array(df_fixed["density"].iloc[0].split(","), dtype=float)
                density_normal = np.array(df_normal["density"].iloc[0].split(","), dtype=float)
                time = np.array(df_fixed["time"].iloc[0].split(","), dtype=float)

                # Visualizing coalescent densitiy estimates
                sns.lineplot(x=time, y=density_fixed, ax=ax, label="fixed", color=color_fixed, lw=1.5)
                sns.lineplot(x=time, y=density_normal, ax=ax, linestyle="--", label="normal", color=color_normal, lw=1.5)

                # visualizing simulated densities as bar plots
                ax.bar(bin_edges[:-1], simprobs, width=bin_widths, alpha=0.35, label="WF_normal", color=color_wf_fixed, align="edge")

                # Axis formatting per subplot
                ax.set_xlim(0, 13000)
                ax.yaxis.set_major_formatter(formatter)
            
                ax.set_xlabel("")
                ax.set_ylabel("")

            # Add Row Titles on the right-hand side of the grid
            right_ax = axes[row, -1]
            right_ax.text(1.05, 0.5, f"s = {s}", transform=right_ax.transAxes, 
                          rotation=-90, va='center', ha='left', fontsize=10)
    
    # if more mutation rates than selection coefficients are provided, we will verify effective mutation rate, otherwise we will verify effective selection coefficient
    else:
        click.echo(f"[INFO] Visualizing coalescent densities for effective mutation rates with population size {pop_size} and selection coefficient {sel_coef[0]}...")
        for row in range(num_rows):
            # When tracking changing mutation rates, selection coefficient is held constant at index 0
            s = sel_coef[0]
            u = mut_rate[row]
            sd = [s * sig for sig in sigma]
        
            for col in range(num_cols):
                ax = axes[row, col]

                # Fetching the effective mutation rate
                u_eff = get_effective_mutation_rate(pop_size, s, u, sd[col])
                # Gather all available fixed files for this N and s
                fixed_dir = results_folder / "fixed"
                # Using a regex-like glob approach or searching for files matching the N and s pattern
                matching_files = list(fixed_dir.glob(f"N{pop_size}_U*_s{s}.out"))

                if not matching_files:
                    raise FileNotFoundError(f"No fixed files found matching N{pop_size} and s{s} in {fixed_dir}")

                # Parse out the 'U' values from the filenames to find the closest one
                available_u = []
                for f in matching_files:
                    # Extract the string between 'N..._U' and '_s'
                    # Filename structure: N{pop_size}_U{mut_rate}_s{sel_coef}.out
                    match = re.search(r"_U([\deE.+-]+)_s", f.name)
                    if match:
                        available_u.append(float(match.group(1)))
                    else:
                        available_u.append(float('inf')) # Safeguard for unexpected formats

                # Find the index of the closest mutation rate
                closest_idx = np.argmin(np.abs(np.array(available_u) - u_eff))
                closest_file_fixed = matching_files[closest_idx]

                # Read the closest fixed file
                df_fixed = pd.read_csv(closest_file_fixed, sep="\t")

                # Gather all available normal files for this N and s
                normal_dir = results_folder / "normal"
                # Using a regex-like glob approach or searching for files matching the N and s pattern
                matching_files = list(normal_dir.glob(f"N{pop_size}_U*_s{s}_sd*.out"))

                if not matching_files:
                    raise FileNotFoundError(f"No normal files found matching N{pop_size} and s{s} in {normal_dir}")

                target_file = None

                # Iterate and match using floating-point safety boundaries
                for f in matching_files:
                # Pattern extracts both the mutation rate (U) and the standard deviation (sd) safely
                    match = re.search(r"N\d+_U([\deE.+-]+)_s[\deE.+-]+_sd([\deE.+-]+)\.out", f.name)
                    if match:
                        file_u = float(match.group(1))
                        file_sd = float(match.group(2))
        
                        # Match file_u exactly and file_sd within floating-point tolerances
                        if file_u == u and np.isclose(file_sd, sd[col]):
                            target_file = f
                            break

                if target_file is None:
                    raise FileNotFoundError(
                    f"No normal file found matching U={u} and sd={sd[col]} (checked with floating-point tolerance) in {normal_dir}"
                    )
                
                out_file_path = target_file
                txt_file_path = target_file.with_suffix(".txt")

                df_normal = pd.read_csv(out_file_path, sep="\t")
                tree_normal = pd.read_csv(txt_file_path, sep="\t")

                normal_tmrca_values = []
                for entry in tree_normal['tmrca_list']:
                    values = [float(x) for x in str(entry).split(',')]
                    normal_tmrca_values.extend(values)

                simprobs, bin_widths, bin_edges = calc_density(normal_tmrca_values, pop_size=2*pop_size, tmax=3)

                density_fixed = np.array(df_fixed["density"].iloc[0].split(","), dtype=float)
                density_normal = np.array(df_normal["density"].iloc[0].split(","), dtype=float)
                time = np.array(df_fixed["time"].iloc[0].split(","), dtype=float)

                # Visualizing coalescent densities
                sns.lineplot(x=time, y=density_fixed, ax=ax, label=f"fixed", color=color_fixed, lw=1.5)
                sns.lineplot(x=time, y=density_normal, ax=ax, linestyle="--", label=f"normal", color=color_normal, lw=1.5)

                # Visualizing simulated densities as bar plots
                ax.bar(bin_edges[:-1], simprobs, width=bin_widths, alpha=0.35, label=f"WF normal", color=color_wf_fixed, align="edge")
                
                # Axis formatting per subplot
                ax.xaxis.set_major_formatter(formatter)
                ax.yaxis.set_major_formatter(formatter)
                ax.set_xlim(0, 13000)
            
                ax.set_xlabel("")
                ax.set_ylabel("")

            # Add Row Titles on the right-hand side of the grid (Labeling U_d instead of s)
            right_ax = axes[row, -1]
            right_ax.text(1.05, 0.5, f"$U_d$ = {u}", transform=right_ax.transAxes, 
                          rotation=-90, va='center', ha='left', fontsize=11)
    
    # Add Column Titles on the top of the grid
    for col in range(num_cols):
        axes[0, col].text(0.5, 1.12, f"$\\sigma = {sigma[col]} \\cdot s$", 
                          transform=axes[0, col].transAxes,
                          ha="center", va="bottom", fontsize=11)

    # Establish global labels centered perfectly across all columns and rows
    fig.supylabel(r"Coalescent Density", fontsize=11, x=0.005)
    fig.supxlabel("Time (Generations)", fontsize=11, y=0.01)
    
    # Extract global legend handles from the first plot
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.08), ncol=num_cols, fontsize=10, frameon=False)

    # Remove individual axis legends to avoid duplicates
    for ax in axes.flat:
        if ax.get_legend() is not None:
            ax.get_legend().remove()

    # Create destination output directory if it doesn't exist
    Path(output).mkdir(parents=True, exist_ok=True)
    
    file_suffix = "s_eff" if s_eff_mode else "u_eff"
    save_path = f"{output}/coalescent_density_{file_suffix}.jpg"
    fig.tight_layout(pad=0.1)
    fig.savefig(save_path, dpi=600, bbox_inches='tight')
    click.echo(f"[INFO] Finished visualizations and saved to {save_path}")

if __name__ == "__main__":
    visualizing_densities()
