from pathlib import Path
import re

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd
import click
import math

def parse_filename(path):
    """
    Function to find all result files and return their parameters.
    """
    match = re.search(r"N(\d+)_U([\deE.+-]+)_s([\deE.+-]+)_sigma([\deE.+-]+)", path.stem)
    if match is None:
        raise ValueError(f"Filename {path.name} does not match expected pattern.")
    pop_size = int(match.group(1))
    mut_rate = float(match.group(2))
    sel_coef = float(match.group(3))
    sigma = float(match.group(4))
    return pop_size, mut_rate, sel_coef, sigma

def read_density_file(path, prefix):
    """
    Function to find all result files with a given prefix and return their paths sorted by parameters.
    """
    file_info = []
    for f in path.glob(f"{prefix}_*.out"):
        pop_size, mut_rate, sel_coef, sigma = parse_filename(f)
        file_info.append((f, pop_size, mut_rate, sel_coef, sigma))

    # sort the files first by population size then mutation rate and finally
    # selection coefficient
    file_info.sort(key=lambda x: (x[3], x[4], x[1], x[2]))

    files = [x[0] for x in file_info]
    return files

def read_tree_file(path, prefix):
    """
    Function to find all tree files with a given prefix and return their paths sorted by parameters.
    """
    file_info = []
    for f in path.glob(f"{prefix}_*.txt"):
        pop_size, mut_rate, sel_coef, sigma = parse_filename(f)
        file_info.append((f, pop_size, mut_rate, sel_coef, sigma))

    # sort the files first by population size then mutation rate and finally
    # selection coefficient
    file_info.sort(key=lambda x: (x[3], x[4], x[1], x[2]))

    files = [x[0] for x in file_info]
    return files

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

    # Scale to generations BEFORE adding 0 or inf to ensure consistent types
    my_breaks = my_breaks * popsize

    if include_0:
        my_breaks = np.insert(my_breaks, 0, 0.0)

    if include_inf:
        my_breaks = np.append(my_breaks, np.inf)

    return my_breaks

def calc_density(tmrca_values, pop_size, nbin=50, tmax=3):
    """
    Calculates the density from a list of TMRCA values.

    Parameters:
        tmrca_values (list or np.ndarray): The raw TMRCA values.
        pop_size (int): The population size (N) used for scaling bins.
        nbin (int): Number of bins.
        tmax (float): Max time in coalescent units (T/N).

    Returns:
        tuple: (simprobs, bin_widths, bin_edges)
    """
    # Convert to numpy array for speed
    tmrca_values = np.array(tmrca_values)
    
    # Define bins using the helper function
    # Note: Using 2*N because of the haploid/diploid sex assumption in your notes
    scaled_N = int(round(pop_size / 2))
    bin_edges = time_intervals(nsam=nbin, max_tmrca=tmax, popsize=2 * scaled_N)
    
    # Calculate widths for normalization (ignoring the inf bin for width calc)
    # We use a finite value for the last width to avoid division by inf
    finite_edges = bin_edges.copy()
    if np.isinf(finite_edges[-1]):
        finite_edges[-1] = finite_edges[-2] * 1.1 # Small buffer for the tail
    
    bin_widths = np.diff(finite_edges)

    # Calculate histogram
    hist_counts, _ = np.histogram(tmrca_values, bins=bin_edges)
    
    # Density calculation: count / (total * width)
    # This ensures the area under the curve integrates to 1
    total_samples = hist_counts.sum()
    if total_samples == 0:
        simprobs = np.zeros_like(bin_widths)
    else:
        simprobs = hist_counts / (total_samples * bin_widths)

    return simprobs, bin_widths, bin_edges

@click.command()
@click.option('--input_folder', '-i', default='results/coalescent_densities',
            help='Input folder for wave results from forward simulations'
            '(default: results/coalescent_densities)')
@click.option('--output', '-o', default='results/escsim_figures',
              help='Output folder for wave plots summary file '
              '(default: results/escsim_figures)')

def visualizing_densities(input_folder, output):
    """
    Takes results from forward simulations as input
    and produces two PDF files visualizing the coalescent
    densities and the effective population size.
    """

    # get path to results folder
    results_folder = Path(input_folder)

    files_fixed = read_density_file(results_folder, "fixed")
    files_normal = read_density_file(results_folder, "normal")

    tree_fixed = read_tree_file(results_folder, "fixed")
    tree_normal = read_tree_file(results_folder, "normal")

    if files_fixed is None:
        click.echo("[INFO] No files with fixed selection coefficient found")
    elif files_normal is None:
        click.echo("[INFO] No files with normally distributed selection coefficient found")
    elif tree_fixed is None:
        click.echo("[INFO] No tree files with fixed selection coefficient found")
    elif tree_normal is None:
        click.echo("[INFO] No tree files with normally distributed selection coefficient found")

    n = len(files_fixed)

    click.echo("[INFO] Starting visualization of coalescent densities"
                " and effective population sizes")
    
    ncols = 4
    nrows = math.ceil(n / ncols)

    # initializing result plot for coalescent densities
    fig1, axes1 = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows),
                               sharey = "row")
    axes1 = axes1.flatten()
    fig1.tight_layout(pad=3.0, rect=[0, 0, 1, 0.95], h_pad=5.0)

    # initializing result plot for effective population sizes
    fig2, axes2 = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows),
                           sharey = "row")
    axes2 = axes2.flatten()
    fig2.tight_layout(pad=3.0, rect=[0, 0, 1, 0.95])

    for i, (file_fixed, file_normal, tree_fixed, tree_normal) in enumerate(zip(files_fixed, files_normal, tree_fixed, tree_normal)):

        # reading in result file with simulation results and parameters
        df_fixed = pd.read_csv(file_fixed, sep="\t")
        df_normal = pd.read_csv(file_normal, sep="\t")

        tree_fixed = pd.read_csv(tree_fixed, sep="\t")
        tree_normal = pd.read_csv(tree_normal, sep="\t")

        fixed_tmrca_values = []
        normal_tmrca_values = []

        for entry in tree_fixed['tmrca_list']:
            # Split the string by comma and convert each piece to a float
            # Use strip() to handle any accidental whitespace
            values = [float(x) for x in str(entry).split(',')]
            fixed_tmrca_values.extend(values)

        for entry in tree_normal['tmrca_list']:
            # Split the string by comma and convert each piece to a float
            # Use strip() to handle any accidental whitespace
            values = [float(x) for x in str(entry).split(',')]
            normal_tmrca_values.extend(values)

        pop_size = df_normal["popsize"].iloc[0]
        mut_rate = df_normal["mutrate"].iloc[0]
        sel_coef = df_normal["selcoef"].iloc[0]
        sigma = df_normal["sigma"].iloc[0]

        simprobs, bin_widths, bin_edges = calc_density(fixed_tmrca_values, pop_size=pop_size)

        plot_edges = bin_edges.copy()
        if np.isinf(plot_edges[-1]):
            plot_edges[-1] = plot_edges[-2] * 1.1

        bin_centers = plot_edges[:-1] + (np.diff(plot_edges) / 2)

        density_fixed = np.array(df_fixed["density"].iloc[0].split(","), dtype=float)
        density_normal = np.array(df_normal["density"].iloc[0].split(","), dtype = float)
        
        effective_pop_size_fixed = np.array(df_fixed["effective_pop_size"].iloc[0].split(","), dtype=float)
        effective_pop_size_normal = np.array(df_normal["effective_pop_size"].iloc[0].split(","), dtype=float)

        time = np.array(df_fixed["time"].iloc[0].split(","), dtype=float)

        axes1[i].hist(fixed_tmrca_values, bins=plot_edges, density=True, alpha=0.3, label="WF_fixed", color='#2ca02c')
        axes1[i].hist(normal_tmrca_values, bins=plot_edges, density=True, alpha=0.3, label="WF_normal", color='gray')
        # visualizing coalescent densities
        sns.lineplot(x=time, y=density_fixed, ax=axes1[i], label="fixed", color='#1f77b4', lw = 2.5)
        sns.lineplot(x=time, y=density_normal, ax=axes1[i], linestyle="--", label = "normal", color="orange", lw = 2.5)
        
        axes1[i].set_title(f"N={pop_size}, U={mut_rate}, s={sel_coef}, sigma={sigma}")
        axes1[i].set_xlabel("Time (Generations)")
        axes1[i].set_ylabel("Coalescent Density")
        axes1[i].set_xlim(0, 5000)

        # converting y-axis-ticks to scientific format
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((0, 0))

        axes1[i].yaxis.set_major_formatter(formatter)

        # visualizing effective population sizes
        sns.lineplot(x=time, y=effective_pop_size_fixed, ax=axes2[i], label="fixed")
        sns.lineplot(x=time, y=effective_pop_size_normal, ax=axes2[i], linestyle="--", label="normal")
        axes2[i].set_title(f"N={pop_size}, U={mut_rate}, s={sel_coef}, sigma={sigma}")
        axes2[i].set_xlabel("Time (Generations)")
        axes2[i].set_ylabel("Effective Population Size ($N_E$)")

    # remove unused axes
    for j in range(n, len(axes1)):
        fig1.delaxes(axes1[j])
        fig2.delaxes(axes2[j])
    
    # create legends
    handles, labels = axes1[0].get_legend_handles_labels()
    fig1.legend(handles, labels, loc="upper right", ncol=2, fontsize=12)

    handles, labels = axes2[0].get_legend_handles_labels()
    fig2.legend(handles, labels, loc="upper right", ncol=2, fontsize=12)

    # removin dublicate legends
    for ax in axes1[:n]:
        ax.legend().remove()
    for ax in axes2[:n]:
        ax.legend().remove()

    # save figures to chosen output folder
    fig1.savefig(f"{output}/coalescent_density.pdf")
    fig2.savefig(f"{output}/effective_population_size.pdf")

    click.echo(f"[INFO] Finished visualizations and saved to "
               f"{output}/coalescent_density.pdf and "
               f"{output}/effective_population_size.pdf")

if __name__ == "__main__":
    visualizing_densities()
