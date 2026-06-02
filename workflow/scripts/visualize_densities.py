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

def read_density_file(path, prefix):
    """
    Function to find all result files with a given prefix and return their paths sorted by parameters.
    """
    file_info = []

    # Select files based on whether they contain 'sd'
    if prefix == "fixed":
        files = [f for f in path.glob("*.out") if "sd" not in f.name]
    else:
        files = [f for f in path.glob("*.out") if "sd" in f.name]

    for f in files:
        pop_size, mut_rate, sel_coef, sigma = parse_filename(f)
        file_info.append((f, pop_size, mut_rate, sel_coef, sigma))

    # Sort files
    if prefix == "fixed":
        # fixed files have no sigma dimension
        file_info.sort(key=lambda x: (x[3], x[1], x[2]))
    else:
        # normal files include sigma
        file_info.sort(key=lambda x: (x[3], x[4], x[1], x[2]))

    return [x[0] for x in file_info]

def read_tree_file(path, prefix):
    """
    Function to find all tree files with a given prefix and return their paths sorted by parameters.
    """
    file_info = []
    for f in path.glob(f"{prefix}_*.txt"):

        pop_size, mut_rate, sel_coef, sigma = parse_filename(f)
        file_info.append((f, pop_size, mut_rate, sel_coef, sigma))

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

    # build lookup tables for fixed files to easily find matching files for normal files
    fixed_density_lookup = build_fixed_lookup(files_fixed)
    fixed_tree_lookup = build_fixed_lookup(tree_fixed)

    extended_files_fixed = []
    extended_tree_fixed = []

    # loop through normal files and find matching fixed files based on parameters
    for normal_file, normal_tree in zip(files_normal, tree_normal):

        pop_size, mut_rate, sel_coef, sigma = parse_filename(normal_file)

        key = (pop_size, mut_rate, sel_coef)

        if key not in fixed_density_lookup:
            raise ValueError(f"No matching fixed density file for {key}")

        if key not in fixed_tree_lookup:
            raise ValueError(f"No matching fixed tree file for {key}")

        extended_files_fixed.append(fixed_density_lookup[key])
        extended_tree_fixed.append(fixed_tree_lookup[key])

    if files_fixed is None:
        click.echo("[INFO] No files with fixed selection coefficient found")
    elif files_normal is None:
        click.echo("[INFO] No files with normally distributed selection coefficient found")
    elif tree_fixed is None:
        click.echo("[INFO] No tree files with fixed selection coefficient found")
    elif tree_normal is None:
        click.echo("[INFO] No tree files with normally distributed selection coefficient found")

    n = len(files_normal)

    click.echo("[INFO] Starting visualization of coalescent densities"
                " and effective population sizes")
    
    sns.set_context("talk")

    ncols = 4
    nrows = math.ceil(n / ncols)

    figsize_scale = 1
    # initializing result plot for coalescent densities
    fig1, axes1 = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows), layout="constrained")
    axes1 = axes1.flatten()

    # initializing result plot for effective population sizes
    fig2, axes2 = plt.subplots(nrows, ncols, figsize=(4*ncols, 4.5*nrows),
                           sharey = "row", layout="constrained")
    axes2 = axes2.flatten()

    for i, (file_fixed, file_normal, tree_fixed, tree_normal) in enumerate(zip(extended_files_fixed,
                                                                               files_normal,
                                                                               extended_tree_fixed,
                                                                               tree_normal)):

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

        simprobs_f, bin_widths_f, bin_edges_f = calc_density(fixed_tmrca_values, pop_size=2*pop_size, tmax=3)
        simprobs_n, bin_widths_n, bin_edges_n = calc_density(normal_tmrca_values, pop_size=2*pop_size, tmax=3)

        density_fixed = np.array(df_fixed["density"].iloc[0].split(","), dtype=float)
        density_normal = np.array(df_normal["density"].iloc[0].split(","), dtype = float)

        density_fixed = density_fixed
        density_normal = density_normal

        effective_pop_size_fixed = np.array(df_fixed["effective_pop_size"].iloc[0].split(","), dtype=float)
        effective_pop_size_normal = np.array(df_normal["effective_pop_size"].iloc[0].split(","), dtype=float)

        time = np.array(df_fixed["time"].iloc[0].split(","), dtype=float)

        # axes1[i].hist(fixed_tmrca_values, bins=30, density=True, alpha=0.3, label="WF_fixed", color='#2ca02c')
        # axes1[i].hist(normal_tmrca_values, bins=30, density=True, alpha=0.3, label="WF_normal", color='gray')

        # visualizing coalescent densities
        sns.lineplot(x=time, y=density_fixed, ax=axes1[i], label="fixed", color='#1f77b4', lw = 2.5)
        sns.lineplot(x=time, y=density_normal, ax=axes1[i], linestyle="--", label = "normal", color="orange", lw = 2.5)

        # visualizing simulated densities as bar plots
        # sns.histplot(fixed_tmrca_values, stat="density", bins=20, ax=axes1[i], label="WF_fixed", color = "grey", alpha=0.3)
        # sns.histplot(normal_tmrca_values, stat="density", bins=20, ax = axes1[i], label="WF_normal", color="#2ca02c", alpha=0.3)
        axes1[i].bar(bin_edges_f[:-1], simprobs_f, width=bin_widths_f, alpha=0.3, label="WF_fixed", color='#2ca02c', align="edge")
        axes1[i].bar(bin_edges_n[:-1], simprobs_n, width=bin_widths_n, alpha=0.3, label="WF_normal", color='gray', align="edge")


        # axes1[i].set_title(rf"N={pop_size}, U={mut_rate}, "
        #                   rf"$s={to_latex_sci(sel_coef)}, "
        #                   rf"\sigma={to_latex_sci(sigma)}$")
        axes1[i].set_xlabel("Time (Generations)")
        axes1[i].set_ylabel("Coalescent Density")
        axes1[i].set_xlim(0, 13000)
        axes1[i].text(0.98, 0.98,
                      rf"$N={pop_size}$" "\n"
                      rf"$U={mut_rate}$" "\n"
                      rf"$s={to_latex_sci(sel_coef)}$" "\n"
                      rf"$\sigma={to_latex_sci(sigma)}$",
                      transform=axes1[i].transAxes,
                      va="top",
                      ha="right",
                      fontsize=10
                      )

        # converting y-axis-ticks to scientific format
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((0, 0))

        axes1[i].yaxis.set_major_formatter(formatter)

        # visualizing effective population sizes
        sns.lineplot(x=time, y=effective_pop_size_fixed, ax=axes2[i], label="fixed")
        sns.lineplot(x=time, y=effective_pop_size_normal, ax=axes2[i], linestyle="--", label="normal")
        # axes2[i].set_title(rf"N={pop_size}, U={mut_rate}, "
        #                   rf"$s={to_latex_sci(sel_coef)}, "
        #                   rf"\sigma={to_latex_sci(sigma)}$")
        axes2[i].set_xlabel("Time (Generations)")
        axes2[i].set_ylabel("Effective Population Size ($N_E$)")
        axes2[i].text(0.98, 0.98,
                      rf"$N={pop_size}$" "\n"
                      rf"$U={mut_rate}$" "\n"
                      rf"$s={to_latex_sci(sel_coef)}$" "\n"
                      rf"$\sigma={to_latex_sci(sigma)}$",
                      transform=axes2[i].transAxes,
                      ha="right",
                      va="top",
                      fontsize=10
                      )

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
    fig1.savefig(f"{output}/coalescent_density.jpg")
    fig2.savefig(f"{output}/effective_population_size.jpg")

    click.echo(f"[INFO] Finished visualizations and saved to "
               f"{output}/coalescent_density.jpg and "
               f"{output}/effective_population_size.jpg")

if __name__ == "__main__":
    visualizing_densities()
