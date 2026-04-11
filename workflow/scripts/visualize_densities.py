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
    file_info = []
    for f in path.glob(f"{prefix}_*.out"):
        pop_size, mut_rate, sel_coef, sigma = parse_filename(f)
        file_info.append((f, pop_size, mut_rate, sel_coef, sigma))

    # sort the files first by population size then mutation rate and finally
    # selection coefficient
    file_info.sort(key=lambda x: (x[3], x[4], x[1], x[2]))

    files = [x[0] for x in file_info]
    return files


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

    if files_fixed is None:
        click.echo("[INFO] No files with fixed selection coefficient found")
    elif files_normal is None:
        click.echo("[INFO] No files with normally distributed selection coefficient found")

    n = len(files_fixed)

    click.echo("[INFO] Starting visualization of coalescent densities"
                " and effective population sizes")
    
    ncols = 4
    nrows = math.ceil(n / ncols)

    # initializing result plot for coalescent densities
    fig1, axes1 = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows),
                               sharey = "row", constrained_layout = True)
    axes1 = axes1.flatten()

    # initializing result plot for effective population sizes
    fig2, axes2 = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows),
                           sharey = "row", constrained_layout = True)
    axes2 = axes2.flatten()

    for i, (file_fixed, file_normal) in enumerate(zip(files_fixed, files_normal)):

        # reading in result file with simulation results and parameters
        df_fixed = pd.read_csv(file_fixed, sep="\t")
        df_normal = pd.read_csv(file_normal, sep="\t")

        pop_size = df_normal["popsize"].iloc[0]
        mut_rate = df_normal["mutrate"].iloc[0]
        sel_coef = df_normal["selcoef"].iloc[0]
        sigma = df_normal["sigma"].iloc[0]

        density_fixed = np.array(df_fixed["density"].iloc[0].split(","), dtype=float)
        density_normal = np.array(df_normal["density"].iloc[0].split(","), dtype = float)
        
        effective_pop_size_fixed = np.array(df_fixed["effective_pop_size"].iloc[0].split(","), dtype=float)
        effective_pop_size_normal = np.array(df_normal["effective_pop_size"].iloc[0].split(","), dtype=float)

        time = np.array(df_fixed["time"].iloc[0].split(","), dtype=float)

        # visualizing coalescent densities
        sns.lineplot(x=time, y=density_fixed, ax=axes1[i], label="fixed")
        sns.lineplot(x=time, y=density_normal, ax=axes1[i], linestyle="--", label = "normal")
        axes1[i].set_title(f"N={pop_size}, U={mut_rate}, s={sel_coef}, sigma={sigma}")
        axes1[i].set_xlabel("Time (Generations)")
        axes1[i].set_ylabel("Coalescent Density")

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
    fig1.legend(handles, labels, loc="upper center", ncol=2, fontsize=12)

    handles, labels = axes2[0].get_legend_handles_labels()
    fig2.legend(handles, labels, loc="upper center", ncol=2, fontsize=12)

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
