from pathlib import Path
import re

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd
import click

def parse_filename(path):
    """
    Function to find all result files and return their parameters.
    """
    match = re.search(r"N(\d+)_U([\deE.+-]+)_s([\deE.+-]+)", path.stem)
    if match is None:
        raise ValueError(f"Filename {path.name} does not match expected pattern.")
    pop_size = int(match.group(1))
    mut_rate = float(match.group(2))
    sel_coef = float(match.group(3))
    return pop_size, mut_rate, sel_coef

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

    # get various parameter configuration for every file
    file_info = []
    for f in results_folder.glob("*.out"):
        pop_size, mut_rate, sel_coef = parse_filename(f)
        file_info.append((f, pop_size, mut_rate, sel_coef))

    # sort the files first by population size then mutation rate and finally
    # selection coefficient
    file_info.sort(key=lambda x: (x[1], x[2], x[3]))

    files = [x[0] for x in file_info]

    click.echo("[INFO] Starting visualization of coalescent densities"
                " and effective population sizes")

    # initializing result plot for coalescent densities
    fig1, axes1 = plt.subplots(3, 4, figsize=(20, 12),
                               sharey = "row", constrained_layout = True)
    axes1 = axes1.flatten()

    # initializing result plot for effective population sizes
    fig2, axes2 = plt.subplots(3, 4, figsize=(16, 10),
                           sharey = "row", constrained_layout = True)
    axes2 = axes2.flatten()

    for i, file_path in enumerate(files):

        # reading in result file with simulation results and parameters
        df = pd.read_csv(file_path, sep="\t")
        pop_size = df["popsize"].iloc[0]
        mut_rate = df["mutrate"].iloc[0]
        sel_coef = df["selcoef"].iloc[0]
        density = np.array(df["density"].iloc[0].split(","), dtype=float)
        effective_pop_size = np.array(df["effective_pop_size"].iloc[0].split(","), dtype=float)
        time = np.array(df["time"].iloc[0].split(","), dtype=float)

        # visualizing coalescent densities
        sns.lineplot(x=time, y=density, ax=axes1[i])
        axes1[i].set_title(f"N={pop_size}, U={mut_rate}, s={sel_coef}")
        axes1[i].set_xlabel("Time (Generations)")
        axes1[i].set_ylabel("Coalescent Density")

        # converting y-axis-ticks to scientific format
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((0, 0))

        axes1[i].yaxis.set_major_formatter(formatter)

        # visualizing effective population sizes
        sns.lineplot(x=time, y=effective_pop_size, ax=axes2[i])
        axes2[i].set_title(f"N={pop_size}, U={mut_rate}, s={sel_coef}")
        axes2[i].set_xlabel("Time (Generations)")
        axes2[i].set_ylabel("Effective Population Size ($N_E$)")

    # save figures to chosen output folder
    fig1.savefig(f"{output}/coalescent_density.pdf")
    fig2.savefig(f"{output}/effective_population_size.pdf")

    click.echo(f"[INFO] Finished visualizations and saved to "
               f"{output}/coalescent_density.pdf and "
               f"{output}/effective_population_size.pdf")

if __name__ == "__main__":
    visualizing_densities()
