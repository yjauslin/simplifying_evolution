import click
import numpy as np
import pandas as pd
import os
import re
import glob
import matplotlib.pyplot as plt
import seaborn as sns


@click.command()
@click.argument('pop_size', required=False, type=int, default=5000)
@click.option('--sel_coef', '-s', required=False, default=(0.001,), type=float, multiple=True, help='Selection coefficient')
@click.option('--mut_rate', '-u', required=False, default=(0.006,), type=float, multiple=True, help='Mutation rate')
@click.option('--input_folder', '-i', default='tmp/results', required=False, type=str,
              help='Input folder containing comparison files.')
@click.option('--output', '-o', default='tmp/results', required=False,
              help='Output folder')
def visualize_sd_vs_s(pop_size, sel_coef, mut_rate, input_folder, output):
    # Create output directory if it doesn't exist
    os.makedirs(output, exist_ok=True)

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
    
    sns.set_context("paper")
    sns.set_style("ticks")

    plt.rcParams.update({
        "text.usetex": False,
        "mathtext.fontset": "cm",        
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
    })

    # fig width corresponds to column with in LateX, 72.27 is the conversion factor from points to inches
    fig_width = 426.79134 / 72.27  
    fig_height = fig_width / 1.618

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Set color palette based on the number of selection coefficients or mutation rates
    colors = sns.color_palette("colorblind", len(sel_coef) if len(sel_coef) > 1 else len(mut_rate))
    markers = ['o', '^', 's']

    # if more selection coefficients than mutation rates are provided, we will plot effective selection coefficient vs sd/s, otherwise we will plot effective mutation rate vs sd/s
    if len(sel_coef) > len(mut_rate) or (len(sel_coef) == len(mut_rate) == 1):
        # for every provided selection coefficient...
        for i in range(len(sel_coef)):
            s = sel_coef[i]
            # Construct the input file path based on the provided parameters
            input_file = os.path.join(
                input_folder,
                f"s_N{pop_size}_U{mut_rate[0]}_s{s}.txt")
            df = pd.read_csv(input_file, sep="\t")

            sns.scatterplot(x='sd/s', y='min_Kolmogorov', data=df, ax=ax, color=colors[i], marker=markers[i], s=25, label=f'$s_{{normal}}={s}$')

        ax.set_xlabel(f"$\sigma/s_{{normal}}$", fontsize=11, labelpad=6)
        ax.set_ylabel(r"Kolmogorov-Smirnov", fontsize=11, labelpad=6)
        ax.tick_params(axis='both', which='major', labelsize=9)
        ax.legend(title=r"Selection coefficient", title_fontsize=10, fontsize=9, loc="best", frameon=False)
        fig.tight_layout(pad=0.1)
        fig.savefig(os.path.join(output, f"s_eff_Kolmogorov.jpg"), dpi=600, bbox_inches='tight')

    else:
        # for every provided mutation rate...
        for i in range(len(mut_rate)):
            u = mut_rate[i]
            input_file = os.path.join(
                input_folder,
                f"U_N{pop_size}_U{u}_s{sel_coef[0]}.txt")
            df = pd.read_csv(input_file, sep="\t")

            sns.scatterplot(x='sd/s', y='min_Kolmogorov', data=df, ax=ax, label=f'$U_d={u}$', color=colors[i], marker=markers[i], s=25)

        ax.set_xlabel(f"$\sigma/s_{{normal}}$", fontsize=11, labelpad=6)
        ax.set_ylabel(r"Kolmogorov-Smirnov", fontsize=11, labelpad=6)
        ax.tick_params(axis='both', which='major', labelsize=9)
        ax.legend(title=r"Mutation rate", title_fontsize=10, fontsize=9, loc="best", frameon=False)
        fig.tight_layout(pad=0.1)
        fig.savefig(os.path.join(output, f"U_eff_Kolmogorov.jpg"), dpi=600, bbox_inches='tight')

if __name__ == "__main__":
    visualize_sd_vs_s()