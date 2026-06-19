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
    os.makedirs(output, exist_ok=True)

    sel_coef = list(sel_coef)
    mut_rate = list(mut_rate)

    # enforce mutual exclusivity
    if len(sel_coef) > 1 and len(mut_rate) > 1:
        raise click.UsageError(
            "You cannot provide multiple --sel_coef AND multiple --mut_rate at the same time."
        )

    if len(sel_coef) == 0 and len(mut_rate) == 0:
        raise click.UsageError(
            "You must provide at least one of --sel_coef or --mut_rate."
        )

    fig, ax = plt.subplots(figsize=(10, 6))


    colors = sns.color_palette("hsv", len(sel_coef) if len(sel_coef) > 1 else len(mut_rate))
    
    if len(sel_coef) > len(mut_rate) or (len(sel_coef) == len(mut_rate) == 1):
        for i in range(len(sel_coef)):
            s = sel_coef[i]
            input_file = os.path.join(
                input_folder,
                f"s_N{pop_size}_U{mut_rate[0]}_s{s}.txt")
            df = pd.read_csv(input_file, sep="\t")

            sns.scatterplot(x='sd/s', y='s', data=df, ax=ax, label=f's={s}', color=colors[i])
        ax.set_xlabel("sd/s")
        ax.set_ylabel("effective selection coefficient")
        ax.legend(title="Selection coefficient")
        fig.savefig(os.path.join(output, f"effective_selection_coefficient.jpg"))
    
    else:
        for i in range(len(mut_rate)):
            u = mut_rate[i]
            input_file = os.path.join(
                input_folder,
                f"U_N{pop_size}_U{u}_s{sel_coef[0]}.txt")
            df = pd.read_csv(input_file, sep="\t")

            sns.scatterplot(x='sd/s', y='U', data=df, ax=ax, label=f'U={u}', color=colors[i])
        ax.set_xlabel("sd/s")
        ax.set_ylabel("effective mutation rate")
        ax.legend(title="Mutation rate")
        fig.savefig(os.path.join(output, f"effective_mutation_rate.jpg"))

if __name__ == "__main__":
    visualize_sd_vs_s()