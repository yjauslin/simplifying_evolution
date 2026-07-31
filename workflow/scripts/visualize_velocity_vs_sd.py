import click
import numpy as np
import pandas as pd
import os
import re
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

# Assuming this custom module exists in your path
from calc_coalescent_density import get_mean_velocity

@click.command()
@click.option('--input-dir', '-i', type=click.Path(exists=True, file_okay=False, dir_okay=True), required=True,
              help='Path to the directory containing the data files.')
@click.option('--output-dir', '-o', type=click.Path(file_okay=False, dir_okay=True), required=True,
              help='Path to the directory where the output files will be saved.')
@click.option('--sel_coef', '-s', required=False, default=(0.001,), type=float, multiple=True, help='Selection coefficient(s)')
@click.option('--mut_rate', '-u', required=False, default=(0.006,), type=float, multiple=True, help='Mutation rate(s)')
@click.argument('pop_size', required=False, type=int, default=5000)
def main(input_dir, output_dir, sel_coef, mut_rate, pop_size):
    """
    Plot grid of velocity vs. standard deviation from escsim output files.
    Generates three figures:
    1. Selection coefficient series (s_velocity.jpg)
    2. Mutation rate series (u_velocity.jpg)
    3. Combined side-by-side plot with panels 'a' and 'b' (combined_velocity.jpg)
    """
    os.makedirs(output_dir, exist_ok=True)
    base_dir = Path(input_dir)

    sel_coef = list(sel_coef)
    mut_rate = list(mut_rate)

    if len(sel_coef) == 0 or len(mut_rate) == 0:
        raise click.UsageError("You must provide at least one --sel_coef AND one --mut_rate.")

    # Styling
    sns.set_context("paper")
    sns.set_style("ticks")

    plt.rcParams.update({
        "text.usetex": True,
        "mathtext.fontset": "cm",
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
    })

    LABEL_FONTSIZE = 11
    LEGEND_FONTSIZE = 11
    TITLE_FONTSIZE = 11
    TICK_FONTSIZE = 9

    # Figure dimensions matching LaTeX column width
    fig_width = 426.79134 / 72.27  
    fig_height = fig_width / 1.618

    # Individual figures
    fig_s, ax_s = plt.subplots(figsize=(fig_width, fig_height))
    fig_u, ax_u = plt.subplots(figsize=(fig_width, fig_height))

    # Combined figure
    fig_comb, (ax_comb_a, ax_comb_b) = plt.subplots(1, 2, figsize=(fig_width, fig_height))

    available_markers = ['o', '^', 's', 'D', 'v', 'p', '*']

    def process_data_series(loop_list, is_s_mode, ax_indiv, ax_comb):
        """Helper to load data, sort parameter lists, and plot to both individual and combined axes."""
        
        # Sort ascending for s, descending for U
        sorted_list = sorted(loop_list) if is_s_mode else sorted(loop_list, reverse=True)
        colors = sns.color_palette("colorblind", len(sorted_list))
        
        total_processed = 0

        for i, val in enumerate(sorted_list):
            marker = available_markers[i % len(available_markers)]
            color = colors[i]

            if is_s_mode:
                pattern = f"escsim_N{pop_size}_U{mut_rate[0]}_s{val}_sd*.out"
                label_text = f'${val}$'
                fixed_s = val
            else:
                pattern = f"escsim_N{pop_size}_U{val}_s{sel_coef[0]}_sd*.out"
                label_text = f'${val}$'
                fixed_s = sel_coef[0]

            matching_files = list(base_dir.glob(pattern))
            x_data, y_data = [], []

            for file in matching_files:
                match = re.search(r"sd(\d+(?:\.\d+)?)\.out", file.name)
                if match:
                    sd = float(match.group(1))
                    sd_s = sd / fixed_s if fixed_s != 0 else 0

                    try:
                        df = pd.read_csv(file, sep="\t")
                        mean_velocity = get_mean_velocity(df)
                        x_data.append(sd_s)
                        y_data.append(mean_velocity)
                        total_processed += 1
                    except Exception as e:
                        click.echo(f"Warning: Failed to process file {file.name}. Error: {e}", err=True)

            if x_data:
                # Plot to individual figure
                sns.scatterplot(
                    x=x_data, y=y_data, ax=ax_indiv,
                    color=color, marker=marker, s=25, label=label_text
                )
                # Plot to combined figure
                sns.scatterplot(
                    x=x_data, y=y_data, ax=ax_comb,
                    color=color, marker=marker, s=20, label=label_text
                )

        # Apply formatting to individual axis
        ax_indiv.set_xlabel(r"$\sigma/s_{normal}$", fontsize=LABEL_FONTSIZE, labelpad=6)
        ax_indiv.set_ylabel(r"Relative Click Rate", fontsize=LABEL_FONTSIZE, labelpad=6)
        ax_indiv.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
        title_text = r"$s_{normal}$" if is_s_mode else r"$U_{normal}$"
        ax_indiv.legend(title=title_text, title_fontsize=TITLE_FONTSIZE, fontsize=LEGEND_FONTSIZE, loc="best", frameon=False)

        # Apply formatting to combined axis panel
        ax_comb.set_xlabel(r"$\sigma/s_{normal}$", fontsize=LABEL_FONTSIZE, labelpad=4)
        ax_comb.set_ylabel(r"Relative Click Rate", fontsize=LABEL_FONTSIZE, labelpad=4)
        ax_comb.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
        ax_comb.legend(title=title_text, title_fontsize=TITLE_FONTSIZE, fontsize=LEGEND_FONTSIZE, loc="best", frameon=False)

        return total_processed

    # 1. Process Selection Coefficient series
    count_s = process_data_series(sel_coef, is_s_mode=True, ax_indiv=ax_s, ax_comb=ax_comb_a)

    # 2. Process Mutation Rate series
    count_u = process_data_series(mut_rate, is_s_mode=False, ax_indiv=ax_u, ax_comb=ax_comb_b)

    if count_s == 0 and count_u == 0:
        click.echo("Error: No matching escsim data files found. Plots were not saved.", err=True)
        return

    # Add panel labels 'a' and 'b' to the combined figure
    ax_comb_a.text(-0.15, 1.05, r'\textbf{a}', transform=ax_comb_a.transAxes, fontsize=12, fontweight='bold', va='top', ha='right')
    ax_comb_b.text(-0.15, 1.05, r'\textbf{b}', transform=ax_comb_b.transAxes, fontsize=12, fontweight='bold', va='top', ha='right')

    # Save Individual Fig 1: Selection Coefficient
    fig_s.tight_layout(pad=0.1)
    fig_s.savefig(os.path.join(output_dir, "s_velocity.jpg"), dpi=600, bbox_inches='tight')

    # Save Individual Fig 2: Mutation Rate
    fig_u.tight_layout(pad=0.1)
    fig_u.savefig(os.path.join(output_dir, "u_velocity.jpg"), dpi=600, bbox_inches='tight')

    # Save Combined Fig 3
    fig_comb.tight_layout(pad=0.5)
    fig_comb.savefig(os.path.join(output_dir, "combined_velocity.jpg"), dpi=600, bbox_inches='tight')

    click.echo(f"Successfully generated 's_velocity.jpg', 'u_velocity.jpg', and 'combined_velocity.jpg' in {output_dir}")

if __name__ == '__main__':
    main()