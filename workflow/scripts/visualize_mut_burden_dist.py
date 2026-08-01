import click
import numpy as np
import pandas as pd
import os
import re
import glob
import copy
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

from calc_coalescent_density import get_mean_profile

def parse_multiple_floats(ctx, param, value):
    """Parses a comma-separated string of numbers into a list of floats."""
    if not value:
        return []
    try:
        return [float(x.strip()) for x in value.split(',')]
    except ValueError:
        raise click.BadParameter("Arguments must be a comma-separated list of floats.")

@click.command()
@click.option('--pop_size', '-N', type=int, default=5000, required=False, help='Population size to match.')
@click.option('--mut_rate', '-u', type=str, required=True, callback=parse_multiple_floats,
              help='Comma-separated list of mutation rates.')
@click.option('--sel_coef', '-s', type=str, required=True, callback=parse_multiple_floats,
              help='Comma-separated list of selection coefficients.')
@click.option('--sigma', '-sd', type=str, required=True, callback=parse_multiple_floats,
              help='Comma-separated list of sigma/s ratios (e.g., 0.25,0.5,0.75,1.0).')
@click.option('--input-dir', '-i', type=click.Path(exists=True, file_okay=False, dir_okay=True), required=True,
              help='Path to the directory containing the data files.')
@click.option('--output-dir', '-o', type=click.Path(file_okay=False, dir_okay=True), required=True,
              help='Path to the directory where the output files will be saved.')
def main(pop_size, mut_rate, sel_coef, sigma, input_dir, output_dir):
    """
    Plot a grid of mutational burden distributions from escsim output files with 2 fixed rows.
    Row 0: Fixed mutation rate, varying selection coefficients across columns.
    Row 1: Fixed selection coefficient, varying mutation rates across columns.
    """
    sel_coef_copy = copy.deepcopy(sel_coef)
    mut_rate_copy = copy.deepcopy(mut_rate)

    # Sort selection coefficients ascending and mutation rates descending
    sel_coef_sorted = sorted(sel_coef)
    mut_rate_sorted = sorted(mut_rate, reverse=True)

    num_rows = 2
    num_cols = max(len(sel_coef_sorted), len(mut_rate_sorted))
    
    colors = [
        "#332288",  # Indigo
        "#88CCEE",  # Cyan
        "#44AA99",  # Teal
        "#117733",  # Green
        "#999933",  # Olive
        "#DDCC77",  # Sand / Yellow
        "#CC6677",  # Rose / Coral
        "#882255",  # Wine
        "#AA4499",  # Purple
        "#661100",  # Dark Brown
    ]

    sns.set_context("paper", font_scale=1.05)
    sns.set_style("ticks")

    plt.rcParams.update({
        "text.usetex": True, 
        "mathtext.fontset": "cm",        
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
    })

    fig_width = 426.79134 / 72.27  
    fig_height = fig_width * (num_rows / num_cols) * 1.25

    # Independent y-axes per subplot
    fig, axes = plt.subplots(figsize=(fig_width, fig_height), nrows=num_rows, ncols=num_cols)
    axes = np.atleast_2d(axes)

    for col in range(num_cols):
        # --- Row 0: Fixed mutation rate, varying selection coefficient (ascending) ---
        s_top = sel_coef_sorted[min(col, len(sel_coef_sorted) - 1)]
        u_top = mut_rate_copy[0]
        ax_top = axes[0, col]
        
        ax_top.set_title(f"$s_{{normal}}$ = {s_top}", fontsize=11, pad=6)
        
        # --- Row 1: Fixed selection coefficient, varying mutation rate (descending) ---
        s_bottom = sel_coef_copy[0]
        u_bottom = mut_rate_sorted[min(col, len(mut_rate_sorted) - 1)]
        ax_bottom = axes[1, col]
        
        ax_bottom.set_title(f"$U_{{normal}}$ = {u_bottom}", fontsize=11, pad=6)

        for row, (current_s, current_u, ax) in enumerate([(s_top, u_top, ax_top), (s_bottom, u_bottom, ax_bottom)]):
            
            # Limit y-tick density to prevent collisions across subplots
            ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4, prune=None))

            glob_pattern = os.path.join(input_dir, f"escsim_N{pop_size}_U{current_u}_s{current_s}_sd*.out")
            matching_files = glob.glob(glob_pattern)
            
            for idx, ratio in enumerate(sigma):
                target_sd = ratio * current_s
                target_filepath = None
                
                for filepath in matching_files:
                    match = re.search(r"_sd([\deE.+-]+)\.out", os.path.basename(filepath))
                    if match:
                        file_sd = float(match.group(1))
                        if np.isclose(file_sd, target_sd, rtol=1e-5, atol=1e-8):
                            target_filepath = filepath
                            break
                
                if target_filepath and os.path.exists(target_filepath):
                    try:
                        df = pd.read_csv(target_filepath, sep='\t')
                        mean_profile = get_mean_profile(df)

                        if mean_profile is not None:
                            x = np.arange(len(mean_profile))
                            ax.plot(
                                x, mean_profile, 
                                marker='s', 
                                markersize=2, 
                                color=colors[idx], 
                                label=f"$\\sigma/s$ = {ratio:.2f}",
                                lw=1.2
                            )
                        else:
                            click.echo(f"Warning: Empty profile data in {os.path.basename(target_filepath)}", err=True)
                    except Exception as e:
                        click.echo(f"Error parsing file {os.path.basename(target_filepath)}: {e}", err=True)
                else:
                    click.echo(f"Skipping profile. No file matches ratio {ratio:.2f} (sd={target_sd}) for s={current_s}", err=True)

            ax.grid(False)

    # Place y-axis label "Frequency" on the leftmost subplots so it sits to the right of plot labels ('a', 'b')
    axes[0, 0].set_ylabel("Frequency", fontsize=11, labelpad=10)
    axes[1, 0].set_ylabel("Frequency", fontsize=11, labelpad=10)

    # Position plot labels ('a', 'b') to the left of the Frequency y-axis label
    axes[0, 0].text(-0.50, 1.18, r'\textbf{a}', transform=axes[0, 0].transAxes, fontsize=11, fontweight='bold', va='top', ha='right')
    axes[1, 0].text(-0.50, 1.18, r'\textbf{b}', transform=axes[1, 0].transAxes, fontsize=11, fontweight='bold', va='top', ha='right')

    # Right-side labels flipped (reading top-to-bottom)
    ax_top_right = axes[0, -1].twinx()
    ax_top_right.set_ylabel(f"$U_{{normal}}$ = {mut_rate_copy[0]}", fontsize=11, labelpad=12, rotation=-90, va="bottom")
    ax_top_right.set_yticks([])

    ax_bottom_right = axes[1, -1].twinx()
    ax_bottom_right.set_ylabel(f"$s_{{normal}}$ = {sel_coef_copy[0]}", fontsize=11, labelpad=12, rotation=-90, va="bottom")
    ax_bottom_right.set_yticks([])

    # Shared X-label
    fig.supxlabel("Mutational burden class $h_k$", fontsize=11, y=0.01)

    # Collect legend handles and labels
    handles, labels = axes[0, 0].get_legend_handles_labels()

    if handles:
        fig.legend(
            handles, labels, 
            loc='upper center', 
            bbox_to_anchor=(0.50, 1.02), 
            ncol=len(sigma), 
            frameon=False, 
            fontsize=11
        )

    os.makedirs(output_dir, exist_ok=True)
    
    fig.tight_layout(pad=0.2, rect=[0.08, 0.04, 0.94, 0.93])
    fig.subplots_adjust(wspace=0.45, hspace=0.45)
    
    output_path = os.path.join(output_dir, "mut_burden_dist.jpg")
    fig.savefig(output_path, dpi=600, bbox_inches='tight')
    click.echo(f"[INFO] Plot layout finalized and written to {output_path}")

if __name__ == "__main__":
    main()