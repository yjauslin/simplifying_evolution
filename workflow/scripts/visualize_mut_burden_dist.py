import click
import numpy as np
import pandas as pd
import os
import re
import glob
import matplotlib.pyplot as plt
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
    
    Accepts sigma as a ratio (e.g., 0.25, 0.5) and calculates the absolute sd on-the-fly,
    matching files robustly despite floating point errors.
    """
    num_rows = 2
    # The number of columns is determined by whichever list has more entries
    num_cols = max(len(sel_coef), len(mut_rate))
    
    # Grab a colorblind-friendly palette using Seaborn based on the ratios
    colors = sns.color_palette("colorblind", n_colors=len(sigma))

    sns.set_context("paper", font_scale=1.05)
    sns.set_style("ticks")

    # Cluster-safe typography using Matplotlib's internal TeX parser engine
    plt.rcParams.update({
        "text.usetex": False, 
        "mathtext.fontset": "cm",        
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
    })

    # Figures correspond to column text boundaries in standard LaTeX templates
    fig_width = 426.79134 / 72.27  
    fig_height = fig_width * (num_rows / num_cols) * 1.15

    # Enforce clear grid alignment by sharing contextual row and column limits
    fig, axes = plt.subplots(figsize=(fig_width, fig_height), sharey="row", nrows=num_rows, ncols=num_cols)
    axes = np.atleast_2d(axes)
    
    # Adjusted left margin to 0.15 to ensure the supylabel has plenty of breathing room
    fig.subplots_adjust(hspace=0.4, wspace=0.15, left=0.15, bottom=0.15)

    for col in range(num_cols):
        # --- Row 0: Fixed mutation rate u[0], varying selection coefficient ---
        s_top = sel_coef[min(col, len(sel_coef) - 1)]
        u_top = mut_rate[0]
        ax_top = axes[0, col]
        
        ax_top.set_title(f"$U_d$ = {u_top}; $s_{{normal}}$ = {s_top}", fontsize=9)
        
        # --- Row 1: Fixed selection coefficient s[0], varying mutation rate ---
        s_bottom = sel_coef[0]
        u_bottom = mut_rate[min(col, len(mut_rate) - 1)]
        ax_bottom = axes[1, col]
        
        ax_bottom.set_title(f"$s_{{normal}}$ = {s_bottom}; $U_d$ = {u_bottom}", fontsize=9)

        # Draw execution threads across both active horizontal frames sequentially
        for row, (current_s, current_u, ax) in enumerate([(s_top, u_top, ax_top), (s_bottom, u_bottom, ax_bottom)]):
            
            # 1. Retrieve all file candidates for this (N, U, s) combination first
            # Pattern helps narrow down search surface to minimize filesystem calls
            glob_pattern = os.path.join(input_dir, f"escsim_N{pop_size}_U{current_u}_s{current_s}_sd*.out")
            matching_files = glob.glob(glob_pattern)
            
            # Loop through the user-specified standard deviation ratios (e.g., 0.25)
            for idx, ratio in enumerate(sigma):
                # Calculate what the standard deviation should mathematically be
                target_sd = ratio * current_s
                
                # 2. Iterate matches and perform floating point-safe check
                target_filepath = None
                for filepath in matching_files:
                    # Match decimal, negative, or scientific notation in the filename
                    match = re.search(r"_sd([\deE.+-]+)\.out", os.path.basename(filepath))
                    if match:
                        file_sd = float(match.group(1))
                        # Float-safe validation check
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

    # Use unified global labels
    fig.supylabel("Frequency", fontsize=11, x=0.01)
    fig.supxlabel("Mutational burden class $h_k$", fontsize=11, y=0.02)

    # Isolated legend declarations to prevent plot occlusion
    axes[0, 0].legend(loc='best', frameon=False, fontsize=8)
    axes[1, 0].legend(loc='best', frameon=False, fontsize=8)

    os.makedirs(output_dir, exist_ok=True)
    fig.tight_layout(pad=0.2)
    
    output_path = os.path.join(output_dir, "mut_burden_dist.jpg")
    fig.savefig(output_path, dpi=600, bbox_inches='tight')
    click.echo(f"[INFO] Plot layout finalized and written to {output_path}")

if __name__ == "__main__":
    main()