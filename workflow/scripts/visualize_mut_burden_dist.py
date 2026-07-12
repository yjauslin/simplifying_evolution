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
        # Splits '0.01,0.05,0.1' into [0.01, 0.05, 0.1]
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
              help='Comma-separated list of standard deviation / sigma values.')
@click.option('--input-dir', '-i', type=click.Path(exists=True, file_okay=False, dir_okay=True), required=True,
              help='Path to the directory containing the data files.')
@click.option('--output-dir', '-o', type=click.Path(file_okay=False, dir_okay=True), required=True,
              help='Path to the directory where the output files will be saved.')

def main(pop_size, mut_rate, sel_coef, sigma, input_dir, output_dir):
    """
    Plot an dynamically-sized grid of mutational burden distributions from escsim output files.
    Calculates the mean across profile repetitions automatically.
    """
    # The number of rows is determined by whichever list has more entries
    num_rows = max(len(sel_coef), len(mut_rate))
    
    # Grab a colorblind-friendly palette using Seaborn
    colors = sns.color_palette("colorblind", n_colors=len(sigma))

    sns.set_context("paper")
    sns.set_style("ticks")

    # fig width corresponds to column with in LateX, 72.27 is the conversion factor from points to inches
    fig_width = 426.79134 / 72.27  
    fig_height = fig_width / 1.618

    fig, axes = plt.subplots(figsize=(fig_width, fig_height), sharey=False, nrows=num_rows, ncols=2)
    fig.subplots_adjust(hspace=0.4, wspace=0.25)

    for row in range(num_rows):
        # Gracefully handle array lookups if parameters have uneven element counts
        s_left = sel_coef[min(row, len(sel_coef) - 1)]
        u_left = mut_rate[0]
        
        s_right = sel_coef[0]
        u_right = mut_rate[min(row, len(mut_rate) - 1)]
        
        # Set dynamic subplot headers
        axes[row, 0].set_title(f"N = {pop_size}; $U_d$ = {u_left}; $s$ = {s_left}", fontsize=11)
        axes[row, 1].set_title(f"N = {pop_size}; $s$ = {s_right}; $U_d$ = {u_right}", fontsize=11)

        for col in range(2):
            ax = axes[row, col]
            current_s = s_left if col == 0 else s_right
            current_u = u_left if col == 0 else u_right
            
            # Loop through individual line parameters (sigma)
            for idx, sig in enumerate(sigma):
                filename = f"escsim_N{pop_size}_U{current_u}_s{current_s}_sd{sig}.out"
                filepath = os.path.join(input_dir, filename)
                
                if os.path.exists(filepath):
                    try:
                        df = pd.read_csv(filepath, sep='\t')
                        mean_profile = get_mean_profile(df)

                        if mean_profile is not None:
                            # Generate the x-axis alignment (0, 1, 2... up to the length of the profile)
                            x = np.arange(len(mean_profile))
                        
                            ax.plot(x, mean_profile, marker='s', markersize=3, color=colors[idx], label=f"$\sigma/s$ = {sig/current_s:.2f}")
                        else:
                            click.echo(f"Warning: Empty profile data in {filename}", err=True)
                    except Exception as e:
                        click.echo(f"Error parsing file {filename}: {e}", err=True)
                else:
                    click.echo(f"Skipping profile. File not found: {filename}", err=True)

            # Axis labels and styling adjustments
            if col == 0:
                ax.set_ylabel("Frequency", fontsize=11)
            if row == num_rows - 1:
                ax.set_xlabel("Mutational burden class $h_k$", fontsize=11)
                
            ax.grid(False)

    # Place clean legend blocks at the top layout panels
    axes[0, 0].legend(loc='best', frameon=True, fontsize=10)
    axes[0, 1].legend(loc='best', frameon=True, fontsize=10)

    fig.tight_layout(pad=0.1)
    fig.savefig(os.path.join(output_dir, f"mut_burden_dist.jpg"), dpi=600, bbox_inches='tight')
    

if __name__ == "__main__":
    main()