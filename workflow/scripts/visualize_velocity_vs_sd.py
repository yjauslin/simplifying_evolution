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
@click.option('--sel_coef', '-s', required=False, default=(0.001,), type=float, multiple=True, help='Selection coefficient')
@click.option('--mut_rate', '-u', required=False, default=(0.006,), type=float, multiple=True, help='Mutation rate')
@click.argument('pop_size', required=False, type=int, default=5000)
def main(input_dir, output_dir, sel_coef, mut_rate, pop_size):
    """
    Plot a grid of velocity vs. standard deviation from escsim output files.
    """
    os.makedirs(output_dir, exist_ok=True)
    base_dir = Path(input_dir)

    sel_coef = list(sel_coef)
    mut_rate = list(mut_rate)

    # Enforce mutual exclusivity
    if len(sel_coef) > 1 and len(mut_rate) > 1:
        raise click.UsageError(
            "You cannot provide multiple --sel_coef AND multiple --mut_rate at the same time."
        )
    
    if len(sel_coef) == 0 and len(mut_rate) == 0:
        raise click.UsageError(
            "You must provide at least one of --sel_coef or --mut_rate."
        )

    # Styling
    sns.set_context("paper")
    sns.set_style("ticks")

    # Column width in LaTeX conversion (points to inches)
    fig_width = 426.79134 / 72.27  
    fig_height = fig_width / 1.618
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Determine execution mode
    plot_by_s = len(sel_coef) > len(mut_rate) or (len(sel_coef) == len(mut_rate) == 1)
    loop_list = sel_coef if plot_by_s else mut_rate
    
    colors = sns.color_palette("colorblind", len(loop_list))
    # Fallback to 'o' if we run out of defined markers
    available_markers = ['o', '^', 's', 'D', 'v', 'p', '*']
    
    total_files_processed = 0

    for i, val in enumerate(loop_list):
        marker = available_markers[i % len(available_markers)]
        
        # Build file pattern depending on mode
        if plot_by_s:
            pattern = f"escsim_N{pop_size}_U{mut_rate[0]}_s{val}_sd*.out"
            label_text = f'$s={val}$'
        else:
            pattern = f"escsim_N{pop_size}_U{val}_s{sel_coef[0]}_sd*.out"
            label_text = f'$U_d={val}$'
            
        matching_files = list(base_dir.glob(pattern))
        
        # Lists to accumulate data per group (ensures clean plotting & legends)
        x_data = []
        y_data = []

        for file in matching_files:
            filename = file.name  # FIXED: changed from Path(file_path).name
            match = re.search(r"sd(\d+(?:\.\d+)?)\.out", filename)
            
            if match:
                sd = float(match.group(1))
                s_current = val if plot_by_s else sel_coef[0]
                
                # Check for zero division just in case
                sd_s = sd / s_current if s_current != 0 else 0 
                
                try:
                    df = pd.read_csv(file, sep="\t")
                    mean_velocity = get_mean_velocity(df)
                    x_data.append(sd_s)
                    y_data.append(mean_velocity)
                    total_files_processed += 1
                except Exception as e:
                    click.echo(f"Warning: Failed to process file {filename}. Error: {e}", err=True)

        # Plot the accumulated data for this variable group all at once
        if x_data:
            sns.scatterplot(
                x=x_data, y=y_data, ax=ax, 
                color=colors[i], marker=marker, s=25, label=label_text
            )

    if total_files_processed == 0:
        click.echo("Error: No matching escsim data files found. Plot not saved.", err=True)
        return

    # Final plot adjustments
    ax.set_xlabel(r"$\sigma/s$", fontsize=11, labelpad=6)
    ax.set_ylabel("Relative Click Rate", fontsize=11, labelpad=6)
    ax.tick_params(axis='both', which='major', labelsize=9)
    
    title_text = "Selection coefficient" if plot_by_s else "Mutation rate"
    ax.legend(title=title_text, title_fontsize=10, fontsize=9, loc="best", frameon=False)
    
    fig.tight_layout(pad=0.1)
    
    out_filename = "s_velocity.jpg" if plot_by_s else "u_velocity.jpg"
    fig.savefig(os.path.join(output_dir, out_filename), dpi=600, bbox_inches='tight')
    click.echo(f"Successfully processed {total_files_processed} files and saved plot to {output_dir}/{out_filename}")

if __name__ == '__main__':
    main()