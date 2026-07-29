import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import click

@click.command()
@click.option('--input', '-i', type=str, required=True, help='Input file')
@click.option('--output', '-o', type=str, required=True, help='Output folder')
@click.option('--axis_type', '-t', is_flag=True, default=False,
              help='If false x-axis corresponds to selection coefficient. If true x-axis corresponds to mutation rate.')
def visualize_rmse(input, output, axis_type):
    # Extract the filename without data-type to use for the output plot
    filename_without_ext = os.path.splitext(os.path.basename(input))[0]

    output_filename = f"{filename_without_ext}.jpg"

    sns.set_context("paper")
    sns.set_style("ticks")

    plt.rcParams.update({
        "text.usetex": True,
        "mathtext.fontset": "cm",        
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
    })

    fig_width = 426.79134 / 72.27 * 0.8
    fig_height = fig_width / 1.618


    df = pd.read_csv(input, sep='\t')

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Set x-axis label based on the axis_type flag
    if axis_type:
        x_col = 'U'
        plt.xlabel(f'$U_{{fixed}}$', fontsize=11)
    else:
        x_col = 's'
        plt.xlabel(f'$s_{{fixed}}$', fontsize=11)

    # Plot RMSE and Kolmogorov Smirnov statistics
    sns.scatterplot(x=x_col, y='RMSE', data=df, ax=ax, label='RMSE')
    sns.scatterplot(x=x_col, y='Kolmogorov', data=df, ax=ax, label='Kolmogorov-Smirnov')

    plt.ylabel('RMSE / D-Statistic', fontsize=11)
    plt.xticks(rotation=45, fontsize=9)
    plt.legend(frameon=False, fontsize = 11)
    plt.tight_layout(pad=0.2)

    os.makedirs(output, exist_ok=True)
    save_path = os.path.join(output, output_filename)
    plt.savefig(save_path, dpi=600, bbox_inches='tight')
    plt.close()
    
    click.echo(f"Saved plot successfully to: {save_path}")

if __name__ == '__main__':
    visualize_rmse()