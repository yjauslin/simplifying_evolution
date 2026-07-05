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


    df = pd.read_csv(input, sep='\t')

    sns.set_context("talk")

    fig, ax = plt.subplots(figsize=(10, 6))

    # Set x-axis label based on the axis_type flag
    if axis_type:
        x_col = 'U'
        plt.xlabel('mutation rate')
    else:
        x_col = 's'
        plt.xlabel('selection coefficient')

    # Plot RMSE and Kolmogorov Smirnov statistics
    sns.scatterplot(x=x_col, y='RMSE', data=df, ax=ax, label='RMSE')
    sns.scatterplot(x=x_col, y='Kolmogorov', data=df, ax=ax, label='Kolmogorov Smirnov')

    plt.ylabel('RMSE / max |D|')
    plt.xticks(rotation=45)
    plt.tight_layout()

    os.makedirs(output, exist_ok=True)
    save_path = os.path.join(output, output_filename)
    plt.savefig(save_path)
    plt.close()
    
    click.echo(f"Saved plot successfully to: {save_path}")

if __name__ == '__main__':
    visualize_rmse()