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
@click.argument('sel_coef', required=True, type=float)
@click.argument('mut_rate', required=True, type=float)
@click.option('--input_folder', '-i', default='tmp/results', required=False, type=str,
              help='Input folder containing comparison files.')
@click.option('--output', '-o', default='tmp/results', required=False,
              help='Output folder')
@click.option('--type', '-t', is_flag=True, default=False,
              help='If false calculates effective selection coefficient. If true calculates effective mutation rate.')

def get_min_kolmogorov(pop_size, sel_coef, mut_rate, input_folder, output, type):
    # Create output directory if it doesn't exist
    os.makedirs(output, exist_ok=True)

    # Determine the prefix for the output file based on the type of comparison
    prefix = "s" if not type else "U"

    # Construct the input file pattern based on the provided parameters
    input_file = os.path.join(
        input_folder,
        f"{prefix}_N{pop_size}_U{mut_rate}_s{sel_coef}_sd*.txt"
    )

    # Use glob to find all files matching the input pattern
    files = glob.glob(input_file)

    if not files:
        click.echo("No files found for the given parameters.")
        return

    # Construct the output file path
    output_file = os.path.join(output, f"{prefix}_N{pop_size}_U{mut_rate}_s{sel_coef}.txt")

    for f in files:
        df = pd.read_csv(f, sep="\t")

        # Extract the standard deviation (sd) from the filename using regex
        match = re.search(r"_sd([0-9.eE+-]+)\.txt$", f)
        sd = float(match.group(1))

        # Find the index of the minimum Kolmogorov value in the dataframe
        idx = df["Kolmogorov"].idxmin()

        if pd.isna(idx):
            click.echo(f"Warning: No valid numbers found in 'Kolmogorov' column for {f}. Skipping.")
            continue

        # Extract the minimum Kolmogorov value and corresponding selection coefficient and mutation rate
        min_kolmogorov = df.loc[idx, "Kolmogorov"]
        s_eff = df.loc[idx, "s"]
        mutation_rate = df.loc[idx, "U"]
        rmse = df.loc[idx, "RMSE"]

        with open(output_file, "a") as f_out:
            # Write header if the file is empty
            if f_out.tell() == 0:
                f_out.write("s\tsd\tsd/s\tU\tmin_Kolmogorov\tRMSE\n")
            # Write the results to the output file
            f_out.write(f"{s_eff}\t{sd}\t{sd/sel_coef}\t{mutation_rate}\t{min_kolmogorov}\t{rmse}\n")
        click.echo(f"File saved to: {output_file}")

if __name__ == "__main__":
    get_min_kolmogorov()
        




