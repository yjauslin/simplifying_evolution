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
    os.makedirs(output, exist_ok=True)

    prefix = "s" if not type else "U"

    input_file = os.path.join(
        input_folder,
        f"{prefix}_N{pop_size}_U{mut_rate}_s{sel_coef}_sd*.txt"
    )

    files = glob.glob(input_file)

    if not files:
        click.echo("No files found for the given parameters.")
        return

    output_file = os.path.join(output, f"{prefix}_N{pop_size}_U{mut_rate}_s{sel_coef}.txt")

    for f in files:
        df = pd.read_csv(f, sep="\t")

        match = re.search(r"_sd([0-9.eE+-]+)\.txt$", f)
        sd = float(match.group(1))


        idx = df["Kolmogorov"].idxmin()

        min_kolmogorov = df.loc[idx, "Kolmogorov"]
        selection_coeff = df.loc[idx, "s"]
        mutation_rate = df.loc[idx, "U"]

        with open(output_file, "a") as f:
            if f.tell() == 0:
                f.write("s\tsd\tsd/s\tU\tmin_Kolmogorov\n")
            f.write(f"{selection_coeff}\t{sd}\t{sd/sel_coef}\t{mutation_rate}\t{min_kolmogorov}\n")
        click.echo(f"File saved to: {output_file}")

if __name__ == "__main__":
    get_min_kolmogorov()
        




