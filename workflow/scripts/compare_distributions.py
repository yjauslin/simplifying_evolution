import numpy as np
import pandas as pd
import os
import click
from scipy.integrate import cumulative_trapezoid
import glob

@click.command()
@click.argument('input_fixed', nargs=-1, required=True)
@click.option('--input_normal', '-i_n', required=True,
              help='Path to sd (normal) coalescent density file.')
@click.option('--output', '-o', required=True,
              help='Output folder.')
def compare_distributions(input_fixed, input_normal, output):
    os.makedirs(output, exist_ok=True)

    df_normal = pd.read_csv(input_normal, sep="\t")

    pop_size = df_normal["popsize"].iloc[0]
    mut_rate = df_normal["mutrate"].iloc[0]
    sigma = df_normal["sigma"].iloc[0]
    sel_coef_normal = df_normal["selcoef"].iloc[0]

    click.echo(f"Comparing distributions for N={pop_size}, U={mut_rate}, s={sel_coef_normal}, sd={sigma}...")

    time = np.array(df_normal["time"].iloc[0].split(","), dtype=float)
    density_normal = np.array(df_normal["density"].iloc[0].split(","), dtype=float)

    cdf_normal = cumulative_trapezoid(density_normal, time, initial=0)

    click.echo(f"Got {len(input_fixed)} files to compare.")


    for f in input_fixed:
        df_fixed = pd.read_csv(f, sep="\t")

        sel_coef = df_fixed["selcoef"].iloc[0]

        density_fixed = np.array(df_fixed["density"].iloc[0].split(","), dtype=float)
        cdf_fixed = cumulative_trapezoid(density_fixed, time, initial=0)

        diff = cdf_fixed - cdf_normal

        rmse = np.sqrt((diff ** 2).mean())
        max_t = np.max(np.abs(diff))

        # OUTPUT uses s from sd file (normal file)
        output_file = os.path.join(
            output,
            f"N{pop_size}_U{mut_rate}_s{sel_coef_normal}_sd{sigma}.txt"
        )

        with open(output_file, "a") as f:
            if f.tell() == 0:
                f.write("s\tRMSE\tT\n")
            f.write(f"{sel_coef}\t{rmse}\t{max_t}\n")
        click.echo("Comparison file saved to: " + output_file)


if __name__ == "__main__":
    compare_distributions()