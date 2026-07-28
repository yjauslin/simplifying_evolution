import numpy as np
import pandas as pd
import os
import click
from scipy.integrate import cumulative_trapezoid
import glob

def format_float_smart(val):
    """ Converts a float to a decimal string without scientific notation, 
    preserving its natural precision without adding unnecessary trailing zeros.
    """
    # Format to a very high maximum precision (here 15 places) to catch small numbers
    # then strip any unnecessary trailing zeros and the decimal point if it's a whole number
    s = f"{val:.15f}".rstrip('0')
    if s.endswith('.'):
        s += '0'
    return s

@click.command()
@click.option('--input_fixed', '-i_f', required=True,
              help='Space-separated paths to fixed parameter density files.')
@click.option('--input_normal', '-i_n', required=True,
              help='Path to sd (normal) coalescent density file.')
@click.option('--output', '-o', required=True,
              help='Output folder.')
@click.option('--type', '-t', is_flag=True, default=False,
              help='If false calculates effective selection coefficient. If true calculates effective mutation rate.')
def compare_distributions(input_fixed, input_normal, output, type):
    
    # Create output directory if it doesn't exist
    os.makedirs(output, exist_ok=True)

    df_normal = pd.read_csv(input_normal, sep="\t")

    # Extract parameters from the normal distribution file
    pop_size = df_normal["popsize"].iloc[0]
    mut_rate_normal = df_normal["mutrate"].iloc[0]
    sigma = df_normal["sigma"].iloc[0]
    sel_coef_normal = df_normal["selcoef"].iloc[0]

    click.echo(f"Comparing distributions for N={pop_size}, U={mut_rate_normal}, s={sel_coef_normal}, sd={sigma}...")

    # Read time and density from the normal distribution file
    time = np.array(df_normal["time"].iloc[0].split(","), dtype=float)
    density_normal = np.array(df_normal["density"].iloc[0].split(","), dtype=float)

    # Compute the cumulative distribution function (CDF) for the normal distribution
    cdf_normal = cumulative_trapezoid(density_normal, time, initial=0)

    # Split the input_fixed string into individual file paths, ensuring to strip whitespace and ignore empty strings
    file_list = [f.strip() for f in input_fixed.split(" ") if f.strip()]
    
    click.echo(f"Got {len(file_list)} files to compare.")

    # Loop through each fixed parameter density file and compute RMSE and Kolmogorov-Smirnov statistics
    for file_path in file_list:
        df_fixed = pd.read_csv(file_path, sep="\t")

        # Extract parameters from the fixed parameter density file
        sel_coef = df_fixed["selcoef"].iloc[0]
        mut_rate = df_fixed["mutrate"].iloc[0]

        # Read density from the fixed parameter density file
        density_fixed = np.array(df_fixed["density"].iloc[0].split(","), dtype=float)
        
        # Compute the cumulative distribution function (CDF) for the fixed parameter distribution
        cdf_fixed = cumulative_trapezoid(density_fixed, time, initial=0)

        # Compute the difference between the two CDFs
        diff = cdf_fixed - cdf_normal

        # Calculate RMSE and maximum absolute difference (Kolmogorov-Smirnov statistic)
        rmse = np.sqrt((diff ** 2).mean())
        max_d = np.max(np.abs(diff))

        # Prepare the output file name based on the type of comparison (s_eff=F or U_eff=T)
        if type:
        # OUTPUT uses s from sd file (normal file) and U from sd file (normal file) as well, because we are comparing to that distribution
            output_file = os.path.join(
                output,
                f"U_N{pop_size}_U{mut_rate_normal}_s{sel_coef_normal}_sd{format_float_smart(sigma)}.txt"
            )
        else:
            output_file = os.path.join(
                output,
                f"s_N{pop_size}_U{mut_rate_normal}_s{sel_coef_normal}_sd{format_float_smart(sigma)}.txt"
            )

        with open(output_file, "a") as file:
            # Write header if the file is empty
            if file.tell() == 0:
                file.write("s\tU\tRMSE\tKolmogorov\n")
            # Write the results to the output file
            file.write(f"{sel_coef}\t{mut_rate}\t{rmse}\t{max_d}\n")
    click.echo("Comparison file saved to: " + output_file)


if __name__ == "__main__":
    compare_distributions()