import click
import numpy as np
import pandas as pd
import os
import re
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress

@click.command()
@click.argument('pop_size', required=False, type=int, default=5000)
@click.option('--sel_coef', '-s', required=False, default=(0.001,), type=float, multiple=True, help='Selection coefficient')
@click.option('--mut_rate', '-u', required=False, default=(0.006,), type=float, multiple=True, help='Mutation rate')
@click.option('--input_folder', '-i', default='tmp/results', required=False, type=str,
              help='Input folder containing comparison files.')
@click.option('--output', '-o', default='tmp/results', required=False,
              help='Output folder')
def visualize_sd_vs_s(pop_size, sel_coef, mut_rate, input_folder, output):
    # Create output directory if it doesn't exist
    os.makedirs(output, exist_ok=True)

    sel_coef = list(sel_coef)
    mut_rate = list(mut_rate)

    # enforce mutual exclusivity, meaning: you cannot provide multiple selection coefficients and multiple mutation rates at the same time
    if len(sel_coef) > 1 and len(mut_rate) > 1:
        raise click.UsageError(
            "You cannot provide multiple --sel_coef AND multiple --mut_rate at the same time."
        )

    if len(sel_coef) == 0 and len(mut_rate) == 0:
        raise click.UsageError(
            "You must provide at least one of --sel_coef or --mut_rate."
        )
    
    sns.set_context("paper")
    sns.set_style("ticks")

    # fig width corresponds to column with in LateX, 72.27 is the conversion factor from points to inches
    fig_width = 426.79134 / 72.27  
    fig_height = fig_width / 1.618

    fig1, ax1 = plt.subplots(figsize=(fig_width, fig_height))
    fig2, ax2 = plt.subplots(figsize=(fig_width, fig_height))

    # Set color palette based on the number of selection coefficients or mutation rates
    colors = sns.color_palette("colorblind", len(sel_coef) if len(sel_coef) > 1 else len(mut_rate))
    markers = ['o', '^', 's']

    # if more selection coefficients than mutation rates are provided, we will plot effective selection coefficient vs sd/s, otherwise we will plot effective mutation rate vs sd/s
    if len(sel_coef) > len(mut_rate) or (len(sel_coef) == len(mut_rate) == 1):
        # for every provided selection coefficient...
        for i in range(len(sel_coef)):
            s = sel_coef[i]
            # Construct the input file path based on the provided parameters
            input_file = os.path.join(
                input_folder,
                f"s_N{pop_size}_U{mut_rate[0]}_s{s}.txt")
            df = pd.read_csv(input_file, sep="\t")

            s_s_eff = df['s'] / s

            # Perform linear regression on the data to find the trend line
            slope, intercept, r_value, p_value, std_err = linregress(df['sd/s'], df['s'])
            # Calculate the coefficient of determination (R-squared) to assess the goodness of fit
            r_squared = r_value**2

            # Generate a range of x values for plotting the trend line
            x_range = np.linspace(df['sd/s'].min(), df['sd/s'].max(), 100)
            # Calculate the corresponding y values for the trend line using the linear regression parameters
            y_trend = slope * x_range + intercept

            sns.scatterplot(x='sd/s', y='s', data=df, ax=ax1, color=colors[i], marker=markers[i], s=25, label=f'$s={s} (R^2={r_squared:.3f})$')
            # Plot the trend line with a dashed black line and include the R-squared value in the label
            ax1.plot(x_range, y_trend, color='black', linestyle='--', linewidth=1)

            sns.scatterplot(x=df['sd/s'], y=s_s_eff, ax=ax2, color=colors[i], marker=markers[i], s=25, label=f'$s={s}$')
        ax1.set_ylim(0, max(sel_coef) * 1.1)
        ax1.set_xlabel(r"$\sigma/s$", fontsize=11, labelpad=6)
        ax1.set_ylabel("effective selection coefficient", fontsize=11, labelpad=6)
        ax1.tick_params(axis='both', which='major', labelsize=9)
        ax1.legend(title="Selection coefficient", title_fontsize=10, fontsize=9, loc="best", frameon=False)
        fig1.tight_layout(pad=0.1)
        fig1.savefig(os.path.join(output, f"effective_selection_coefficient.jpg"), dpi=600, bbox_inches='tight')

        ax2.set_ylim(0, 1)
        ax2.set_xlabel(r"$\sigma/s$", fontsize=11, labelpad=6)
        ax2.set_ylabel("effective selection coefficient / s", fontsize=11, labelpad=6)
        ax2.tick_params(axis='both', which='major', labelsize=9)
        ax2.legend(title="Selection coefficient", title_fontsize=10, fontsize=9, loc="best", frameon=False)
        fig2.tight_layout(pad=0.1)
        fig2.savefig(os.path.join(output, f"relative_effective_selection_coefficient.jpg"), dpi=600, bbox_inches='tight')

    else:
        # for every provided mutation rate...
        for i in range(len(mut_rate)):
            u = mut_rate[i]
            input_file = os.path.join(
                input_folder,
                f"U_N{pop_size}_U{u}_s{sel_coef[0]}.txt")
            df = pd.read_csv(input_file, sep="\t")

            u_eff_u = df['U'] / u

            # Perform linear regression on the data to find the trend line
            slope, intercept, r_value, p_value, std_err = linregress(df['sd/s'], df['U'])
            # Calculate the coefficient of determination (R-squared) to assess the goodness of fit
            r_squared = r_value**2

            # Generate a range of x values for plotting the trend line
            x_range = np.linspace(df['sd/s'].min(), df['sd/s'].max(), 100)
            # Calculate the corresponding y values for the trend line using the linear regression parameters
            y_trend = slope * x_range + intercept

            sns.scatterplot(x='sd/s', y='U', data=df, ax=ax1, label=f'$U_d={u} (R^2={r_squared:.3f})$', color=colors[i], marker=markers[i], s=25)
            # Plot the trend line with a dashed black line and include the R-squared value in the label
            ax1.plot(x_range, y_trend, color='black', linestyle='--')
            sns.scatterplot(x=df['sd/s'], y=u_eff_u, ax=ax2, color=colors[i], marker=markers[i], s=25, label=f'$U_d={u}$')


        ax1.set_ylim(0, max(mut_rate) * 1.1)
        ax1.set_xlabel(r"$\sigma/s$", fontsize=11, labelpad=6)
        ax1.set_ylabel("effective mutation rate", fontsize=11, labelpad=6)
        ax1.tick_params(axis='both', which='major', labelsize=9)
        ax1.legend(title="Mutation rate", title_fontsize=10, fontsize=9, loc="best", frameon=False)
        fig1.tight_layout(pad=0.1)
        fig1.savefig(os.path.join(output, f"effective_mutation_rate.jpg"), dpi=600, bbox_inches='tight')

        ax2.set_ylim(0, 1)
        ax2.set_xlabel(r"$\sigma/s$", fontsize=11, labelpad=6)
        ax2.set_ylabel("effective mutation rate / $U_d$", fontsize=11, labelpad=6)
        ax2.tick_params(axis='both', which='major', labelsize=9)
        ax2.legend(title="Mutation rate", title_fontsize=10, fontsize=9, loc="best", frameon=False)
        fig2.tight_layout(pad=0.1)
        fig2.savefig(os.path.join(output, f"relative_effective_mutation_rate.jpg"), dpi=600, bbox_inches='tight')

if __name__ == "__main__":
    visualize_sd_vs_s()