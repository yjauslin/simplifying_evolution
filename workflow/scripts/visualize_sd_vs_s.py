import click
import numpy as np
import pandas as pd
import os
import re
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
from matplotlib.ticker import ScalarFormatter

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

    # Enforce mutual exclusivity
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

    plt.rcParams.update({
        "text.usetex": True,
        "mathtext.fontset": "cm",        
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
    })

    # Base figure setup (matching LaTeX column width)
    fig_width = 426.79134 / 72.27  
    fig_height = fig_width / 1.618

    # Individual figures
    fig1, ax1 = plt.subplots(figsize=(fig_width, fig_height))
    fig2, ax2 = plt.subplots(figsize=(fig_width, fig_height))

    # Combined Figure Setup (2 side-by-side subplots fitting in the same figure dimensions)
    fig3, (ax3_left, ax3_right) = plt.subplots(1, 2, figsize=(fig_width, fig_height))

    markers = ['o', '^', 's']
    is_sel_coef_mode = len(sel_coef) > len(mut_rate) or (len(sel_coef) == len(mut_rate) == 1)

    # Fontsize configuration
    LABEL_FONTSIZE = 11
    LEGEND_FONTSIZE = 9
    TITLE_FONTSIZE = 10
    TICK_FONTSIZE = 9

    if is_sel_coef_mode:
        # Sort ascending for Selection Coefficients
        sorted_sel_coef = sorted(sel_coef)
        colors = sns.color_palette("colorblind", len(sorted_sel_coef))

        for i, s in enumerate(sorted_sel_coef):
            input_file = os.path.join(input_folder, f"s_N{pop_size}_U{mut_rate[0]}_s{s}.txt")
            df = pd.read_csv(input_file, sep="\t")
            s_s_eff = df['s'] / s

            label = f'${s}$'
            marker = markers[i % len(markers)]
            color = colors[i]

            # Individual figures
            sns.scatterplot(x='sd/s', y='s', data=df, ax=ax1, color=color, marker=marker, s=25, label=label)
            sns.scatterplot(x=df['sd/s'], y=s_s_eff, ax=ax2, color=color, marker=marker, s=25, label=label)

            # Combined figure panels
            sns.scatterplot(x='sd/s', y='s', data=df, ax=ax3_left, color=color, marker=marker, s=20, label=label)
            sns.scatterplot(x=df['sd/s'], y=s_s_eff, ax=ax3_right, color=color, marker=marker, s=20, label=label)

        # Formatting Individual Fig 1
        ax1.set_ylim(-0.05 * min(sorted_sel_coef), max(sorted_sel_coef) * 1.15)
        ax1.set_xlabel(r"$\sigma/s_{normal}$", fontsize=LABEL_FONTSIZE, labelpad=6)
        ax1.set_ylabel("Effective Selection Coefficient", fontsize=LABEL_FONTSIZE, labelpad=6)
        ax1.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
        ax1.legend(title=r"$s_{normal}$", title_fontsize=TITLE_FONTSIZE, loc="upper right", ncol=1, frameon=False, fontsize=LEGEND_FONTSIZE)
        fig1.tight_layout(pad=0.1)
        fig1.savefig(os.path.join(output, "effective_selection_coefficient.jpg"), dpi=600, bbox_inches='tight')

        # Formatting Individual Fig 2
        ax2.set_ylim(0, 1.05)
        ax2.set_xlabel(r"$\sigma/s_{normal}$", fontsize=LABEL_FONTSIZE, labelpad=6)
        ax2.set_ylabel(r"Effective Selection Coefficient / $s_{normal}$", fontsize=LABEL_FONTSIZE, labelpad=6)
        ax2.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
        ax2.legend(title=r"$s_{normal}$", title_fontsize=TITLE_FONTSIZE, fontsize=LEGEND_FONTSIZE, loc="upper right", frameon=False)
        fig2.tight_layout(pad=0.1)
        fig2.savefig(os.path.join(output, "relative_effective_selection_coefficient.jpg"), dpi=600, bbox_inches='tight')

        # Formatting Combined Fig 3 (Selection Coefficient Mode)
        ax3_left.set_ylim(-0.05 * min(sorted_sel_coef), max(sorted_sel_coef) * 1.15)
        ax3_left.set_xlabel(r"$\sigma/s_{normal}$", fontsize=LABEL_FONTSIZE, labelpad=4)
        ax3_left.set_ylabel("Effective Selection Coefficient", fontsize=LABEL_FONTSIZE, labelpad=4)
        ax3_left.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
        ax3_left.legend(title=r"$s_{normal}$", title_fontsize=TITLE_FONTSIZE, loc="upper right", ncol=1, frameon=False, fontsize=LEGEND_FONTSIZE)

        ax3_right.set_ylim(0, 1.05)
        ax3_right.set_xlabel(r"$\sigma/s_{normal}$", fontsize=LABEL_FONTSIZE, labelpad=4)
        ax3_right.set_ylabel(r"Effective Selection Coefficient / $s_{normal}$", fontsize=LABEL_FONTSIZE, labelpad=4)
        ax3_right.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
        ax3_right.legend(title=r"$s_{normal}$", title_fontsize=TITLE_FONTSIZE, loc="upper right", frameon=False, fontsize=LEGEND_FONTSIZE)

    else:
        # Sort descending for Mutation Rates
        sorted_mut_rate = sorted(mut_rate, reverse=True)
        colors = sns.color_palette("colorblind", len(sorted_mut_rate))
        
        lowest_u = sorted_mut_rate[-1]
        df_lowest_u = None

        for i, u in enumerate(sorted_mut_rate):
            input_file = os.path.join(input_folder, f"U_N{pop_size}_U{u}_s{sel_coef[0]}.txt")
            df = pd.read_csv(input_file, sep="\t")
            u_eff_u = df['U'] / u

            if u == lowest_u:
                df_lowest_u = df

            label = f'${u}$'
            marker = markers[i % len(markers)]
            color = colors[i]

            # Individual figures
            sns.scatterplot(x='sd/s', y='U', data=df, ax=ax1, label=label, color=color, marker=marker, s=25)
            sns.scatterplot(x=df['sd/s'], y=u_eff_u, ax=ax2, color=color, marker=marker, s=25, label=label)

            # Combined figure panels
            sns.scatterplot(x='sd/s', y='U', data=df, ax=ax3_left, color=color, marker=marker, s=20, label=label)
            sns.scatterplot(x=df['sd/s'], y=u_eff_u, ax=ax3_right, color=color, marker=marker, s=20, label=label)

        # Scientific Notation Formatter for Insets
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((0, 0))

        # --- Inset for Individual Fig 1 ---
        ax1_inset = ax1.inset_axes([0.55, 0.58, 0.40, 0.36])
        sns.scatterplot(x='sd/s', y='U', data=df_lowest_u, ax=ax1_inset, color=colors[-1], marker=markers[(len(sorted_mut_rate)-1) % len(markers)], s=15)
        ax1_inset.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE - 2)
        ax1_inset.yaxis.set_major_formatter(formatter)
        ax1_inset.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax1_inset.yaxis.get_offset_text().set_fontsize(TICK_FONTSIZE - 3)
        ax1_inset.set_xlabel(r"$\sigma/s_{normal}$", fontsize=LABEL_FONTSIZE - 2, labelpad=1)
        ax1_inset.set_ylabel(r"Effective $U$", fontsize=LABEL_FONTSIZE - 2, labelpad=1)

        # --- Inset for Combined Fig 3 ---
        ax3_inset = ax3_left.inset_axes([0.55, 0.52, 0.42, 0.38])
        sns.scatterplot(x='sd/s', y='U', data=df_lowest_u, ax=ax3_inset, color=colors[-1], marker=markers[(len(sorted_mut_rate)-1) % len(markers)], s=12)
        ax3_inset.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE - 3)
        ax3_inset.yaxis.set_major_formatter(formatter)
        ax3_inset.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax3_inset.yaxis.get_offset_text().set_fontsize(TICK_FONTSIZE - 4)
        ax3_inset.set_xlabel(r"$\sigma/s_{normal}$", fontsize=LABEL_FONTSIZE - 3, labelpad=1)
        ax3_inset.set_ylabel(r"Effective Mutation Rate", fontsize=LABEL_FONTSIZE - 3, labelpad=1)

        # Formatting Individual Fig 1
        ax1.set_ylim(-0.001, max(sorted_mut_rate) * 1.25)
        ax1.set_xlabel(r"$\sigma/s_{normal}$", fontsize=LABEL_FONTSIZE, labelpad=6)
        ax1.set_ylabel("Effective Mutation Rate", fontsize=LABEL_FONTSIZE, labelpad=6)
        ax1.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
        ax1.legend(title=r"$U_{normal}$", title_fontsize=TITLE_FONTSIZE, fontsize=LEGEND_FONTSIZE, loc="lower right", bbox_to_anchor=(0.98, 0.10), frameon=False)
        fig1.tight_layout(pad=0.1)
        fig1.savefig(os.path.join(output, "effective_mutation_rate.jpg"), dpi=600, bbox_inches='tight')

        # Formatting Individual Fig 2
        ax2.set_ylim(-0.05, 1.40)
        ax2.set_xlabel(r"$\sigma/s_{normal}$", fontsize=LABEL_FONTSIZE, labelpad=6)
        ax2.set_ylabel(r"Effective Mutation Rate / $U_{normal}$", fontsize=LABEL_FONTSIZE, labelpad=6)
        ax2.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
        ax2.legend(title=r"$U_{normal}$", title_fontsize=TITLE_FONTSIZE, loc="upper right", frameon=False, fontsize=LEGEND_FONTSIZE)
        fig2.tight_layout(pad=0.1)
        fig2.savefig(os.path.join(output, "relative_effective_mutation_rate.jpg"), dpi=600, bbox_inches='tight')

        # Formatting Combined Fig 3 (Effective Mutation Rate Mode)
        ax3_left.set_ylim(-0.001, max(sorted_mut_rate) * 1.30)
        ax3_left.set_xlabel(r"$\sigma/s_{normal}$", fontsize=LABEL_FONTSIZE, labelpad=4)
        ax3_left.set_ylabel("Effective Mutation Rate", fontsize=LABEL_FONTSIZE, labelpad=4)
        ax3_left.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
        # Position legend cleanly above bottom baseline points
        ax3_left.legend(title=r"$U_{normal}$", title_fontsize=TITLE_FONTSIZE, loc="lower right", bbox_to_anchor=(0.98, 0.12), ncol=1, frameon=False, fontsize=LEGEND_FONTSIZE)

        # Extend y-limit on right panel to create top-right clearance above peak data (~1.25)
        ax3_right.set_ylim(-0.05, 1.40)
        ax3_right.set_xlabel(r"$\sigma/s_{normal}$", fontsize=LABEL_FONTSIZE, labelpad=4)
        ax3_right.set_ylabel(r"Effective Mutation Rate / $U_{normal}$", fontsize=LABEL_FONTSIZE, labelpad=4)
        ax3_right.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)
        ax3_right.legend(title=r"$U_{normal}$", title_fontsize=TITLE_FONTSIZE, loc="upper right", bbox_to_anchor=(0.98, 0.98), frameon=False, fontsize=LEGEND_FONTSIZE)

    # Panel labels 'a' and 'b' accurately aligned above their respective subplots
    ax3_left.text(-0.20, 1.05, 'a', transform=ax3_left.transAxes, fontsize=12, fontweight='bold', va='top', ha='right')
    ax3_right.text(-0.15, 1.05, 'b', transform=ax3_right.transAxes, fontsize=12, fontweight='bold', va='top', ha='right')

    # Save combined figure cleanly without excess whitespace
    fig3.tight_layout(pad=0.2)
    if is_sel_coef_mode:
        fig3.savefig(os.path.join(output, "combined_effective_selection_coefficient.jpg"), dpi=600, bbox_inches='tight')
    else:
        fig3.savefig(os.path.join(output, "combined_effective_mutation_rate.jpg"), dpi=600, bbox_inches='tight')


if __name__ == "__main__":
    visualize_sd_vs_s()