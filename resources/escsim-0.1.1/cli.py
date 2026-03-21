#!/usr/bin/env python

import sys
import os
import timeit
import re

import click
from importlib.metadata import version as get_version
from core import run_external
from core import create_seeds
from core import calc_phi
from core import calc_popsize_sc
from core import estimate_coaldens_from_popsize
from core import calc_popsize_esc
import hashlib
import datetime

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns


def format_sci(value):
    """Format a float in scientific notation with three significant figures."""
    if value == 0:
        return "0"
    # Use regular notation for values >= 0.001 and < 100000
    if abs(value) >= 0.001 and abs(value) < 100000:
        return f"{value:.3g}"
    exponent = int(np.floor(np.log10(abs(value))))
    mantissa = value / (10 ** exponent)
    # Round mantissa to get 3 significant figures total
    mantissa_rounded = round(mantissa, 2)
    # Remove trailing zeros and unnecessary decimal point
    mantissa_str = f"{mantissa_rounded:.2f}".rstrip('0').rstrip('.')
    if mantissa_str == "0":
        return "0"
    return f"{mantissa_str} \\cdot 10^{{{exponent}}}"


def print_version(ctx, param, value):
    """Callback to print version and exit."""
    if not value or ctx.resilient_parsing:
        return
    try:
        ver = get_version("escsim")
    except Exception:
        ver = "unknown"
    click.echo(f"escsim version {ver}")
    ctx.exit()


@click.group()
@click.option('--version', '-v', is_flag=True, callback=print_version, 
              expose_value=False, is_eager=True, help='Show version and exit.')
def main():
    """escsim - Evolution simulation command line tool."""
    pass


@main.command()
@click.argument('popsize', type=float)
@click.argument('selcoef', type=float)
@click.argument('mutrate', type=float)
@click.argument('sigma', type=float, required=False)
@click.argument('chrmlen', type=int, required=False)
@click.argument('burnin', type=int, required=False)
@click.argument('gens', type=int, required=False)
@click.option('--jobs', '-j', default=1, help='Number of parallel jobs to run')
@click.option('--workers', '-w', default=1, help='Maximum number of worker processes')
@click.option('--folder', '-f', default='tmp/results', help='Output folder for simulation results; will use existing simulations if found (default: tmp/results)')
@click.option('--mode', '-m', default='f', help='Determines whether to run forward simulations with a fixed selection coefficient (f), ' \
'a normally distributed selection coefficient with mean and standard deviation (n) or one drawn out of discrete bins (d)')
def run(popsize, selcoef, sigma, mutrate, chrmlen, burnin, gens, jobs, workers, folder, mode):
    """Run evolution simulations.
    
    Arguments:
        POPSIZE: Population size
        SELCOEF: Selection coefficient
        MUTRATE: Deleterious mutation rate (per chromosome per generation)
        CHRMLEN: Simulated chromosome size in base pairs (optional, default: 15,000)
        BURNIN: Number of burn-in generations (optional, default: 1N)
        GENS: Number of generations to simulate after burn-in (optional, default: 1N)
    """
    click.echo(f"[INFO] Running {jobs} job(s) with {workers} worker(s)...")
    click.echo(f"[INFO] Output folder: {folder}")

    # Add your simulation logic here
    start_time = timeit.default_timer()

    # Ensure output folder exists
    if os.path.exists(folder):
        click.echo(f"[INFO] Output folder already exists.")
    else:
        os.makedirs(folder, exist_ok=True)
        click.echo(f"[INFO] Created output folder.")

    # Set default values for optional parameters
    if chrmlen is None:
        chrmlen = 15000
    if burnin is None:
        burnin = int(popsize)
    if gens is None:
        gens = int(popsize)+burnin
    if sigma is None:
        sigma = 0.0001
    

    # Run simulations
    seeds = create_seeds(n=workers, base_value=str(timeit.default_timer()))
    params = dict(
        popsize=popsize,
        selcoef=selcoef,
        mutrate=mutrate,
        sigma=sigma,
        chrmlen=chrmlen,
        burnin=burnin,
        gens=gens,
        jobs=jobs,
        mode=mode
    )
    results = run_external(seeds=seeds, **params)

    # Test seeds
    seeds_observed = [res["seed"] for res in results]
    if set(seeds) != set(seeds_observed):
        click.echo(f"[WARNING] Mismatch in seeds! Expected: {seeds}, Observed: {seeds_observed}")


    # Write to output folder (Folder name should be escsim_YYYY-MM-DD_rndomstr.out)
    # date = datetime.datetime.now().strftime('%Y-%m-%d')
    # rndstr = hashlib.md5(''.join(map(str, seeds)).encode()).hexdigest()[:8]
    
    escsim_file = os.path.join(folder, f"escsim_N{int(popsize)}_U{mutrate}_s{selcoef}_sigma{sigma}.out")
    wave_file = os.path.join(folder, f"wave_N{int(popsize)}_U{mutrate}_s{selcoef}_sigma{sigma}.out")

    
    with open(escsim_file, 'w') as f_esc, open(wave_file, 'w') as f_wave:
        # Write header
        header_esc = ["sim_id", "seed", "popsize", "selcoef", "mutrate", "velocity", "profile"]
        header_wave = ["time", "wave"]
        
        f_esc.write("\t".join(header_esc) + "\n")
        f_wave.write("\t".join(header_wave) + "\n")
        
    for res in run_external(seeds=seeds, **params):
        # 1. Write metadata to escsim
        profile_str = ",".join(map(str, res["profile"]))
        f_esc.write(f"{res['sim_id']}\t{res['seed']}\t{res['popsize']}\t{profile_str}\n")
        
        # 2. Stream wave data line-by-line
        # This prevents storing the entire wave history in a string before writing
        for t, wave_row in zip(res["time"], res["wave"]):
            wave_str = ",".join(map(str, wave_row))
            f_wave.write(f"{t}\t{wave_str}\n")
        
        # 3. Explicitly clear local reference to large objects
        del res


    # Print summary of parameters
    if mode == "n":
        click.echo(f"[INFO] Parameters: popsize={popsize}, selcoef={selcoef}, mutrate={mutrate}, sigma={sigma}")
    else:
        click.echo(f"[INFO] Parameters: popsize={popsize}, selcoef={selcoef}, mutrate={mutrate}")
    click.echo(f"[INFO] chrmlen={chrmlen}, burnin={burnin}, gens={gens}")
    click.echo(f"[INFO] Number of simulations run: {len(results)}")
    click.echo(f"[INFO] Simulations were run in mode {mode}")
    click.echo(f"[INFO] Results written to: {os.path.basename(escsim_file)}")
    click.echo(f"[INFO] Wave file written to: {os.path.basename(wave_file)}")
    

    # Final message
    end_time = timeit.default_timer()
    elapsed = end_time - start_time
    hours = int(elapsed // 3600)
    minutes = int((elapsed % 3600) // 60)
    seconds = elapsed % 60
    click.echo(f"[INFO] Total run time: {hours:02d}:{minutes:02d}:{seconds:05.2f}")
    click.echo("[INFO] Simulation complete!")


@main.command()
@click.argument('figure_pdf', type=click.Path())
@click.option('--input_folder', '-i', default='tmp/results', help='Folder containing simulation results')
@click.option('--output', '-o', default='results/escsim_figures', help='Output folder for figures')
@click.option('--no-sep-sumplot', '-n', is_flag=True, help='Do not create separate summary plot for mean velocities')
def summarize(figure_pdf, input_folder, output, no_sep_sumplot):
    """Plot and summarize simulation results."""
    matplotlib.use('agg')  # Use non-interactive backend
    click.echo(f"[INFO] Plotting results from folder: {input_folder}")

    # Get all files that have the form escsim_YYYY-MM-DD_*.out
    # sim_files = [f for f in os.listdir(input_folder) if re.match(r"escsim_\d{4}-\d{2}-\d{2}_.+\.out", f)]
    # Get all files that have the form escsim_N{N}_U{U}_s{s}.out
    num = r"\d+(?:\.\d+)?(?:e-?\d+)?"
    sim_files = [f for f in os.listdir(input_folder)
                 if re.match(rf"escsim_N\d+_U{num}_s{num}_sigma{num}\.out",
                 f)]
    click.echo(f"[INFO] Found {len(sim_files)} simulation result file(s).")
    df_list = []
    for _, sim_file in enumerate(sim_files):
        click.echo(f"[INFO] Processing file: {sim_file} ({_ + 1}/{len(sim_files)})")
        filepath = os.path.join(input_folder, sim_file)
        df = pd.read_csv(filepath, sep="\t")
        df_list.append(df)
    df = pd.concat(df_list, ignore_index=True).drop("sim_id", axis=1)

    # Assert seeds are unique
    if df['seed'].nunique() != len(df):
        click.echo(f"[WARNING] Duplicate seeds found in the combined data!")


    # Begin plotting
    figures = []

    ## Summarize all mean velocities into one figure
    sns.set(style="ticks", context="paper")
    fig, ax = plt.subplots(figsize=(8, 6))

    velocity_summary = df.groupby(["popsize", "selcoef", "mutrate"])['velocity'].agg(['mean', 'std', 'count']).reset_index()
    velocity_summary["Ns"] = velocity_summary['popsize'] * velocity_summary['selcoef']
    velocity_summary["lineid"] = velocity_summary.apply(lambda row: f"$N={format_sci(row['popsize'])}$, $U_d={format_sci(row['mutrate'])}$", axis=1)
    velocity_summary["Popsize $N$"] = velocity_summary['popsize'].apply(lambda x: format_sci(x))
    velocity_summary["Mutation rate $U_d$"] = velocity_summary['mutrate'].apply(lambda x: format_sci(x))


    # As the x scale will be log, make the zero to be on the xlimits as if it weren't zero
    velocity_summary_complete = velocity_summary.copy()
    velocity_summary = velocity_summary.loc[velocity_summary['selcoef'] != 0]

    ax.axhline(0, color='darkgray', linestyle='--')
    ax.axhline(1, color='darkgray', linestyle='--')

    print(velocity_summary.head())

    sns.lineplot(
        data=velocity_summary,
        x='Ns',
        y='mean',
        hue='lineid',
        ax=ax,
        palette='tab10',
        legend=False,
    )
    ax.set_xscale('log')

    # Get x limit and save them for overplotting
    xlim_left, xlim_right = ax.get_xlim()
    
    # Modify the zero selcoef points to be at xlim_left
    velocity_summary_complete.loc[velocity_summary_complete['selcoef'] == 0, 'Ns'] = xlim_left
    velocity_summary = velocity_summary_complete

    sns.lineplot(
        data=velocity_summary,
        x='Ns',
        y='mean',
        hue='lineid',
        ax=ax,
        palette='tab10',
        linestyle='dashed',
        legend=False,
        )
    sns.scatterplot(
        data=velocity_summary,
        x='Ns',
        y='mean',
        hue='lineid',
        style='Popsize $N$',
        ax=ax,
        palette='tab10',
        legend=True,
        s=50,
    )

    ax.set_yticks(np.arange(-.1, 1.1, 0.1), minor=False)
    ax.set_yticks(np.arange(-.1, 1.1, 0.05), minor=True)
    ax.set_ylim(-0.05, 1.05)

    ax.set_xscale('log')
    ax.set_xlim(xlim_left, xlim_right)

    ax.set_xlabel("$N \\cdot s$")
    ax.set_ylabel("Estimated mean velocity")
    ax.set_title("Mean velocity across parameter combinations")
    ax.legend(bbox_to_anchor=(0.975, 1), loc='upper right', framealpha=1.0)
    figures.append(fig)
    plt.close()

    # Create name of separate summary plot if needed
    if not no_sep_sumplot:
        base, ext = os.path.splitext(figure_pdf)
        sumplot_name = f"{base}_mean_velocity_summary{ext}"

        # Creating output directory
        os.makedirs(output, exist_ok=True)
        click.echo(f"[INFO] Saving separate mean velocity summary plot to: {output}/{sumplot_name}")
        with PdfPages(f"{output}/{sumplot_name}.pdf") as pdf:
            pdf.savefig(fig)
            plt.close(fig)


    ## Figure logic of individual parameter combinations
    # Loop through the unique parameter combinations
    param_cols = ['popsize', 'selcoef', 'mutrate']
    df_sorted = df.sort_values(by=['popsize', 'mutrate', 'selcoef'], ascending=[False, False, False])
    grouped = df_sorted.groupby(param_cols, sort=False)
    for params, group in grouped:
        popsize, selcoef, mutrate = params
        click.echo(f"[INFO] Plotting for parameters: popsize={popsize}, selcoef={selcoef}, mutrate={mutrate}")

        # Calculate phi
        phi = calc_phi(popsize, selcoef, mutrate)

        # Create a figure for this parameter set
        sns.set(style="ticks", context="paper")
        fig, axs = plt.subplots(2, 2, figsize=(8, 6))
        ax = axs[0, 0]

        ## Plot A
        # Plot velocity distribution as a histogram
        sns.histplot(group['velocity'], bins=np.arange(-0.15, 1.15, 0.02), kde=False, ax=ax, color='gray', edgecolor='none')
        ax.set_xlabel("Velocity")
        ax.set_ylabel("Count")

        ax.axvline(np.mean(group['velocity']), color='red', linestyle='--', label='Mean Velocity')
        
        ax.set_xlim(-0.135, 1.135)
        ax.set_xticks([0.0, 0.5, 1.0], minor=False)
        ax.set_xticks(np.arange(-0.1, 1.1, 0.1), minor=True)

        ax.set_title(f"Velocity Distribution")

        
        ## Plot B
        # Plot mutational burden profile
        ax = axs[0, 1]
        profiles = group['profile'].apply(lambda x: np.fromstring(x, sep=','))
        profiles_len = max([len(p) for p in profiles])
        profile_matrix = np.zeros((len(profiles), profiles_len))
        for i, p in enumerate(profiles):
            profile_matrix[i, :len(p)] = p
        mean_profile = np.mean(profile_matrix, axis=0)
        ci_lower = np.percentile(profile_matrix, 2.5, axis=0)
        ci_upper = np.percentile(profile_matrix, 97.5, axis=0)
        click.echo(f"[INFO] Original mutational burden profile length: {len(mean_profile)}")
        threshold = 0.001 * np.max(mean_profile)
        last_idx = np.where(mean_profile >= threshold)[0][-1] + 1
        mean_profile = mean_profile[:last_idx]
        mean_profile /= np.sum(mean_profile)  # Normalize
        ci_lower = ci_lower[:last_idx]
        ci_upper = ci_upper[:last_idx]
        click.echo(f"[INFO] Mutational burden profile length after thresholding: {len(mean_profile)}")
        ax.plot(range(len(mean_profile)), mean_profile, color='gray')
        ax.fill_between(range(len(mean_profile)), ci_lower, ci_upper, color='gray', alpha=0.3)
        ax.scatter(range(len(mean_profile)), mean_profile, color='black', marker="s")

        ax.set_xlim(0, None)
        ax.set_ylim(0, None)

        # Set major ticks (at most 10 labels)
        step = max(1, len(mean_profile) // 10)
        ax.set_xticks(np.arange(0, len(mean_profile), step=step), minor=False)
        # Set minor ticks at all integers
        ax.set_xticks(np.arange(0, len(mean_profile), step=1), minor=True)

        ax.set_xlabel("Mutational burden class")
        ax.set_ylabel("Density")

        ax.set_title("Mutational burden profile")


        ## Plot C
        ax = axs[1, 0]
        # Coalescent density of SC (Nicolaisen and Desai 2012) vs ESC (Strütt et al. 2024)
        # First popsize over time for SC
        params = dict(popsize=popsize, selcoef=selcoef, mutrate=mutrate, tmin=0, tmax=2*popsize, ntimes=300)
        sc_times, sc_popsize = calc_popsize_sc(**params)
        sc_times, sc_density = estimate_coaldens_from_popsize(sc_times, sc_popsize, nsam=2)
        ax.plot(sc_times, sc_density, label='SC', color='pink', linewidth=4)

        yticks = ax.get_yticks()
        ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(yticks))
        ax.set_yticklabels([f"${format_sci(ytick)}$" for ytick in yticks])

        # Second popsize over time for ESC
        params = dict(times=sc_times, popsize=popsize, mutrate=mutrate, velocity=np.mean(group['velocity']), profile=mean_profile)
        esc_times, esc_popsize = calc_popsize_esc(**params)
        _, esc_density = estimate_coaldens_from_popsize(esc_times, esc_popsize, nsam=2)
        ax.plot(esc_times, esc_density, label='ESC', color='darkgreen')

        ax.set_xlim(0, 1.15*popsize)
        ax.set_ylim(0, 1.15*np.max(esc_density)*1.1)
        ax.set_xlabel("Time ago (gen)")
        ax.set_ylabel("Coalescent density")

        ax.set_title("Coalescent density predictions")


        ## Plot D
        # Population sizes of SC (Nicolaisen and Desai 2012) vs ESC (Strütt et al. 2024)
        ax = axs[1, 1]
        ax.plot(sc_times, sc_popsize, label='SC', color='pink', linewidth=4)
        ax.plot(esc_times, esc_popsize, label='ESC', color='darkgreen')

        ax.set_xlim(0, 1.15*popsize)
        ax.set_ylim(0, 1.15*popsize)
        ax.set_xlabel("Time ago (gen)")
        ax.set_ylabel("Population size")

        ax.legend()
        ax.set_title("Population size predictions")

        ## Figure adjustments
        # Add labels for subplots
        for i, axi in enumerate(axs.flat):
            axi.text(axi.get_xlim()[0]*0.9, axi.get_ylim()[1]*1.05, chr(65 + i), fontsize=12)


        fig.suptitle(
            f"\n$N={format_sci(popsize)}$, $s={format_sci(selcoef)}$, $U_d={format_sci(mutrate)}$, $\\phi={format_sci(phi)}$"+
            f"\nMean Velocity={np.mean(group['velocity']):.4f} ± {np.std(group['velocity']):.4f}"+
            f", $n={len(group)}$ simulations"
        )
        
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        figures.append(fig)
        plt.close()

    # Creatting output directory if not done prior
    os.makedirs(output, exist_ok=True)
    ## Plotting logic to make one pdf page per figure
    click.echo(f"[INFO] Saving figures to PDF: {output}/{figure_pdf}")
    with PdfPages(f"{output}/{figure_pdf}.pdf") as pdf:
        for figid, fig in enumerate(figures):
            click.echo(f"[INFO] Saving figure {figid + 1}/{len(figures)} to PDF.")
            pdf.savefig(fig)
            plt.close(fig)


if __name__ == "__main__":
    main()
