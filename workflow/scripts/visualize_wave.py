from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import click


def load_wave_file(filename):
    times = []
    waves = []

    with open(filename) as f:
        next(f)  # skip header
        for line in f:
            t, wave_str = line.rstrip().split("\t")
            times.append(int(t))
            waves.append(np.fromstring(wave_str, sep=",", dtype=int))

    return np.array(times), waves


def compute_density_matrix(wave_matrix):
    max_load = wave_matrix.max()
    T = wave_matrix.shape[0]

    density = np.zeros((T, max_load + 1), dtype=int)

    for i in range(T):
        density[i] = np.bincount(
            wave_matrix[i],
            minlength=max_load + 1
        )

    return density


@click.command()
@click.option('--input_folder', '-i', default='results/escsim',
            help='Input folder for wave results from forward simulations'
            '(default: results/escsim)')
@click.option('--output', '-o', default='results/escsim_figures',
              help='Output folder for wave plots summary file '
              '(default: results/escsim_figures)')

def summarize_waves(input_folder, output):

    input_path = Path(input_folder)
    files = sorted(input_path.glob("wave_*.out"))

    if len(files) == 0:
        click.echo("No wave files found.")
        return

    n = len(files)

    cols = 4
    rows = int(np.ceil(n / cols))

    fig = plt.figure(figsize=(16, 4 * rows))

    # grid with extra column for colorbar
    gs = plt.GridSpec(
        rows,
        cols + 1,
        width_ratios=[1]*cols + [0.05],
        figure=fig
    )

    # create subplot axes
    axes = []
    for r in range(rows):
        for c in range(cols):
            axes.append(fig.add_subplot(gs[r, c]))

    # colorbar axis (spans all rows)
    cax = fig.add_subplot(gs[:, -1])

    cmap = plt.cm.inferno
    cmap.set_under("white")

    vmax_global = 0
    densities = []
    times_list = []

    # first pass: compute densities
    for file in files:
        times, waves = load_wave_file(file)
        wave_matrix = np.vstack(waves)
        density = compute_density_matrix(wave_matrix)

        densities.append(density)
        times_list.append(times)

        vmax_global = max(vmax_global, density.max())

    # second pass: plot
    for ax, file, density, times in zip(axes, files, densities, times_list):

        im = ax.imshow(
            density.T,
            aspect="auto",
            origin="lower",
            extent=[times.min(), times.max(), 0, density.shape[1]],
            norm=LogNorm(vmin=1, vmax=vmax_global),
            cmap=cmap
        )

        ax.set_title(file.stem.replace("wave_", ""))
        ax.set_xlabel("Time")
        ax.set_ylabel("Mutational load")

    # hide unused axes
    for ax in axes[n:]:
        ax.axis("off")

    # add colorbar in dedicated axis
    fig.colorbar(im, cax=cax, label="Individuals (log scale)")

    plt.tight_layout()

    os.makedirs(output, exist_ok=True)
    plt.savefig(f"{output}/wave_summary.pdf", bbox_inches="tight")

    click.echo(f"[INFO] Saved summary plot to {output}")


if __name__ == "__main__":
    summarize_waves()