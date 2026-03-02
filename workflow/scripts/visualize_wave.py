from pathlib import Path
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
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
@click.argument('pop_size', type=int)
@click.argument('sel_coef', type=float)
@click.argument('mut_rate', type=float)

@click.option('--input_folder', '-i', default='tmp/results',
            help='Input folder for simulation results (default: tmp/results)')
@click.option('--output', '-o', default='results/escsim_figures',
              help='Output folder for coalescent densities '
              '(default: results/coalescent_densities)')

def visualize_wave(input_folder, output, pop_size, mut_rate, sel_coef):
    
    file_path = Path(input_folder)
    if not file_path.exists():
        click.echo(f"Error: File {file_path} does not exist.")
        return
    
    # Ensure output folder exists
    if os.path.exists(output):
        click.echo("[INFO] Output folder already exists.")
    else:
        os.makedirs(output, exist_ok=True)
        click.echo("[INFO] Created output folder.")

    times, waves = load_wave_file(f"{input_folder}/wave_N{pop_size}_U{mut_rate}_s{sel_coef}.out")
    wave_matrix = np.vstack(waves)

    density = compute_density_matrix(wave_matrix)

    plt.figure(figsize=(10, 6))

    # Use LogNorm to map 1..max to colors, zeros will be white
    cmap = plt.cm.inferno 
    cmap.set_under("white")  # values below vmin are white

    plt.imshow(
        density.T,               # transpose so mutational load is vertical
        aspect="auto",
        origin="lower",
        extent=[times.min(), times.max(), 0, density.shape[1]],
        norm=LogNorm(vmin=1, vmax=density.max()),
        cmap=cmap
        )

    plt.colorbar(label="Number of individuals (log scale)")
    plt.xlabel("Time")
    plt.ylabel("Mutational load")
    plt.title(f"N = {pop_size} U = {mut_rate} s = {sel_coef}")
    plt.tight_layout()

    plt.savefig(f"{output}/wave_N{pop_size}_U{mut_rate}_s{sel_coef}.pdf")
    click.echo(f"[INFO] Saved plot to {output}.")


if __name__ == "__main__":
    visualize_wave()