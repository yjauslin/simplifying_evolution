import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages
import click

def load_multi_wave_file(filename):
    """
    Yields (times, waves) tuples. 
    Detects a new matrix when the 'time' value decreases (reset).
    """
    current_times, current_waves = [], []
    last_t = -1

    with open(filename) as f:
        next(f)  # skip the very first header
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 2: continue
            
            t = int(parts[0])
            wave_data = np.fromstring(parts[1], sep=",", dtype=int)

            # If time resets, yield the completed matrix
            if t < last_t and current_times:
                yield np.array(current_times), current_waves
                current_times, current_waves = [], []
            
            current_times.append(t)
            current_waves.append(wave_data)
            last_t = t
            
        if current_times:
            yield np.array(current_times), current_waves

def compute_density_matrix(wave_matrix):
    max_load = wave_matrix.max()
    T = wave_matrix.shape[0]
    density = np.zeros((T, max_load + 1), dtype=int)
    for i in range(T):
        density[i] = np.bincount(wave_matrix[i], minlength=max_load + 1)
    return density

@click.command()
@click.option('--input_folder', '-i', default='results/escsim')
@click.option('--output', '-o', default='results/escsim_figures')
@click.option('--mode', '-m', default='f')
def summarize_waves(input_folder, output, mode):
    input_path = Path(input_folder)
    files = sorted(input_path.glob("wave_*.out"))

    if not files:
        click.echo("No wave files found.")
        return

    os.makedirs(output, exist_ok=True)
    out_name = "wave_normal_summary.pdf" if mode == 'n' else "wave_fixed_summary.pdf"
    pdf_path = Path(output) / out_name

    cmap = plt.cm.inferno
    cmap.set_under("white")

    with PdfPages(pdf_path) as pdf:
        for file in files:
            # Load all matrices for this specific file
            matrices = list(load_multi_wave_file(file))
            
            # Start first page for this file
            fig = None
            plot_count = 0

            for i, (times, waves) in enumerate(matrices):
                # Start a new page if 12 plots reached OR first matrix of a file
                if plot_count % 12 == 0:
                    if fig:
                        pdf.savefig(fig, bbox_inches="tight")
                        plt.close(fig)
                    
                    fig = plt.figure(figsize=(18, 12))
                    # Title only on the first page of the file, or every page? 
                    # Prompt says: "start a new page for the new file and write a title again"
                    fig.suptitle(f"File: {file.name}", fontsize=16, fontweight='bold')
                    gs = plt.GridSpec(3, 5, width_ratios=[1, 1, 1, 1, 0.05])
                    cax = fig.add_subplot(gs[:, -1]) # Global colorbar for this page
                
                row, col = divmod(plot_count % 12, 4)
                ax = fig.add_subplot(gs[row, col])
                
                wave_matrix = np.vstack(waves)
                density = compute_density_matrix(wave_matrix)
                
                im = ax.imshow(
                    density.T, aspect="auto", origin="lower",
                    extent=[times.min(), times.max(), 0, density.shape[1]],
                    norm=LogNorm(vmin=1, vmax=max(2, density.max())), cmap=cmap
                )

                ax.set_title(f"Matrix {i+1}")
                ax.set_xlabel("Time")
                # Label y-axis only for the leftmost column
                if col == 0:
                    ax.set_ylabel("Mutational load")

                # Add/Update colorbar for the current page based on the latest plot
                fig.colorbar(im, cax=cax, label="Individuals (log scale)")
                
                plot_count += 1

            if fig:
                plt.tight_layout(rect=[0, 0, 0.95, 0.95])
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

    click.echo(f"[INFO] Multi-page PDF saved to {pdf_path}")

if __name__ == "__main__":
    summarize_waves()