import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages
import click
import re

def load_multi_wave_file(filename):
    """
    Yields (times, waves) tuples. 
    Detects a new matrix when the 'time' value decreases (reset).
    """
    current_times, current_waves = [], []
    last_t = -1

    with open(filename) as f:
        next(f)  # skip the header
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 2: continue
            
            t = int(parts[0])
            wave_data = np.fromstring(parts[1], sep=",", dtype=int)

            if t < last_t and current_times:
                yield np.array(current_times), np.array(current_waves)
                current_times, current_waves = [], []
            
            current_times.append(t)
            current_waves.append(wave_data)
            last_t = t
            
        if current_times:
            yield np.array(current_times), np.array(current_waves)

def compute_density_matrix(wave_matrix):
    max_load = wave_matrix.max()
    T = wave_matrix.shape[0]
    density = np.zeros((T, max_load + 1), dtype=int)
    for i in range(T):
        density[i] = np.bincount(wave_matrix[i], minlength=max_load + 1)
    return density

def create_new_page(file_name, cmap):
    """Helper to initialize a new figure page."""
    fig = plt.figure(figsize=(18, 12))
    fig.suptitle(f"File: {file_name}", fontsize=16, fontweight='bold')
    gs = plt.GridSpec(3, 5, width_ratios=[1, 1, 1, 1, 0.05])
    cax = fig.add_subplot(gs[:, -1]) 
    return fig, gs, cax

@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('--output_folder', '-o', default='results/escsim_figures', help="Directory to save the PDF to.")
@click.option('--mode', '-m', default='f', help="Mode 'n' for normal, 'f' for fixed")
def summarize_waves(input_file, output_folder, mode):
    file_path = Path(input_file)
    os.makedirs(output_folder, exist_ok=True)
    
    # Generate output filename based on the input filename and mode
    match = re.search(r"(N\d+.*)", file_path.stem)
    if match:
        params_suffix = match.group(1)
    else:
        # Fallback if the pattern isn't found
        params_suffix = file_path.stem
    out_name = f"wave_summary_{params_suffix}.pdf"
    pdf_path = Path(output_folder) / out_name

    cmap = plt.cm.inferno
    cmap.set_under("white")

    with PdfPages(pdf_path) as pdf:
        fig = None
        plot_count = 0
        
        # Process the single specified file
        for i, (times, wave_matrix) in enumerate(load_multi_wave_file(file_path)):
            
            # If we hit 12 plots or it's the very first matrix
            if plot_count % 12 == 0:
                if fig:
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close(fig)
                
                fig, gs, cax = create_new_page(file_path.name, cmap)

            row, col = divmod(plot_count % 12, 4)
            ax = fig.add_subplot(gs[row, col])
            
            density = compute_density_matrix(wave_matrix)
            
            im = ax.imshow(
                density.T, aspect="auto", origin="lower",
                extent=[times.min(), times.max(), 0, density.shape[1]],
                norm=LogNorm(vmin=1, vmax=max(2, density.max())), cmap=cmap
            )

            ax.set_title(f"Matrix {i+1}")
            ax.set_xlabel("Time")
            if col == 0:
                ax.set_ylabel("Mutational load")

            plt.colorbar(im, cax=cax, label="Individuals (log scale)")
            
            plot_count += 1
            del wave_matrix
            del density

        # Save and Close the last page
        if fig:
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    click.echo(f"[INFO] Multi-page PDF saved to {pdf_path}")

if __name__ == "__main__":
    summarize_waves()