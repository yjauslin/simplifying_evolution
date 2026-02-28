import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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

times, waves = load_wave_file("tmp/results/wave_N5000_U0.006_s0.0004.out")
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
plt.title("Mutation wave dynamics")
plt.tight_layout()
plt.show()

