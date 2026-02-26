from pathlib import Path

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

results_folder = Path("results/coalescent_densities")

files = sorted(results_folder.glob("*.out"))

fig1, axes1 = plt.subplots(3, 4, figsize=(16, 10))
axes1 = axes1.flatten()

fig2, axes2 = plt.subplots(3, 4, figsize=(16, 10))
axes2 = axes2.flatten()

for i, file_path in enumerate(files):
    df = pd.read_csv(file_path, sep="\t")
    N = df["popsize"].iloc[0]
    U = df["mutrate"].iloc[0]
    sel_coef = df["selcoef"].iloc[0]
    v = df["velocity"].iloc[0]
    density = np.array(df["density"].iloc[0].split(","), dtype=float)
    effective_pop_size = np.array(df["effective_pop_size"].iloc[0].split(","), dtype=float)
    time = np.array(df["time"].iloc[0].split(","), dtype=float)
    
    generation_time = N * time

    sns.lineplot(x=time, y=density, ax=axes1[i])
    axes1[i].set_title(f"N={N}, U={U}, s={sel_coef}")
    axes1[i].set_xlabel("Time (Generations)")
    axes1[i].set_ylabel("Coalescent Density")

    sns.lineplot(x=time, y=effective_pop_size, ax=axes2[i])
    axes2[i].set_title(f"N={N}, U={U}, s={sel_coef}")
    axes2[i].set_xlabel("Time (Generations)")
    axes2[i].set_ylabel("Effective Population Size (N_e)")

plt.tight_layout()

fig1.savefig("results/escsim_figures/coalescent_density.pdf")
fig2.savefig("results/escsim_figures/effective_population_size.pdf")