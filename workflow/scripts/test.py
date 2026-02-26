from pathlib import Path

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

results_folder = Path("results/coalescent_densities")

df = pd.read_csv(f"{results_folder}/N5000_U0.006_s0.01.out", sep="\t")
N = df["popsize"].iloc[0]
U = df["mutrate"].iloc[0]
sel_coef = df["selcoef"].iloc[0]
v = df["velocity"].iloc[0]
density = np.array(df["density"].iloc[0].split(","), dtype=float)
effective_pop_size = np.array(df["effective_pop_size"].iloc[0].split(","), dtype=float)
time = np.array(df["time"].iloc[0].split(","), dtype=float)
    
generation_time = N * time

sns.lineplot(x=time, y=density)
plt.title(f"N={N}, U={U}, s={sel_coef}")
plt.xlabel("Time (Generations)")
plt.ylabel("Coalescent Density")
# plt.ylim(0, 0.0003)

plt.show()

profile_df =  pd.read_csv("results/escsim/escsim_N5000_U0.006_s0.01.out", sep="\t")
profile = np.array(profile_df["profile"].iloc[0].split(","), dtype=float)
plt.bar(x=np.arange(len(profile)), height=profile)
plt.xlabel("Mutation Burden Class")
plt.ylabel("Frequency")
plt.title(f"N={N}, U={U}, s={sel_coef}")
plt.show()