import glob
import os

S_FIXED = config["experiments"]["RMSE"]["s_fixed"][0]

comparison_outputs = [
    f"results/comparisons/N{p['N']}_U{p['U']}_s{S_FIXED}_sd{clean_sigma(p['sigma'])}.txt"
    for p in NORMAL_PARAM_COMBINATIONS]

comparison_plots = [
    f"results/escsim_figures/comparison_plots/N{p['N']}_U{p['U']}_s{S_FIXED}_sd{clean_sigma(p['sigma'])}.jpg"
    for p in NORMAL_PARAM_COMBINATIONS]

rule compare_estimates:
    wildcard_constraints: 
        N = r"\d+",
        U = r"[\deE.+-]+",
        sigma = r"[\deE.+-]+"
    input:
        fixed=lambda wildcards: glob.glob(f"results/coalescent_densities/fixed/N{wildcards.N}_U{wildcards.U}_s*.out"),
        normal=f"results/coalescent_densities/normal/N{{N}}_U{{U}}_s{S_FIXED}_sd{{sigma}}.out",
    output:
        f"results/comparisons/N{{N}}_U{{U}}_s{S_FIXED}_sd{{sigma}}.txt"
    log:
        f"logs/comparisons/N{{N}}_U{{U}}_s{S_FIXED}_sd{{sigma}}.log"
    threads: 1
    resources:
        mem_mb=20*1000,
        runtime=30
    shell:
        """
        python workflow/scripts/compare_distributions.py \
            -i_n {input.normal} \
            -o results/comparisons/ \
            {input.fixed} 2>&1 | tee {log}
        """

rule visualize_comparisons:
    input:
        f"results/comparisons/N{{N}}_U{{U}}_s{S_FIXED}_sd{{sigma}}.txt"
    output:
        f"results/escsim_figures/comparison_plots/N{{N}}_U{{U}}_s{S_FIXED}_sd{{sigma}}.jpg"
    log:
        f"logs/comparison_plots/N{{N}}_U{{U}}_s{S_FIXED}_sd{{sigma}}.log"
    threads: 1
    resources:
        mem_mb=10*1000,
        runtime=10,
    shell:
        """
        python workflow/scripts/visualize_rmse.py \
        -i results/comparisons/N{wildcards.N}_U{wildcards.U}_s{S_FIXED}_sd{wildcards.sigma}.txt \
        -o results/escsim_figures/comparison_plots > {log} 2>&1
        """
