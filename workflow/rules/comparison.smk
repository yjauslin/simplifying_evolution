import glob
import os

S_FIXED = config["experiments"]["RMSE"]["s_fixed"][0]

comparison_outputs_s = [
    f"results/comparisons/s_N{p['N']}_U{p['U']}_s{S_FIXED}_sd{clean_sigma(p['sigma'])}.txt"
    for p in NORMAL_PARAM_COMBINATIONS]

comparison_plots_s = [
    f"results/escsim_figures/comparison_plots/s_N{p['N']}_U{p['U']}_s{S_FIXED}_sd{clean_sigma(p['sigma'])}.jpg"
    for p in NORMAL_PARAM_COMBINATIONS]

min_values_s = [
    f"results/min_values/s_N{p['N']}_U{p['U']}_s{S_FIXED}.txt"
    for p in NORMAL_PARAM_COMBINATIONS]

sigma_values = sorted({x['sigma'] for x in NORMAL_PARAM_COMBINATIONS})

rule compare_estimates_s:
    wildcard_constraints: 
        N = r"\d+",
        U = r"[\deE.+-]+",
        sigma = r"[\deE.+-]+"
    input:
        fixed=lambda wildcards: glob.glob(f"results/coalescent_densities/fixed/N{wildcards.N}_U{wildcards.U}_s*.out"),
        normal=f"results/coalescent_densities/normal/N{{N}}_U{{U}}_s{S_FIXED}_sd{{sigma}}.out",
    output:
        f"results/comparisons/s_N{{N}}_U{{U}}_s{S_FIXED}_sd{{sigma}}.txt"
    log:
        f"logs/comparisons/s_N{{N}}_U{{U}}_s{S_FIXED}_sd{{sigma}}.log"
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

rule visualize_comparisons_s:
    input:
        f"results/comparisons/s_N{{N}}_U{{U}}_s{S_FIXED}_sd{{sigma}}.txt"
    output:
        f"results/escsim_figures/comparison_plots/s_N{{N}}_U{{U}}_s{S_FIXED}_sd{{sigma}}.jpg"
    log:
        f"logs/comparison_plots/s_N{{N}}_U{{U}}_s{S_FIXED}_sd{{sigma}}.log"
    threads: 1
    resources:
        mem_mb=10*1000,
        runtime=10,
    shell:
        """
        python workflow/scripts/visualize_rmse.py \
        -i results/comparisons/s_N{wildcards.N}_U{wildcards.U}_s{S_FIXED}_sd{wildcards.sigma}.txt \
        -o results/escsim_figures/comparison_plots 2>&1 | tee {log}
        """

rule get_min_value_s:
    input:
        lambda wildcards: expand(
            f"results/comparisons/s_N{{N}}_U{{U}}_s{S_FIXED}_sd{{sigma}}.txt",
            N=wildcards.N,
            U=wildcards.U,
            sigma=sigma_values
        )
    output:
        f"results/min_values/s_N{{N}}_U{{U}}_s{S_FIXED}.txt"
    log:
        f"logs/min_values/s_N{{N}}_U{{U}}_s{S_FIXED}.log"
    threads: 1
    resources:
        mem_mb=10*1000,
        runtime=10,
    shell:
        """
        python workflow/scripts/get_min_kolmogorov.py \
        {wildcards.N} {S_FIXED} {wildcards.U} \
        -i results/comparisons -o results/min_values 2>&1 | tee {log}
        """

rule effective_selection_coefficient:
    input:
        min_values_s
    output:
        f"results/escsim_figures/effective_selection_coefficient.jpg"
    log:
        f"logs/effective_selection_coefficient.log"
    threads: 1
    resources:
        mem_mb=5*1000,
        runtime=10,
    shell:
        """
        python workflow/scripts/visualize_sd_vs_s.py \
        -s {S_FIXED} \
        -i results/min_values -o results/escsim_figures
        """
