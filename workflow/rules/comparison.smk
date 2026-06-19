import glob
import os
import numpy as np

def get_fixed_by_N_and_U_seff(wildcards):
    """
    Finds ONLY the N_SEL_COEF fixed files generated for the specific 
    baseline selection coefficient matching the current job.
    """
    target_N = int(wildcards.N)
    target_U_str = format_fs_flat(wildcards.U)
    
    # The wildcard {s} for the current job tells us exactly which 
    # baseline selection coefficient we are evaluating right now.
    target_s_normal = float(wildcards.s)

    # 1. Step 1: Re-calculate the specific slice of 's' values for THIS baseline
    # We mirror the exact np.linspace calculation logic from your dataframe setup script
    n_sel_coef = config["constants"]["N_SEL_COEF"]
    s_range_factor = config["experiments"]["s_eff"]["s_range"][0]
    
    # Generate ONLY the sub-values belonging to the current baseline s
    specific_s_values = np.linspace(s_range_factor * target_s_normal, target_s_normal, n_sel_coef)
    
    # 2. Step 2: Cross-reference with the DataFrame to ensure they exist
    df = exp_data["s_eff"]["fixed"]
    
    # Filter for rows matching N, U, and containing only our localized s-slice
    matched = df[
        (df["N"].astype(int) == target_N) & 
        (df["U"].astype(str) == target_U_str) & 
        (df["s"].isin(specific_s_values))
    ]
    
    matching_s_vals = list(set(matched["s"].tolist()))
    
    # 3. Return the isolated cluster of files to Snakemake
    return expand(
        "results/coalescent_densities/fixed/N{N}_U{U}_s{s}.out",
        N=wildcards.N,
        U=target_U_str,
        s=[format_fs_flat(s_val) for s_val in matching_s_vals]
    )

def get_fixed_by_N_and_s_ueff(wildcards):
    """
    Finds ONLY the N_SEL_COEF fixed files generated for the specific 
    baseline mutation rate (U) matching the current job.
    """
    target_N = int(wildcards.N)
    target_s_str = format_fs_flat(wildcards.s)
    target_mut_rate_normal = float(wildcards.U)

    # 1. Step 1: Re-calculate the specific slice of 'U' values
    n_sel_coef = config["constants"]["N_SEL_COEF"]
    u_range_factor = config["experiments"]["U_eff"]["mut_rate_range"][0]
    
    specific_U_values = np.linspace(
        u_range_factor * target_mut_rate_normal,
        (1 + u_range_factor) * target_mut_rate_normal,
        n_sel_coef
    )
    
    # CRITICAL FIX: Convert our calculated floats into the exact string 
    # format your pipeline uses for filenames and records
    specific_U_strs = [format_fs_flat(u_val) for u_val in specific_U_values]
    
    # 2. Step 2: Cross-reference with the U_eff dataframe using string comparison
    df = exp_data["U_eff"]["fixed"]
    
    # Convert df["U"] to string via format_fs_flat or astype(str) depending on how it's stored
    # Assuming standard string conversion matches your format_fs_flat output:
    matched = df[
        (df["N"].astype(int) == target_N) & 
        (df["s"].astype(str) == target_s_str) & 
        (df["U"].astype(float).apply(format_fs_flat).isin(specific_U_strs))
    ]
    
    matching_U_vals = list(set(matched["U"].tolist()))
    
    # If the dataframe filtering failed due to a structural mismatch, 
    # fallback to using the calculated range directly so the job doesn't crash empty
    if not matching_U_vals:
        return expand(
            "results/coalescent_densities/fixed/N{N}_U{U}_s{s}.out",
            N=wildcards.N,
            U=specific_U_strs,
            s=target_s_str
        )
    
    # 3. Return the isolated cluster of files to Snakemake
    return expand(
        "results/coalescent_densities/fixed/N{N}_U{U}_s{s}.out",
        N=wildcards.N,
        U=[format_fs_flat(u_val) for u_val in matching_U_vals],
        s=target_s_str
    )

def get_all_comparison_outputs_seff(wildcards):
    outputs = []
    df = exp_data["s_eff"]["normal"]
    
    matched = df[
        (df["N"].astype(str) == str(wildcards.N)) &
        (df["U"].astype(str) == str(wildcards.U)) &
        (df["s"].astype(str) == str(wildcards.s))
    ]
    
    for _, row in matched.iterrows():
        outputs.append(f"results/comparisons/s_eff/s_N{row['N']}_U{row['U']}_s{row['s']}_sd{row['sigma']}.txt")
    return outputs

def get_all_comparison_outputs_ueff(wildcards):
    outputs = []
    df = exp_data["U_eff"]["normal"]
    
    matched = df[
        (df["N"].astype(str) == str(wildcards.N)) &
        (df["U"].astype(str) == str(wildcards.U)) &
        (df["s"].astype(str) == str(wildcards.s))
    ]
    
    for _, row in matched.iterrows():
        outputs.append(f"results/comparisons/U_eff/U_N{row['N']}_U{row['U']}_s{row['s']}_sd{row['sigma']}.txt")
    return outputs

rule compare_estimates_s:
    wildcard_constraints: 
        N = r"\d+",
        U = r"[0-9eE.+-]+",
        s = r"[0-9eE.+-]+",
        sigma = r"[0-9eE.+-]+"
    input:
        fixed=get_fixed_by_N_and_U_seff,
        normal=f"results/coalescent_densities/normal/N{{N}}_U{{U}}_s{{s}}_sd{{sigma}}.out",
    output:
        f"results/comparisons/s_eff/s_N{{N}}_U{{U}}_s{{s}}_sd{{sigma}}.txt"
    log:
        f"logs/comparisons/s_eff/s_N{{N}}_U{{U}}_s{{s}}_sd{{sigma}}.log"
    threads: 1
    resources:
        mem_mb=20*1000,
        runtime=30
    shell:
        """
        python workflow/scripts/compare_distributions.py \
            -i_f '{input.fixed}' \
            -i_n {input.normal} \
            -o results/comparisons/s_eff > {log} 2>&1
        """

rule visualize_comparisons_s:
    input:
        f"results/comparisons/s_eff/s_N{{N}}_U{{U}}_s{{s}}_sd{{sigma}}.txt"
    output:
        f"results/escsim_figures/comparison_plots/s_eff/s_N{{N}}_U{{U}}_s{{s}}_sd{{sigma}}.jpg"
    log:
        f"logs/comparison_plots/s_N{{N}}_U{{U}}_s{{s}}_sd{{sigma}}.log"
    threads: 1
    resources:
        mem_mb=10*1000,
        runtime=10,
    shell:
        """
        python workflow/scripts/visualize_rmse.py \
        -i {input} \
        -o results/escsim_figures/comparison_plots/s_eff 2>&1 | tee {log}
        """

rule get_min_value_s:
    wildcard_constraints: 
        N = r"\d+",
        U = r"[0-9eE.+-]+",
        s = r"[0-9eE.+-]+"
    input:
        get_all_comparison_outputs_seff
    output:
        f"results/min_values/s_eff/s_N{{N}}_U{{U}}_s{{s}}.txt"
    log:
        f"logs/min_values/s_N{{N}}_U{{U}}_s{{s}}.log"
    threads: 1
    resources:
        mem_mb=10*1000,
        runtime=10,
    shell:
        """
        python workflow/scripts/get_min_kolmogorov.py \
        {wildcards.N} {wildcards.s} {wildcards.U} \
        -i results/comparisons/s_eff -o results/min_values/s_eff 2>&1 | tee {log}
        """

rule effective_selection_coefficient:
    input:
        # Explicitly pull ONLY from the normal track to guarantee exactly 3 files
        files=expand(
            "results/min_values/s_eff/s_N{N}_U{U}_s{s}.txt",
            zip,
            N=exp_data["s_eff"]["normal"]["N"],
            U=exp_data["s_eff"]["normal"]["U"],
            s=exp_data["s_eff"]["normal"]["s"],
        )
    output:
        "results/escsim_figures/effective_selection_coefficient.jpg"
    log:
        "logs/effective_selection_coefficient.log"
    params:
        # Dynamically extract unique baseline population size and mutation rate
        pop_size = lambda wildcards: exp_data["s_eff"]["normal"]["N"].iloc[0],
        mut_rate = lambda wildcards: exp_data["s_eff"]["normal"]["U"].iloc[0],
        # Turns your unique s values cleanly into a string like "-s 0.0004 -s 0.001 -s 0.004"
        s_flags = lambda wildcards: " ".join([f"-s {s}" for s in exp_data["s_eff"]["normal"]["s"].unique()])
    threads: 1
    resources:
        mem_mb=5000,
        runtime=10
    shell:
        """
        python workflow/scripts/visualize_sd_vs_s.py \
            {params.pop_size} \
            -u {params.mut_rate} \
            {params.s_flags} \
            -i results/min_values/s_eff \
            -o results/escsim_figures > {log} 2>&1
        """

# ==============================================================================
# U_eff
# ==============================================================================

rule compare_estimates_U:
    wildcard_constraints: 
        N = r"\d+",
        U = r"[0-9eE.+-]+",
        s = r"[0-9eE.+-]+",
        sigma = r"[0-9eE.+-]+"
    input:
        fixed=get_fixed_by_N_and_s_ueff,
        normal="results/coalescent_densities/normal/N{N}_U{U}_s{s}_sd{sigma}.out"
    output:
        "results/comparisons/U_eff/U_N{N}_U{U}_s{s}_sd{sigma}.txt"
    log:
        "logs/comparisons/U_eff/U_N{N}_U{U}_s{s}_sd{sigma}.log"
    threads: 1
    resources:
        mem_mb=20*1000,
        runtime=30
    shell:
        """
        python workflow/scripts/compare_distributions.py \
            -i_f '{input.fixed}' \
            -i_n {input.normal} \
            -o results/comparisons/U_eff \
            --type > {log} 2>&1
        """

rule visualize_comparisons_U:
    input:
        "results/comparisons/U_eff/U_N{N}_U{U}_s{s}_sd{sigma}.txt"
    output:
        "results/escsim_figures/comparison_plots/U_eff/U_N{N}_U{U}_s{s}_sd{sigma}.jpg"
    log:
        "logs/comparison_plots/U_eff/U_N{N}_U{U}_s{s}_sd{sigma}.log"
    threads: 1
    resources:
        mem_mb=10*1000,
        runtime=10
    shell:
        """
        python workflow/scripts/visualize_rmse.py \
            -i {input} \
            -o results/escsim_figures/comparison_plots/U_eff 2>&1 | tee {log}
        """

rule get_min_value_U:
    wildcard_constraints: 
        N = r"\d+",
        U = r"[0-9eE.+-]+",
        s = r"[0-9eE.+-]+"
    input:
        get_all_comparison_outputs_ueff 
    output:
        "results/min_values/U_eff/U_N{N}_U{U}_s{s}.txt"
    log:
        "logs/min_values/U_eff/U_N{N}_U{U}_s{s}.log"
    threads: 1
    resources:
        mem_mb=10*1000,
        runtime=10
    shell:
        """
        python workflow/scripts/get_min_kolmogorov.py \
            {wildcards.N} {wildcards.s} {wildcards.U} \
            -i results/comparisons/U_eff \
            -o results/min_values/U_eff \
            --type > {log} 2>&1
        """

rule effective_mutation_rate:
    input:
        # Require all 3 baseline min_value files before running the plot
        files=expand(
            "results/min_values/U_eff/U_N{N}_U{U}_s{s}.txt",
            zip,
            N=exp_data["U_eff"]["normal"]["N"],
            U=exp_data["U_eff"]["normal"]["U"],
            s=exp_data["U_eff"]["normal"]["s"],
        )
    output:
        "results/escsim_figures/effective_mutation_rate.jpg"
    log:
        "logs/effective_mutation_rate.log"
    params:
        # Helper to extract unique population sizes and selection coefficients
        pop_size = lambda wildcards: exp_data["U_eff"]["normal"]["N"].iloc[0],
        sel_coef = lambda wildcards: exp_data["U_eff"]["normal"]["s"].iloc[0],
        # Turns the U values list cleanly into a string like "-u 0.012 -u 0.006 -u 0.003"
        u_flags = lambda wildcards: " ".join([f"-u {u}" for u in exp_data["U_eff"]["normal"]["U"].unique()])
    threads: 1
    resources:
        mem_mb=5000,
        runtime=10
    shell:
        """
        python workflow/scripts/visualize_sd_vs_s.py \
            {params.pop_size} \
            -s {params.sel_coef} \
            {params.u_flags} \
            -i results/min_values/U_eff \
            -o results/escsim_figures > {log} 2>&1
        """

# ==============================================================================
# BATCH MARKER AGGREGATIONS
# ==============================================================================

rule gather_comparisons_s_eff_normal:
    input:
        expand(
            "results/comparisons/s_eff/s_N{N}_U{U}_s{s}_sd{sigma}.txt",
            zip,
            N=exp_data["s_eff"]["normal"]["N"],
            U=exp_data["s_eff"]["normal"]["U"],
            s=exp_data["s_eff"]["normal"]["s"],
            sigma=exp_data["s_eff"]["normal"]["sigma"],
        )
    output:
        "results/markers/comparisons_s_eff_normal.done"
    shell:
        "touch {output}"

rule gather_comparisons_U_eff_normal:
    input:
        expand(
            "results/comparisons/U_eff/U_N{N}_U{U}_s{s}_sd{sigma}.txt",
            zip,
            N=exp_data["U_eff"]["normal"]["N"],
            U=exp_data["U_eff"]["normal"]["U"],
            s=exp_data["U_eff"]["normal"]["s"],
            sigma=exp_data["U_eff"]["normal"]["sigma"],
        )
    output:
        "results/markers/comparisons_U_eff_normal.done"
    shell:
        "touch {output}"

rule gather_min_values_s_eff_fixed:
    input:
        expand(
            "results/min_values/s_eff/s_N{N}_U{U}_s{s}.txt",
            zip,
            N=exp_data["s_eff"]["normal"]["N"].unique().tolist() * len(exp_data["s_eff"]["normal"]["s"].unique()),
            U=exp_data["s_eff"]["normal"]["U"].repeat(len(exp_data["s_eff"]["normal"]["s"].unique())),
            s=exp_data["s_eff"]["normal"]["s"].tolist() * len(exp_data["s_eff"]["normal"]["U"].unique()),
        )
    output:
        "results/markers/min_values_s_eff_fixed.done"
    shell:
        "touch {output}"

rule gather_min_values_U_eff_fixed:
    input:
        expand(
            "results/min_values/U_eff/U_N{N}_U{U}_s{s}.txt",
            zip,
            N=exp_data["U_eff"]["normal"]["N"].unique().tolist() * len(exp_data["U_eff"]["normal"]["s"].unique()),
            U=exp_data["U_eff"]["normal"]["U"].repeat(len(exp_data["U_eff"]["normal"]["s"].unique())),
            s=exp_data["U_eff"]["normal"]["s"].tolist() * len(exp_data["U_eff"]["normal"]["U"].unique()),
        )
    output:
        "results/markers/min_values_U_eff_fixed.done"
    shell:
        "touch {output}"