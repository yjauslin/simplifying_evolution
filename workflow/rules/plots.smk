rule visualize_wave_fixed:
    wildcard_constraints:
        N = r"\d+",
        U = r"[\deE.+-]+",
        s = r"[\deE.+-]+"
    input:
        "results/escsim/fixed/wave_N{N}_U{U}_s{s}.out"
    output:
        pdf = "results/escsim_figures/fixed/wave_summary_N{N}_U{U}_s{s}.pdf"
    log:
        "logs/visualize_wave/fixed/N{N}_U{U}_s{s}.log"
    params:
    resources:
        mem_mb = 12 * 1000,
        runtime = 60,
    threads: 1
    shell:
        """
        python workflow/scripts/visualize_wave.py \
            {input} \
            --output_folder results/escsim_figures/fixed \
            --limit 20 > {log} 2>&1
        """

rule visualize_wave_normal:
    wildcard_constraints:
        N = r"\d+",
        U = r"[\deE.+-]+",
        s = r"[\deE.+-]+"
    input:
        "results/escsim/normal/wave_N{N}_U{U}_s{s}_sd{sigma}.out"
    output:
        pdf = "results/escsim_figures/normal/wave_summary_N{N}_U{U}_s{s}_sd{sigma}.pdf"
    log:
        "logs/visualize_wave/normal/N{N}_U{U}_s{s}_sd{sigma}.log"
    params:
    resources:
        mem_mb = 12 * 1000,
        runtime = 60,
    threads: 1
    shell:
        """
        python workflow/scripts/visualize_wave.py \
            {input} \
            --output_folder results/escsim_figures/normal \
            --limit 20 > {log} 2>&1
        """

rule visualize_densities:
    input:
        # OPTIMIZATION: Instead of forcing the DAG engine to parse millions of strings
        # via massive expand configurations, we anchor this rule to the cluster aggregation points.
        # This prevents the master thread from experiencing memory exhaustion.
        "results/markers/calc_coalescent_density_fixed.done",
        "results/markers/calc_coalescent_density_normal.done"
    output:
        "results/escsim_figures/coalescent_density.jpg",
        "results/escsim_figures/effective_population_size.jpg"
    log:
        "logs/visualize_densities.log"
    resources:
        mem_mb=10*1000,        
        runtime=30,
    threads: 10
    shell:
        """
        python workflow/scripts/visualize_densities.py 2>&1 | tee {log}
        """

rule effective_selection_coefficient:
    input:
        files=expand(
            "results/min_values/s_eff/s_N{N}_U{U}_s{s}.txt",
            zip,
            N=exp_data["s_eff"]["normal"]["N"],
            U=exp_data["s_eff"]["normal"]["U"],
            s=exp_data["s_eff"]["normal"]["s"],
        )
    output:
        "results/escsim_figures/effective_selection_coefficient.jpg",
        "results/escsim_figures/relative_effective_selection_coefficient.jpg"
    log:
        "logs/effective_selection_coefficient.log"
    params:
        # Dynamically extract unique baseline population size and mutation rate
        pop_size = lambda wildcards: exp_data["s_eff"]["normal"]["N"].iloc[0],
        mut_rate = lambda wildcards: exp_data["s_eff"]["normal"]["U"].iloc[0],
        # Turns unique s values cleanly into a string like "-s 0.0004 -s 0.001 -s 0.004"
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
        "results/escsim_figures/effective_mutation_rate.jpg",
        "results/escsim_figures/relative_effective_mutation_rate.jpg"
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

rule s_eff_RMSE:
    input:
        files=expand(
            "results/min_values/s_eff/s_N{N}_U{U}_s{s}.txt",
            zip,
            N=exp_data["s_eff"]["normal"]["N"],
            U=exp_data["s_eff"]["normal"]["U"],
            s=exp_data["s_eff"]["normal"]["s"],
        )
    output:
        "results/escsim_figures/s_eff_RMSE.jpg"
    log:
        "logs/s_eff_RMSE.log"
    params:
        # Dynamically extract unique baseline population size and mutation rate
        pop_size = lambda wildcards: exp_data["s_eff"]["normal"]["N"].iloc[0],
        mut_rate = lambda wildcards: exp_data["s_eff"]["normal"]["U"].iloc[0],
        # Turns unique s values cleanly into a string like "-s 0.0004 -s 0.001 -s 0.004"
        s_flags = lambda wildcards: " ".join([f"-s {s}" for s in exp_data["s_eff"]["normal"]["s"].unique()])
    threads: 1
    resources:
        mem_mb=5000,
        runtime=10
    shell:
        """
        python workflow/scripts/visualize_sd_vs_rmse.py \
            {params.pop_size} \
            -u {params.mut_rate} \
            {params.s_flags} \
            -i results/min_values/s_eff \
            -o results/escsim_figures > {log} 2>&1
        """

rule U_eff_RMSE:
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
        "results/escsim_figures/U_eff_RMSE.jpg"
    log:
        "logs/U_eff_RMSE.log"
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
        python workflow/scripts/visualize_sd_vs_rmse.py \
            {params.pop_size} \
            -s {params.sel_coef} \
            {params.u_flags} \
            -i results/min_values/U_eff \
            -o results/escsim_figures > {log} 2>&1
        """

def get_specific_mut_burden_inputs(wildcards):
    """
    selects the input files to generate the mut_burden_dist plot
    """
    selected_files = []
    
    # Define the precise fractions we want to display on the plot
    target_fractions = [0.0, 0.5, 1.0]
    
    # 1. Filter s_eff experiment files
    for s in config["experiments"]["s_eff"]["s_normal"]:
        for frac in target_fractions:
            # Replicate the exact math + string formatting your simulation pipeline used
            sigma_val = format_fs_flat(frac * s)
            u_val = format_fs_flat(config["experiments"]["s_eff"]["mut_rate"][0])
            s_val = format_fs_flat(s)
            N_val = config["experiments"]["s_eff"]["pop_size"][0]
            
            selected_files.append(
                f"results/escsim/normal/escsim_N{N_val}_U{u_val}_s{s_val}_sd{sigma_val}.out"
            )
            
    # 2. Filter U_eff experiment files
    for U in config["experiments"]["U_eff"]["mut_rate_normal"]:
        s_base = config["experiments"]["U_eff"]["sel_coef"][0]
        for frac in target_fractions:
            sigma_val = format_fs_flat(frac * s_base)
            u_val = format_fs_flat(U)
            s_val = format_fs_flat(s_base)
            N_val = config["experiments"]["U_eff"]["pop_size"][0]
            
            selected_files.append(
                f"results/escsim/normal/escsim_N{N_val}_U{u_val}_s{s_val}_sd{sigma_val}.out"
            )
            
    return selected_files

rule mut_burden_dist:
    input:
        files=get_specific_mut_burden_inputs
    output:
        "results/escsim_figures/mut_burden_dist.jpg"
    log:
        "logs/mut_burden_dist.log"
    params:
        # Extract unique values and convert them into comma-separated strings
        pop_size = int(config["experiments"]["s_eff"]["pop_size"][0]),
        
        # Pull distinct selection coefficients and mutation rates to layout the grid axes
        s_eff_s = lambda w: ",".join(map(str, config["experiments"]["s_eff"]["s_normal"])),
        s_eff_u = lambda w: ",".join(map(str, config["experiments"]["s_eff"]["mut_rate"])),
        
        ueff_s  = lambda w: ",".join(map(str, config["experiments"]["U_eff"]["sel_coef"])),
        ueff_u  = lambda w: ",".join(map(str, config["experiments"]["U_eff"]["mut_rate_normal"])),
        
        # Collect all specific formatted sigmas required for the line files parsing
        sigmas  = lambda w: ",".join(sorted(list(set(
            [format_fs_flat(f * s) for s in config["experiments"]["s_eff"]["s_normal"] for f in [0.0, 0.25, 0.5, 0.75, 1.0]] +
            [format_fs_flat(f * config["experiments"]["U_eff"]["sel_coef"][0]) for f in [0.0, 0.25, 0.5, 0.75, 1.0]]
        )))),
    threads: 1
    resources:
        mem_mb=5000,
        runtime=10,
    shell:
        """
        python workflow/scripts/visualize_mut_burden_dist.py \
        -i results/escsim/normal \
        --pop_size {params.pop_size} \
        --mut_rate {params.s_eff_u},{params.ueff_u} \
        --sel_coef {params.s_eff_s},{params.ueff_s} \
        --sigma {params.sigmas} \
        -o results/escsim_figures > {log} 2>&1
        """

rule sel_coeff_velocity:
    input:
        files=expand(
            "results/escsim/normal/escsim_N{N}_U{U}_s{s}_sd{sigma}.out",
            zip,
            N=exp_data["s_eff"]["normal"]["N"],
            U=exp_data["s_eff"]["normal"]["U"],
            s=exp_data["s_eff"]["normal"]["s"],
            sigma=exp_data["s_eff"]["normal"]["sigma"],
        )
    output:
        "results/escsim_figures/s_velocity.jpg"
    log:
        "logs/s_velocity.log"
    params:
        # Dynamically extract unique baseline population size and mutation rate
        pop_size = lambda wildcards: exp_data["s_eff"]["normal"]["N"].iloc[0],
        mut_rate = lambda wildcards: exp_data["s_eff"]["normal"]["U"].iloc[0],
        # Turns unique s values cleanly into a string like "-s 0.0004 -s 0.001 -s 0.004"
        s_flags = lambda wildcards: " ".join([f"-s {s}" for s in exp_data["s_eff"]["normal"]["s"].unique()])
    threads: 1
    resources:
        mem_mb=5000,
        runtime=10
    shell:
        """
        python workflow/scripts/visualize_velocity_vs_sd.py \
            {params.pop_size} \
            {params.s_flags} \
            -u {params.mut_rate} \
            -i results/escsim/normal \
            -o results/escsim_figures > {log} 2>&1
        """

rule mut_rate_velocity:
    input:
        files=expand(
            "results/escsim/normal/escsim_N{N}_U{U}_s{s}_sd{sigma}.out",
            zip,
            N=exp_data["U_eff"]["normal"]["N"],
            U=exp_data["U_eff"]["normal"]["U"],
            s=exp_data["U_eff"]["normal"]["s"],
            sigma=exp_data["U_eff"]["normal"]["sigma"],
        )
    output:
        "results/escsim_figures/u_velocity.jpg"
    log:
        "logs/u_velocity.log"
    params:
        # Dynamically extract unique baseline population size and mutation rate
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
        python workflow/scripts/visualize_velocity_vs_sd.py \
            {params.pop_size} \
            -s {params.sel_coef} \
             {params.u_flags} \
            -i results/escsim/normal \
            -o results/escsim_figures > {log} 2>&1
        """

def get_required_files_seff(wildcards):
    """Gathers all required .out and .txt files for the s_eff experiment."""
    files = []
    # Fixed files (.out and .txt)
    df_normal = exp_data["s_eff"]["normal"]

    sd_multipliers = [0.0, 0.5, 1.0]
    for _, row in df_normal.iterrows():
        # Estimates
        files.append(f"results/coalescent_densities/fixed/N{row['N']}_U{row['U']}_s{row['s']}.out")
        files.append(f"results/coalescent_densities/normal/N{row['N']}_U{row['U']}_s{row['s']}_sd{row['sigma']}.out")

        for multiplier in sd_multipliers:
            # Calculate the scaled standard deviation value
            calculated_sd = multiplier * float(row['s'])

            sd_str = f"{calculated_sd}"
            # WF-Simulations
            files.append(f"results/coalescent_densities/normal/N{row['N']}_U{row['U']}_s{row['s']}_sd{sd_str}.txt")
        
    return files

def get_required_files_ueff(wildcards):
    """Gathers all required .out and .txt files for the U_eff experiment."""
    files = []
    # Fixed files (.out and .txt)
    df_normal = exp_data["U_eff"]["normal"]

    sd_multipliers = [0.0, 0.5, 1.0]
    for _, row in df_normal.iterrows():
        # Estimates
        files.append(f"results/coalescent_densities/fixed/N{row['N']}_U{row['U']}_s{row['s']}.out")
        files.append(f"results/coalescent_densities/normal/N{row['N']}_U{row['U']}_s{row['s']}_sd{row['sigma']}.out")

        for multiplier in sd_multipliers:
            # Calculate the scaled standard deviation value
            calculated_sd = multiplier * float(row['s'])

            sd_str = f"{calculated_sd}"
            # WF-Simulations
            files.append(f"results/coalescent_densities/normal/N{row['N']}_U{row['U']}_s{row['s']}_sd{sd_str}.txt")
        
    return files

rule verify_s_eff:
    input:
        "results/markers/min_values_s_eff_fixed.done",
        data_files = get_required_files_seff
    output:
        "results/escsim_figures/coalescent_density_s_eff.jpg"
    log:
        "logs/verify_s_eff.log"
    params:
        pop_size = lambda wildcards: exp_data["s_eff"]["normal"]["N"].iloc[0],
        mut_rate = lambda wildcards: exp_data["s_eff"]["normal"]["U"].iloc[0],
        # Expands selection values cleanly into CLI tokens: "-s 0.0004 -s 0.001 -s 0.004"
        s_flags  = lambda wildcards: "-s " + ",".join([str(s) for s in exp_data["s_eff"]["normal"]["s"].unique()]),
        # Formats list into distinct parameter tokens: "-sd 0.0 -sd 0.5 -sd 1.0"
        sd_flags = "-sd " + ",".join([str(sig) for sig in [0.0, 0.5, 1.0]])
    threads: 1
    resources:
        mem_mb=5000,
        runtime=10
    shell:
        """
        python workflow/scripts/visualize_densities.py \
            -N {params.pop_size} \
            {params.s_flags} \
            -u {params.mut_rate} \
            {params.sd_flags} \
            -i results/coalescent_densities \
            -o results/escsim_figures > {log} 2>&1
        """

rule verify_U_eff:
    input:
        "results/markers/min_values_U_eff_fixed.done",
        data_files=get_required_files_ueff
    output:
        "results/escsim_figures/coalescent_density_u_eff.jpg"
    log:
        "logs/verify_U_eff.log"
    params:
        pop_size  = lambda wildcards: exp_data["U_eff"]["normal"]["N"].iloc[0],
        sel_coef  = lambda wildcards: exp_data["U_eff"]["normal"]["s"].iloc[0],
        # Expands mutation values cleanly into CLI tokens: "-u 0.006 -u 0.01 -u 0.0001"
        u_flags   = lambda wildcards: "-u " + ",".join([str(u) for u in exp_data["U_eff"]["normal"]["U"].unique()]),
        sd_flags = "-sd " + ",".join([str(sig) for sig in [0.0, 0.5, 1.0]])
    threads: 1
    resources:
        mem_mb=5000,
        runtime=10
    shell:
        """
        python workflow/scripts/visualize_densities.py \
            -N {params.pop_size} \
            -s {params.sel_coef} \
            {params.u_flags} \
            {params.sd_flags} \
            -i results/coalescent_densities \
            -o results/escsim_figures > {log} 2>&1
        """

# ==============================================================================
# BATCH MARKER AGGREGATIONS (Prevents Master DAG Memory Bloat)
# ==============================================================================

rule gather_wave_summary_s_eff_fixed:
    input:
        expand(
            "results/escsim_figures/fixed/wave_summary_N{N}_U{U}_s{s}.pdf",
            zip,
            N=exp_data["s_eff"]["fixed"]["N"],
            U=exp_data["s_eff"]["fixed"]["U"],
            s=exp_data["s_eff"]["fixed"]["s"],
        )
    output:
        "results/markers/wave_summary_s_eff_fixed.done"
    shell:
        "touch {output}"

rule gather_wave_summary_U_eff_fixed:
    input:
        expand(
            "results/escsim_figures/fixed/wave_summary_N{N}_U{U}_s{s}.pdf",
            zip,
            N=exp_data["U_eff"]["fixed"]["N"],
            U=exp_data["U_eff"]["fixed"]["U"],
            s=exp_data["U_eff"]["fixed"]["s"],
        )
    output:
        "results/markers/wave_summary_U_eff_fixed.done"
    shell:
        "touch {output}"

rule gather_wave_summary_s_eff_normal:
    input:
        expand(
            "results/escsim_figures/normal/wave_summary_N{N}_U{U}_s{s}_sd{sigma}.pdf",
            zip,
            N=exp_data["s_eff"]["normal"]["N"],
            U=exp_data["s_eff"]["normal"]["U"],
            s=exp_data["s_eff"]["normal"]["s"],
            sigma=exp_data["s_eff"]["normal"]["sigma"],
        )
    output:
        "results/markers/wave_summary_s_eff_normal.done"
    shell:
        "touch {output}"

rule gather_wave_summary_U_eff_normal:
    input:
        expand(
            "results/escsim_figures/normal/wave_summary_N{N}_U{U}_s{s}_sd{sigma}.pdf",
            zip,
            N=exp_data["U_eff"]["normal"]["N"],
            U=exp_data["U_eff"]["normal"]["U"],
            s=exp_data["U_eff"]["normal"]["s"],
            sigma=exp_data["U_eff"]["normal"]["sigma"],
        )
    output:
        "results/markers/wave_summary_U_eff_normal.done"
    shell:
        "touch {output}"