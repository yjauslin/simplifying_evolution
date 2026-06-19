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
            --mode f > {log} 2>&1
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
            --mode n > {log} 2>&1
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