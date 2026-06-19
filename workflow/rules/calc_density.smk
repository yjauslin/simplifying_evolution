rule calc_coalescent_density_fixed:
    wildcard_constraints:
        N = r"\d+",
        U = r"[\deE.+-]+",
        s = r"[\deE.+-]+"
    input:
        "results/escsim/fixed/escsim_N{N}_U{U}_s{s}.out"
    output:
        "results/coalescent_densities/fixed/N{N}_U{U}_s{s}.out"
    log:
        "logs/coalescent_density/fixed/N{N}_U{U}_s{s}.log"
    params:
    resources:
        mem_mb=5*1000,        
        runtime=60,
    threads: 4
    shell:
        """
        python workflow/scripts/calc_coalescent_density.py \
            --input_folder {input} \
            --output results/coalescent_densities/fixed \
            --mode f \
            {wildcards.N} {wildcards.s} {wildcards.U} 2>&1 | tee {log}
        """

rule calc_coalescent_density_normal:
    wildcard_constraints:
        N = r"\d+",
        U = r"[\deE.+-]+",
        s = r"[\deE.+-]+",
        sigma = r"[\deE.+-]+"
    input:
        "results/escsim/normal/escsim_N{N}_U{U}_s{s}_sd{sigma}.out"
    output:
        "results/coalescent_densities/normal/N{N}_U{U}_s{s}_sd{sigma}.out"
    log:
        "logs/coalescent_density/normal/N{N}_U{U}_s{s}_sd{sigma}.log"
    params:
        # Truncates to max 15 decimals, strips trailing zeros to preserve short values
        sigma_formatted=lambda wc: f"{float(wc.sigma):.15f}".rstrip('0').rstrip('.')
    resources:
        mem_mb=5000,
        runtime=60,
    threads: 4
    shell:
        """
        python workflow/scripts/calc_coalescent_density.py \
            --input_folder {input} \
            --output results/coalescent_densities/normal \
            --mode n \
            {wildcards.N} {wildcards.s} {wildcards.U} {params.sigma_formatted} \
            2>&1 | tee {log}
        """

# ==============================================================================
# BATCH MARKER AGGREGATIONS
# ==============================================================================

rule gather_coalescent_density_fixed:
    input:
        # Resolves both s_eff and U_eff parameters within localized scope
        s_eff_files = expand(
            "results/coalescent_densities/fixed/N{N}_U{U}_s{s}.out",
            zip,
            N=exp_data["s_eff"]["fixed"]["N"],
            U=exp_data["s_eff"]["fixed"]["U"],
            s=exp_data["s_eff"]["fixed"]["s"],
        ),
        U_eff_files = expand(
            "results/coalescent_densities/fixed/N{N}_U{U}_s{s}.out",
            zip,
            N=exp_data["U_eff"]["fixed"]["N"],
            U=exp_data["U_eff"]["fixed"]["U"],
            s=exp_data["U_eff"]["fixed"]["s"],
        )
    output:
        "results/markers/calc_coalescent_density_fixed.done"
    shell:
        "touch {output}"

rule gather_coalescent_density_normal:
    input:
        s_eff_files = expand(
            "results/coalescent_densities/normal/N{N}_U{U}_s{s}_sd{sigma}.out",
            zip,
            N=exp_data["s_eff"]["normal"]["N"],
            U=exp_data["s_eff"]["normal"]["U"],
            s=exp_data["s_eff"]["normal"]["s"],
            sigma=exp_data["s_eff"]["normal"]["sigma"],
        ),
        U_eff_files = expand(
            "results/coalescent_densities/normal/N{N}_U{U}_s{s}_sd{sigma}.out",
            zip,
            N=exp_data["U_eff"]["normal"]["N"],
            U=exp_data["U_eff"]["normal"]["U"],
            s=exp_data["U_eff"]["normal"]["s"],
            sigma=exp_data["U_eff"]["normal"]["sigma"],
        )
    output:
        "results/markers/calc_coalescent_density_normal.done"
    shell:
        "touch {output}"