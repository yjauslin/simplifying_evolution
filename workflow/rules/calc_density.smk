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
        sigma = r"[\deE.+-]+"
    input:
        "results/escsim/normal/escsim_N{N}_U{U}_s{s}_sd{sigma}.out"
    output:
        "results/coalescent_densities/normal/N{N}_U{U}_s{s}_sd{sigma}.out"
    log:
        "logs/coalescent_density/normal/N{N}_U{U}_s{s}_sd{sigma}.log"
    params:
        s_fixed=config["experiments"]["RMSE"]["s_fixed"][0],
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
            {wildcards.N} {params.s_fixed} {wildcards.U} {params.sigma_formatted} \
            2>&1 | tee {log}
        """