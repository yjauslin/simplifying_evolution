rule calc_coalescent_density_fixed:
    wildcard_constraints:
        N = r"\d+",
        U = r"[\deE.+-]+",
        s = r"[\deE.+-]+"
    input:
        "results/escsim/fixed/escsim_N{N}_U{U}_s{s}.out"
    output:
        "results/coalescent_densities/N{N}_U{U}_s{s}.out"
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
            --output results/coalescent_densities \
            --mode f \
            {wildcards.N} {wildcards.s} {wildcards.U} 2>&1 | tee {log}
        """

rule calc_coalescent_density_normal:
    wildcard_constraints:
        N = r"\d+",
        U = r"[\deE.+-]+",
        s = r"[\deE.+-]+"
    input:
        "results/escsim/normal/escsim_N{N}_U{U}_s{s}_sd{sigma}.out"
    output:
        "results/coalescent_densities/N{N}_U{U}_s{s}_sd{sigma}.out"
    log:
        "logs/coalescent_density/normal/N{N}_U{U}_s{s}_sd{sigma}.log"
    params:
    resources:
        mem_mb=5*1000,        
        runtime=60,
    threads: 4
    shell:
        """
        python workflow/scripts/calc_coalescent_density.py \
            --input_folder {input} \
            --output results/coalescent_densities \
            --mode n \
            {wildcards.N} {wildcards.s} {wildcards.U} {wildcards.sigma} 2>&1 | tee {log}
        """