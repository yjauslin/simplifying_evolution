rule calc_coalescent_density:
    input:
        "results/escsim/{mode}/escsim_N{N}_U{U}_s{s}_sigma{sigma}.out"
    output:
        "results/coalescent_densities/{mode}_N{N}_U{U}_s{s}_sigma{sigma}.out"
    log:
        "logs/coalescent_density/{mode}/N{N}_U{U}_s{s}_sigma{sigma}.log"
    params:
        mode_flag=lambda wc: MODES[wc.mode]
    resources:
        mem_mb=10*1024,        
        runtime=120,
    threads: 4
    shell:
        """
        python workflow/scripts/calc_coalescent_density.py \
            --input_folder {input} \
            --output results/coalescent_densities \
            --mode {params.mode_flag} \
            {wildcards.N} {wildcards.s} {wildcards.U} {wildcards.sigma} 2>&1 | tee {log}
        """