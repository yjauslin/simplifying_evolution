rule calc_coalescent_density:
    input:
        "results/escsim/escsim_normal_N{N}_U{U}_s{s}_sigma{sigma}.out"
    output:
        "results/coalescent_densities/N{N}_U{U}_s{s}_sigma{sigma}.out"
    log:
        "logs/coalescent_density/N{N}_U{U}_s{s}_sigma{sigma}.log"
    shell:
        """
        python workflow/scripts/calc_coalescent_density.py \
            --input_folder {input} \
            --output results/coalescent_densities \
            {wildcards.N} {wildcards.s} {wildcards.U} {wildcards.sigma} 2>&1 | tee {log}
        """