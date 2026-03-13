density_outputs = [
    f"results/coalescent_densities/N{p['N']}_U{p['U']}_s{p['s']}_sigma{p['sigma']}.out"
    for p in PARAM_COMBINATIONS
]

wave_outputs = [
    f"results/escsim/wave_{mode}_N{p['N']}_U{p['U']}_s{p['s']}_sigma{p['sigma']}.out"
    for mode in MODES
    for p in PARAM_COMBINATIONS
]

rule visualize_wave:
    input:
        lambda wc: [w for w in wave_outputs if f"wave_{wc.mode}_" in w]
    output:
        "results/escsim_figures/wave_{mode}_summary.pdf"
    log:
        "logs/visualize_wave_{mode}.log"
    params:
        mode_flag=lambda wc: MODES[wc.mode]
    shell:
        """
        python workflow/scripts/visualize_wave.py \
            --input_folder results/escsim \
            --mode {params.mode_flag} 2>&1 | tee {log}
        """

rule visualize_densities:
    input:
        density_outputs
    output:
        "results/escsim_figures/coalescent_density.pdf",
        "results/escsim_figures/effective_population_size.pdf"
    log:
        "logs/visualize_densities.log"
    shell:
        """
        python workflow/scripts/visualize_densities.py 2>&1 | tee {log}
        """