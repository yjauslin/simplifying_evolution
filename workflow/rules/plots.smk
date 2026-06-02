density_outputs_fixed = [
    f"results/coalescent_densities/N{p['N']}_U{p['U']}_s{p['s']}.out"
    for p in PARAM_COMBINATIONS
]

density_outputs_normal = [
    f"results/coalescent_densities/N{p['N']}_U{p['U']}_s{p['s']}_sd{p['sigma']}.out"
    for p in PARAM_COMBINATIONS
]

wave_outputs_fixed = [
    f"results/escsim_figures/fixed/wave_summary_N{p['N']}_U{p['U']}_s{p['s']}.pdf"
    for p in PARAM_COMBINATIONS
]

wave_outputs_normal = [
    f"results/escsim_figures/normal/wave_summary_N{p['N']}_U{p['U']}_s{p['s']}_sd{p['sigma']}.pdf"
    for p in PARAM_COMBINATIONS
] 

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
        "logs/visualize_wave_fixed_N{N}_U{U}_s{s}.log"
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
        "logs/visualize_wave_fixed_N{N}_U{U}_s{s}_sd{sigma}.log"
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
        density_outputs_fixed,
        density_outputs_normal,
        frequency_outputs_fixed,
        frequency_outputs_normal
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