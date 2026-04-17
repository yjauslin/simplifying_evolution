density_outputs = [
    f"results/coalescent_densities/{mode}_N{p['N']}_U{p['U']}_s{p['s']}_sigma{p['sigma']}.out"
    for mode in MODES
    for p in PARAM_COMBINATIONS
]

def get_wave_input(wildcards):
    return f"results/escsim/{wildcards.mode}/wave_N{wildcards.N}_U{wildcards.U}_s{wildcards.s}_sigma{wildcards.sigma}.out"

wave_outputs = [
    f"results/escsim_figures/{mode}/wave_summary_N{p['N']}_U{p['U']}_s{p['s']}_sigma{p['sigma']}.pdf"
    for mode in MODES
    for p in PARAM_COMBINATIONS
]

rule visualize_wave:
    input:
        # Marking the input as temporary here ensures that once this specific 
        # rule finishes, the massive .out file is deleted.
        wave_file = get_wave_input
    output:
        pdf = "results/escsim_figures/{mode}/wave_summary_N{N}_U{U}_s{s}_sigma{sigma}.pdf"
    log:
        "logs/visualize_wave_{mode}_N{N}_U{U}_s{s}_sigma{sigma}.log"
    params:
        mode_flag = lambda wc: wc.mode[0] # e.g., 'f' or 'n'
    resources:
        mem_mb = 12 * 1000,
        runtime = 60,
    threads: 1
    shell:
        """
        python workflow/scripts/visualize_wave.py \
            {input.wave_file} \
            --output_folder results/escsim_figures/{wildcards.mode} \
            --mode {params.mode_flag} > {log} 2>&1
        """

rule visualize_densities:
    input:
        density_outputs,
        frequency_outputs
    output:
        "results/escsim_figures/coalescent_density.pdf",
        "results/escsim_figures/effective_population_size.pdf"
    log:
        "logs/visualize_densities.log"
    resources:
        mem_mb=50*1000,        
        runtime=120,
    threads: 10
    shell:
        """
        python workflow/scripts/visualize_densities.py 2>&1 | tee {log}
        """