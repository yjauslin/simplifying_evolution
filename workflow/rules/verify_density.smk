frequency_outputs = [
    f"results/coalescent_densities/{mode}_N{p['N']}_U{p['U']}_s{p['s']}_sigma{p['sigma']:.1e}.txt"
    for mode in MODES
    for p in PARAM_COMBINATIONS
]

rule generate_trees:
    output:
        trees = temp(expand(
            "results/trees/{{mode}}_N{{N}}_U{{U}}_s{{s}}_sigma{{sigma}}_{i}.trees",
            i=range(config["constants"]["N_ITER"])
        ))
    log:
        "logs/trees/{mode}/N{N}_U{U}_s{s}_sigma{sigma}_batch.log"
    params:
        mode_flag=lambda wc: MODES[wc.mode],
        n_iter=config["constants"]["N_ITER"],
    resources:
        mem_mb= 12 * 1000,
        runtime = 180,
    threads: 5
    shell:
        """
        escsim run \
        -t -w {params.n_iter} -j 5 -m {params.mode_flag} \
        --folder results/trees/ \
        {wildcards.N} {wildcards.s} {wildcards.U} {wildcards.sigma} \
        """

rule generate_coalescent_frequencies:
    input:
        expand(
            "results/trees/{{mode}}_N{{N}}_U{{U}}_s{{s}}_sigma{{sigma}}_{i}.trees",
            i=range(config["constants"]["N_ITER"])
        )
    output:
        "results/coalescent_densities/{mode}_N{N}_U{U}_s{s}_sigma{sigma}.txt"
    log:
        "logs/frequencies/{mode}/N{N}_U{U}_s{s}_sigma{sigma}.log"
    params:
        n_iter=config["constants"]["N_ITER"],
        n_sam=100
    resources:
        mem_mb = 10 * 1000,
        runtime = 60,
    threads: 1
    shell:
        """
        python workflow/scripts/tree.py \
        -t {wildcards.mode}_N{wildcards.N}_U{wildcards.U}_s{wildcards.s}_sigma{wildcards.sigma} \
        -i results/trees/ -o results/coalescent_densities \
        {params.n_iter} {params.n_sam} > {log} 2>&1
        """
        

