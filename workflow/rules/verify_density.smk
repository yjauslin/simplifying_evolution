frequency_outputs = [
    f"results/coalescent_densities/{mode}/N{p['N']}_U{p['U']}_s{p['s']}_sigma{p['sigma']}.txt"
    for mode in MODES
    for p in PARAM_COMBINATIONS
]

rule generate_trees:
    output:
        trees = temp(expand(
            "results/trees/{{mode}}/N{{N}}_U{{U}}_s{{s}}_sigma{{sigma}}_{i}.trees",
            i=range(100)
    log:
        "logs/trees/{mode}/N{N}_U{U}_s{s}_sigma{sigma}_batch.log"
    params:
        mode_flag=lambda wc: MODES[wc.mode]
        n_iter=100
    resources:
        mem_mb= 12 * 1000,
        runtime = 120,
    threads: 5
    shell:
    """
        escsim run \
        -t -w {params.n_iter} -j 5 -m {params.mode_flag} \
        --folder results/trees/{wildcards.mode} \
        {wildcards.N} {wildcards.s} {wildcards.U} {wildcards.sigma} \
    """

rule generate_coalescent_frequencies:
    input:
        expand(
            "results/trees/{{mode}}/N{{N}}_U{{U}}_s{{s}}_sigma{{sigma}}_{i}.trees",
            i=range(100)
        )
    output:
        results/coalescent_densities/N{N}_U{U}_s{s}_sigma{sigma}.txt
    logs:
    params:
        n_iter=config["constants"]["N_ITER"]
        n_sam=1000
    resources:
        mem_mb = 10 * 1000,
        runtime = 60,
    threads: 1
    shell:
    """
        python workflow/scripts/tree.py \
        -t N{wildcards.N}_U{wildcards.U}_s{wildcards.s}_sigma{wildcards.sigma} \
        -i results/trees/{wildcards.mode} {params.n_iter} {params.n_sam} > {log} 2>&1
    """
        

