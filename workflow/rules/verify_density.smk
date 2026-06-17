frequency_outputs_fixed = [
    f"results/coalescent_densities/fixed/N{p['N']}_U{p['U']}_s{p['s']}.txt"
    for p in FIXED_PARAM_COMBINATIONS
]

frequency_outputs_normal = [
    f"results/coalescent_densities/normal/N{p['N']}_U{p['U']}_s{p['s']}_sd{p['sigma']:.1e}.txt"
    for p in NORMAL_PARAM_COMBINATIONS
]


rule generate_trees_fixed:
    input:
        "workflow/scripts/tree.slim"
    output:
        "results/trees/fixed_N{N}_U{U}_s{s}/N{N}_U{U}_s{s}_{sim_id}.trees"
    log:
        "logs/trees/fixed/N{N}_U{U}_s{s}/{sim_id}.log"
    params:
        CHRMLEN=config["constants"]["CHRMLEN"],
        BURNIN=config["constants"]["BURNIN"],
        folder='results/trees/fixed_N{N}_U{U}_s{s}/'
    resources:
        mem_mb= 1000,
        runtime = 60,
    threads: 1
    shell:
        """
        slim -d popsize={wildcards.N} \
             -d selcoef={wildcards.s} \
             -d mutrate={wildcards.U} \
             -d seqlen={params.CHRMLEN} \
             -d burnin={params.BURNIN} \
             -d OUTPUT_FOLDER='"{params.folder}"' \
             -d SIM_ID={wildcards.sim_id} \
             {input}
        """

rule generate_trees_normal:
    input:
        "workflow/scripts/tree_normal.slim"
    output:
        "results/trees/normal_N{N}_U{U}_s{s}_sd{sigma}/N{N}_U{U}_s{s}_sd{sigma}_{sim_id}.trees"
    log:
        "logs/trees/normal/N{N}_U{U}_s{s}_sd{sigma}/{sim_id}.log"
    params:
        CHRMLEN=config["constants"]["CHRMLEN"],
        BURNIN=config["constants"]["BURNIN"],
        folder="results/trees/normal_N{N}_U{U}_s{s}_sd{sigma}/"
    resources:
        mem_mb= 1000,
        runtime = 60,
    threads: 1
    shell:
        """
        slim -d popsize={wildcards.N} \
             -d selcoef={wildcards.s} \
             -d sigma={wildcards.sigma} \
             -d mutrate={wildcards.U} \
             -d seqlen={params.CHRMLEN} \
             -d burnin={params.BURNIN} \
             -d OUTPUT_FOLDER='"{params.folder}"' \
             -d SIM_ID={wildcards.sim_id} \
             {input}
        """

rule generate_coalescent_frequencies_fixed:
    wildcard_constraints:
        N = r"\d+",
        U = r"[\deE.+-]+",
        s = r"[\deE.+-]+"
    input:
        expand(
            "results/trees/fixed_N{{N}}_U{{U}}_s{{s}}/N{{N}}_U{{U}}_s{{s}}_{i}.trees",
            i=range(config["constants"]["N_SIM"])
        )
    output:
        "results/coalescent_densities/fixed/N{N}_U{U}_s{s}.txt"
    log:
        "logs/frequencies/fixed/N{N}_U{U}_s{s}.log"
    params:
        n_sim=config["constants"]["N_SIM"],
        n_sam=config["constants"]["N_SAM"]
    resources:
        mem_mb = 10 * 1000,
        runtime = 60,
    threads: 1
    shell:
        """
        python workflow/scripts/tree.py \
        -t N{wildcards.N}_U{wildcards.U}_s{wildcards.s} \
        -i results/trees/fixed_N{wildcards.N}_U{wildcards.U}_s{wildcards.s}/ \
        -o results/coalescent_densities/fixed \
        {params.n_sim} {params.n_sam} 2>&1 | tee {log}
        """

rule generate_coalescent_frequencies_normal:
    wildcard_constraints:
        N = r"\d+",
        U = r"[\deE.+-]+",
        s = r"[\deE.+-]+"
    input:
        expand(
            "results/trees/normal_N{{N}}_U{{U}}_s{{s}}_sd{{sigma}}/N{{N}}_U{{U}}_s{{s}}_sd{{sigma}}_{i}.trees",
            i=range(config["constants"]["N_SIM"])
        )
    output:
        "results/coalescent_densities/normal/N{N}_U{U}_s{s}_sd{sigma}.txt"
    log:
        "logs/frequencies/normal/N{N}_U{U}_s{s}_sigma{sigma}.log"
    params:
        n_iter=config["constants"]["N_SIM"],
        n_sam=config["constants"]["N_SAM"]
    resources:
        mem_mb = 10 * 1000,
        runtime = 60,
    threads: 1
    shell:
        """
        python workflow/scripts/tree.py \
        -t N{wildcards.N}_U{wildcards.U}_s{wildcards.s}_sd{wildcards.sigma} \
        -i results/trees/normal_N{wildcards.N}_U{wildcards.U}_s{wildcards.s}_sd{wildcards.sigma}/ \
        -o results/coalescent_densities/normal \
        {params.n_iter} {params.n_sam} 2>&1 | tee {log}
        """
        

