MODES = {
    "fixed": "f",
    "normal": "n"
}

escsim_outputs = [
    f"results/escsim/{mode}/escsim_N{p['N']}_U{p['U']}_s{p['s']}_sigma{p['sigma']}.out"
    for mode in MODES
    for p in PARAM_COMBINATIONS
]

rule escsim_run:
    output:
        sim="results/escsim/{mode}/escsim_N{N}_U{U}_s{s}_sigma{sigma}.out",
        wave="results/escsim/{mode}/wave_N{N}_U{U}_s{s}_sigma{sigma}.out"
    log:
        "logs/escsim/{mode}/escsim_N{N}_U{U}_s{s}_sigma{sigma}.log"
    params:
        CHRMLEN=config["constants"]["CHRMLEN"],
        BURNIN=config["constants"]["BURNIN"],
        mode_flag=lambda wc: MODES[wc.mode]
    shell:
        """
        escsim run \
            -f results/escsim/{wildcards.mode} \
            -m {params.mode_flag} \
            {wildcards.N} {wildcards.s} {wildcards.U} {wildcards.sigma} \
            {params.CHRMLEN} {params.BURNIN} 2>&1 | tee {log}
        """

rule escsim_summarize:
    input:
        escsim_outputs
    output:
        "results/escsim_figures/{mode}/FIGURE_PDF.pdf",
        "results/escsim_figures/{mode}/FIGURE_PDF_mean_velocity_summary.pdf"
    log:
        "logs/escsim_summarize_{mode}.log"
    shell:
        """
        escsim summarize \
            -i results/escsim/{wildcards.mode} \
            -o results/escsim_figures/{wildcards.mode} \
            FIGURE_PDF 2>&1 | tee {log}
        """