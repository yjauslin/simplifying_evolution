MODES = {
    "fixed": "f",
    "normal": "n"
}

escsim_outputs = [
    f"results/escsim/escsim_{mode}_N{p['N']}_U{p['U']}_s{p['s']}_sigma{p['sigma']}.out"
    for mode in MODES
    for p in PARAM_COMBINATIONS
]

rule escsim_run:
    output:
        sim="results/escsim/escsim_{mode}_N{N}_U{U}_s{s}_sigma{sigma}.out",
        wave="results/escsim/wave_{mode}_N{N}_U{U}_s{s}_sigma{sigma}.out"
    log:
        "logs/escsim/escsim_{mode}_N{N}_U{U}_s{s}_sigma{sigma}.log"
    params:
        CHRMLEN=config["constants"]["CHRMLEN"],
        BURNIN=config["constants"]["BURNIN"],
        mode_flag=lambda wc: MODES[wc.mode]
    shell:
        """
        escsim run \
            -f results/escsim \
            -m {params.mode_flag} \
            {wildcards.N} {wildcards.s} {wildcards.U} {wildcards.sigma} \
            {params.CHRMLEN} {params.BURNIN} 2>&1 | tee {log}
        """

rule escsim_summarize:
    input:
        escsim_outputs
    output:
        "results/escsim_figures/FIGURE_PDF.pdf",
        "results/escsim_figures/FIGURE_PDF_mean_velocity_summary.pdf"
    log:
        "logs/escsim_summarize.log"
    shell:
        """
        escsim summarize \
            -i results/escsim \
            -o results/escsim_figures \
            FIGURE_PDF 2>&1 | tee {log}
        """