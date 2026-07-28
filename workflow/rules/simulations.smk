import math
import re

# Fallback definition for custom math helper logic if 's' wildcard is omitted
S_FIXED = config["experiments"]["U_eff"]["sel_coef"][0]

def get_mem_mb(wildcards):
    # wildcards.N is a string, convert to int
    n_val = int(wildcards.N)
    
    if n_val >= 10000:
        return 210 * 1000  # 210 GB for large populations
    else:
        return 100 * 1000  # 100 GB default


def get_num_sim(wildcards):
    """
    Determine number of simulations from N, U, s, and optionally sigma.

    Uses:
    - wildcards.s if available, otherwise global S_FIXED
    - wildcards.sigma if available
    """

    float_pattern = re.compile(r"^-?\d*\.?\d+(?:[eE][+-]?\d+)?$")

    def parse_float(name, value):
        if value is None:
            return None
        if not float_pattern.match(str(value)):
            raise ValueError(f"Invalid wildcard '{name}' = '{value}'")
        return float(value)

    # required
    N = parse_float("N", getattr(wildcards, "N"))
    U = parse_float("U", getattr(wildcards, "U"))

    # optional: s may be missing
    s_raw = getattr(wildcards, "s", None)
    s = parse_float("s", s_raw) if s_raw is not None else S_FIXED

    # optional: sigma may be missing
    sigma_raw = getattr(wildcards, "sigma", None)
    sigma = parse_float("sigma", sigma_raw)

    # safety check for s = 0
    if s == 0:
        return config["constants"]["ESTIMATE_N_SIM"]

    # If sigma is present, check its ratio to s
    if sigma is not None:
        if sigma / s > 0.25:
            return config["constants"]["ESTIMATE_N_SIM"]

    # original calculation
    phi = N * s * math.exp(-U / s)

    n_sim = config["constants"]["ESTIMATE_N_SIM"] if phi < 1 or U <= 0.001 else 100

    return n_sim

rule escsim_run_fixed:
    wildcard_constraints:
        N = r"\d+",
        U = r"[\deE.+-]+",
        s = r"[\deE.+-]+"
    output:
        sim="results/escsim/fixed/escsim_N{N}_U{U}_s{s}.out",
        wave=temp("results/escsim/fixed/wave_N{N}_U{U}_s{s}.out")
    log:
        "logs/escsim/fixed/escsim_N{N}_U{U}_s{s}.log"
    params:
        CHRMLEN=config["constants"]["CHRMLEN"],
        BURNIN=config["constants"]["BURNIN"],
        w_val=get_num_sim,
    threads: 20
    shadow: "minimal"
    resources:
        mem_mb=get_mem_mb,        
        runtime=180
    shell:
        """
        escsim run \
            -f results/escsim/fixed \
            -m f \
            -w {params.w_val} \
            -j {threads} \
            {wildcards.N} {wildcards.s} {wildcards.U} \
            {params.CHRMLEN} {params.BURNIN} 2>&1 | tee {log}
        """

rule escsim_summarize_fixed:
    input:
        "results/markers/escsim_run_s_eff_fixed.done",
        "results/markers/escsim_run_U_eff_fixed.done"
    output:
        "results/escsim_figures/fixed/FIGURE_PDF.pdf",
        "results/escsim_figures/fixed/FIGURE_PDF_mean_velocity_summary.pdf"
    log:
        "logs/escsim_summarize_fixed.log"
    resources:
        mem_mb=50*1000,        
        runtime=4320,
    threads: 20
    shell:
        """
        escsim summarize \
            -i results/escsim/fixed \
            -o results/escsim_figures/fixed \
            -m f FIGURE_PDF 2>&1 | tee {log}
        """

rule escsim_run_normal:
    wildcard_constraints:
        N = r"\d+",
        U = r"[\deE.+-]+",
        s = r"[\deE.+-]+",
        sigma = r"[\deE.+-]+"
    params:
        CHRMLEN=config["constants"]["CHRMLEN"],
        BURNIN=config["constants"]["BURNIN"],
        w_val=get_num_sim,
    output:
        sim=f"results/escsim/normal/escsim_N{{N}}_U{{U}}_s{{s}}_sd{{sigma}}.out",
        wave=temp("results/escsim/normal/wave_N{N}_U{U}_s{s}_sd{sigma}.out")
    log:
        f"logs/escsim/normal/escsim_N{{N}}_U{{U}}_s{{s}}_sd{{sigma}}.log"
    threads: 20
    shadow: "minimal"
    resources:
        mem_mb=get_mem_mb,        
        runtime=180
    shell:
        """
        escsim run \
            -f results/escsim/normal \
            -m n \
            -w {params.w_val} \
            -j {threads} \
            {wildcards.N} {wildcards.s} {wildcards.U} {wildcards.sigma} \
            {params.CHRMLEN} {params.BURNIN} 2>&1 | tee {log}
        """

rule escsim_summarize_normal:
    input:
        "results/markers/escsim_run_s_eff_normal.done",
        "results/markers/escsim_run_U_eff_normal.done"
    output:
        "results/escsim_figures/normal/FIGURE_PDF.pdf",
        "results/escsim_figures/normal/FIGURE_PDF_mean_velocity_summary.pdf"
    log:
        "logs/escsim_summarize_normal.log"
    resources:
        mem_mb=50*1000,        
        runtime=4320,
    threads: 20
    shell:
        """
        escsim summarize \
            -i results/escsim/normal \
            -o results/escsim_figures/normal \
            -m n FIGURE_PDF 2>&1 | tee {log}
        """

# ==============================================================================
# BATCH MARKER AGGREGATIONS (Prevents Master DAG Memory Bloat)
# ==============================================================================

rule gather_escsim_run_s_eff_fixed:
    input:
        expand(
            "results/escsim/fixed/escsim_N{N}_U{U}_s{s}.out",
            zip,
            N=exp_data["s_eff"]["fixed"]["N"],
            U=exp_data["s_eff"]["fixed"]["U"],
            s=exp_data["s_eff"]["fixed"]["s"],
        )
    output:
        "results/markers/escsim_run_s_eff_fixed.done"
    shell:
        "touch {output}"

rule gather_escsim_run_U_eff_fixed:
    input:
        expand(
            "results/escsim/fixed/escsim_N{N}_U{U}_s{s}.out",
            zip,
            N=exp_data["U_eff"]["fixed"]["N"],
            U=exp_data["U_eff"]["fixed"]["U"],
            s=exp_data["U_eff"]["fixed"]["s"],
        )
    output:
        "results/markers/escsim_run_U_eff_fixed.done"
    shell:
        "touch {output}"

rule gather_escsim_run_s_eff_normal:
    input:
        expand(
            "results/escsim/normal/escsim_N{N}_U{U}_s{s}_sd{sigma}.out",
            zip,
            N=exp_data["s_eff"]["normal"]["N"],
            U=exp_data["s_eff"]["normal"]["U"],
            s=exp_data["s_eff"]["normal"]["s"],
            sigma=exp_data["s_eff"]["normal"]["sigma"],
        )
    output:
        "results/markers/escsim_run_s_eff_normal.done"
    shell:
        "touch {output}"

rule gather_escsim_run_U_eff_normal:
    input:
        expand(
            "results/escsim/normal/escsim_N{N}_U{U}_s{s}_sd{sigma}.out",
            zip,
            N=exp_data["U_eff"]["normal"]["N"],
            U=exp_data["U_eff"]["normal"]["U"],
            s=exp_data["U_eff"]["normal"]["s"],
            sigma=exp_data["U_eff"]["normal"]["sigma"],
        )
    output:
        "results/markers/escsim_run_U_eff_normal.done"
    shell:
        "touch {output}"