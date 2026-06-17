import math
import re

S_FIXED = config["experiments"]["RMSE"]["s_fixed"][0]

# Helper to format sigma up to 15 decimals without trailing zero padding
def clean_sigma(val):
    # Formats to 15 decimals, then strips unnecessary right-side zeros
    s = f"{float(val):.15f}".rstrip('0').rstrip('.')
    # Handle edge case where it truncates to exactly '0' if the number was too small
    return s if s != "" else "0"

def get_mem_mb(wildcards):
    # wildcards.N is a string, convert to int
    n_val = int(wildcards.N)
    
    if n_val >= 10000:
        return 210 * 1000  # 210 GB for large populations
    else:
        return 100 * 1000  # 100 GB default


def get_num_sim(wildcards):
    """
    Determine number of simulations from N, U and s.

    Uses:
    - wildcards.s if available
    - otherwise global S_FIXED
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

    # safety check
    if s == 0:
        return 200

    phi = N * s * math.exp(-U / s)

    return 100 if phi > 1 else 200

escsim_outputs_fixed = [
    f"results/escsim/fixed/escsim_N{p['N']}_U{p['U']}_s{p['s']}.out"
    for p in FIXED_PARAM_COMBINATIONS
]

escsim_outputs_normal = [
    f"results/escsim/normal/escsim_N{p['N']}_U{p['U']}_s{p['s']}_sd{clean_sigma(p['sigma'])}.out"
    for p in NORMAL_PARAM_COMBINATIONS
]

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
    threads: 16
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
            -j 16 \
            {wildcards.N} {wildcards.s} {wildcards.U} \
            {params.CHRMLEN} {params.BURNIN} 2>&1 | tee {log}
        """

rule escsim_summarize_fixed:
    input:
        escsim_outputs_fixed
    output:
        "results/escsim_figures/fixed/FIGURE_PDF.pdf",
        "results/escsim_figures/fixed/FIGURE_PDF_mean_velocity_summary.pdf"
    log:
        "logs/escsim_summarize_fixed.log"
    resources:
        mem_mb=5*1000,        
        runtime=60,
    threads: 4
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
        sigma = r"[\deE.+-]+"
    params:
        CHRMLEN=config["constants"]["CHRMLEN"],
        BURNIN=config["constants"]["BURNIN"],
        w_val=get_num_sim,
        # Truncates to max 15 decimals, strips trailing zeros to preserve short values
        sigma_formatted=lambda wc: f"{float(wc.sigma):.15f}".rstrip('0').rstrip('.')
    output:
        sim=f"results/escsim/normal/escsim_N{{N}}_U{{U}}_s{S_FIXED}_sd{{sigma}}.out",
        wave=temp(f"results/escsim/normal/wave_N{{N}}_U{{U}}_s{S_FIXED}_sd{{sigma}}.out")
    log:
        f"logs/escsim/normal/escsim_N{{N}}_U{{U}}_s{S_FIXED}_sd{{sigma}}.log"
    threads: 16
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
            -j 16 \
            {wildcards.N} {S_FIXED} {wildcards.U} {params.sigma_formatted} \
            {params.CHRMLEN} {params.BURNIN} 2>&1 | tee {log}
        """

rule escsim_summarize_normal:
    input:
        escsim_outputs_normal
    output:
        "results/escsim_figures/normal/FIGURE_PDF.pdf",
        "results/escsim_figures/normal/FIGURE_PDF_mean_velocity_summary.pdf"
    log:
        "logs/escsim_summarize_normal.log"
    resources:
        mem_mb=5*1000,        
        runtime=60,
    threads: 4
    shell:
        """
        escsim summarize \
            -i results/escsim/normal \
            -o results/escsim_figures/normal \
            -m n FIGURE_PDF 2>&1 | tee {log}
        """