import numpy as np

# Collect all parameter combinations
FIXED_PARAM_COMBINATIONS = []
NORMAL_PARAM_COMBINATIONS = []


for exp_name, exp in config["experiments"].items():
    base_s = exp["sel_coef"][0]
    s_fixed = exp["s_fixed"][0]

    num_simulations = config["constants"]["N_SIM"]

    # --------------------------------------------------
    # Fixed model: vary s
    # --------------------------------------------------
    s_values = np.linspace(0.2*base_s, base_s, 100)

    for sel_coef in s_values:
        for i in range(num_simulations):
            FIXED_PARAM_COMBINATIONS.append({
                "N": exp["pop_size"][0],
                "U": exp["mut_rate"][0],
                "s": sel_coef,
                "n_sim": i,
                "exp": exp_name,
            })

    # --------------------------------------------------
    # Normal model: keep s fixed, vary sigma
    # --------------------------------------------------
    sigma_values = np.linspace(
        0.1 * s_fixed,
        0.25 * s_fixed,
        12
    )
    
    sigma_values = np.round(sigma_values, 15)

    for sigma in sigma_values:
        for i in range(num_simulations):
            NORMAL_PARAM_COMBINATIONS.append({
                "N": exp["pop_size"][0],
                "U": exp["mut_rate"][0],
                "s": s_fixed,
                "sigma": sigma,
                "n_sim": i,
                "exp": exp_name,
            })  