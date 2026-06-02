import numpy as np

# Collect all parameter combinations
PARAM_COMBINATIONS = []

for exp_name, exp in config["experiments"].items():
    s = exp["sel_coef"][0]
    sigma_values = np.linspace(0.1 * s, 0.25 * s, 12)
    for sigma in sigma_values:
        # Extract n_iter from config (e.g., 100)
        num_simulations = config["constants"]["N_SIM"]

        
        for i in range(num_simulations):
            PARAM_COMBINATIONS.append({
                "N": exp["pop_size"][0],
                "U": exp["mut_rate"][0],
                "s": exp["sel_coef"][0],
                "sigma": sigma,
                "n_sim": i,  # This will be 0, 1, 2... 99
                "exp": exp_name
            })
    
   