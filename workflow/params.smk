# Collect all parameter combinations
PARAM_COMBINATIONS = []

for exp_name, exp in config["experiments"].items():
    for sigma in exp["sigma"]:
        # Extract n_iter from config (e.g., 100)
        num_iterations = config["constants"]["N_ITER"]
        
        for i in range(num_iterations):
            PARAM_COMBINATIONS.append({
                "N": exp["pop_size"][0],
                "U": exp["mut_rate"][0],
                "s": exp["sel_coef"][0],
                "sigma": sigma,
                "n_iter": i,  # This will be 0, 1, 2... 99
                "exp": exp_name
            })
    
   