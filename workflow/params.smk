# Collect all parameter combinations
PARAM_COMBINATIONS = []

for exp_name, exp in config["experiments"].items():
    
    # Sweep pop_size
    if exp_name == "vary_n":
        for N in exp["pop_size"]:
            PARAM_COMBINATIONS.append({
                "N": N,
                "U": exp["mut_rate"][0],
                "s": exp["sel_coef"][0],
                "sigma": exp["sigma"][0],
                "exp": exp_name
            })
    
    # Sweep sel_coef
    elif exp_name == "vary_s":
        for s in exp["sel_coef"]:
            PARAM_COMBINATIONS.append({
                "N": exp["pop_size"][0],
                "U": exp["mut_rate"][0],
                "s": s,
                "sigma": exp["sigma"][0],
                "exp": exp_name
            })

    # Sweep sigma
    elif exp_name == "vary_sigma":
        for sigma in exp["sigma"]:
            PARAM_COMBINATIONS.append({
                "N": exp["pop_size"][0],
                "U": exp["mut_rate"][0],
                "s": exp["sel_coef"][0],
                "sigma": sigma,
                "exp": exp_name
            })