# Collect all parameter combinations
PARAM_COMBINATIONS = []

for exp_name, exp in config["experiments"].items():
    
    # Sweep pop_size
    if exp_name == "low_phi":
        for sigma in exp["sigma"]:
            PARAM_COMBINATIONS.append({
                "N": exp["pop_size"][0],
                "U": exp["mut_rate"][0],
                "s": exp["sel_coef"][0],
                "sigma": sigma,
                "exp": exp_name
            })
    
    # Sweep sel_coef
    elif exp_name == "intermediate_phi":
        for sigma in exp["sigma"]:
            PARAM_COMBINATIONS.append({
                "N": exp["pop_size"][0],
                "U": exp["mut_rate"][0],
                "s": exp["sel_coef"][0],
                "sigma": sigma,
                "exp": exp_name
            })

    # Sweep sigma
    elif exp_name == "high_phi":
        for sigma in exp["sigma"]:
            PARAM_COMBINATIONS.append({
                "N": exp["pop_size"][0],
                "U": exp["mut_rate"][0],
                "s": exp["sel_coef"][0],
                "sigma": sigma,
                "exp": exp_name
            })