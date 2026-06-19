import numpy as np
import pandas as pd
import yaml

# Helper to format values uniformly to 15 decimal places for Snakemake wildcards
def format_fs_flat(value):
    v = float(value)

    if v == 0:
        return "0.0"

    # keep 15 decimals but avoid scientific drift, then strip zeros
    s = f"{v:.15f}".rstrip('0').rstrip('.')

    # ensure at least one decimal for consistency
    return s if '.' in s else s + ".0"


def get_clean_dataframes_by_experiment(config):
    sigma_range = config["constants"]["SIGMA_RANGE"]
    n_sigma = config["constants"]["N_SIGMA"]
    n_sel_coef = config["constants"]["N_SEL_COEF"]

    # Dictionary to hold dataframes for each experiment
    experiment_dfs = {}

    for exp_name, exp in config["experiments"].items():
        pop_size = exp["pop_size"][0]
        fixed_rows = []
        normal_rows = []

        if exp_name == "s_eff":
            mut_rate = exp["mut_rate"][0]
            for s in exp["s_normal"]:
                s_values = np.linspace(exp["s_range"][0] * s, s, n_sel_coef)
                for sel_coef in s_values:
                    fixed_rows.append([pop_size, format_fs_flat(mut_rate), sel_coef])

                sigma_values = np.round(
                    np.linspace(sigma_range[0] * s, sigma_range[1] * s, n_sigma), 15
                )
                for sigma in sigma_values:
                    normal_rows.append([pop_size, format_fs_flat(mut_rate), s, format_fs_flat(sigma)])

        elif exp_name == "U_eff":
            sel_coef = exp["sel_coef"][0]
            for U in exp["mut_rate_normal"]:
                U_values = np.linspace(
                    exp["mut_rate_range"][0] * U,
                    (1 + exp["mut_rate_range"][0]) * U,
                    n_sel_coef,
                )
                for mut_rate in U_values:
                    fixed_rows.append([pop_size, format_fs_flat(mut_rate), sel_coef])

                sigma_values = np.round(
                    np.linspace(sigma_range[0] * sel_coef, sigma_range[1] * sel_coef, n_sigma), 15
                )
                for sigma in sigma_values:
                    normal_rows.append([pop_size, format_fs_flat(U), sel_coef, format_fs_flat(sigma)])

        # Create, deduplicate, and store dataframes specifically for THIS experiment
        df_fixed = pd.DataFrame(fixed_rows, columns=["N", "U", "s"]).drop_duplicates().reset_index(drop=True)
        df_normal = pd.DataFrame(normal_rows, columns=["N", "U", "s", "sigma"]).drop_duplicates().reset_index(drop=True)

        experiment_dfs[exp_name] = {
            "fixed": df_fixed,
            "normal": df_normal
        }

    return experiment_dfs

# Load data
with open("config/config.yml", "r") as file:
    config = yaml.safe_load(file)

# exp_data["s_eff"]["fixed"] will give you the specific dataframe
exp_data = get_clean_dataframes_by_experiment(config)

# Create master lists combining all experiments for the summary rules
all_fixed_N = exp_data["s_eff"]["fixed"]["N"].tolist() + exp_data["U_eff"]["fixed"]["N"].tolist()
all_fixed_U = exp_data["s_eff"]["fixed"]["U"].tolist() + exp_data["U_eff"]["fixed"]["U"].tolist()
all_fixed_s = exp_data["s_eff"]["fixed"]["s"].tolist() + exp_data["U_eff"]["fixed"]["s"].tolist()

all_normal_N = exp_data["s_eff"]["normal"]["N"].tolist() + exp_data["U_eff"]["normal"]["N"].tolist()
all_normal_U = exp_data["s_eff"]["normal"]["U"].tolist() + exp_data["U_eff"]["normal"]["U"].tolist()
all_normal_s = exp_data["s_eff"]["normal"]["s"].tolist() + exp_data["U_eff"]["normal"]["s"].tolist()
all_normal_sigma = exp_data["s_eff"]["normal"]["sigma"].tolist() + exp_data["U_eff"]["normal"]["sigma"].tolist()

# Create combined lists for the effective selection coefficient aggregation rule
seff_combined_N = exp_data["s_eff"]["fixed"]["N"].tolist() + exp_data["U_eff"]["fixed"]["N"].tolist()
seff_combined_U = exp_data["s_eff"]["fixed"]["U"].tolist() + exp_data["U_eff"]["fixed"]["U"].tolist()
seff_combined_s = exp_data["s_eff"]["fixed"]["s"].tolist() + exp_data["U_eff"]["fixed"]["s"].tolist()