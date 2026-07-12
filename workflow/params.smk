import numpy as np
import pandas as pd
import yaml

# Helper to format values uniformly to 15 decimal places for Snakemake wildcards
def format_fs_flat(value):
    # Convert value to float
    v = float(value)

    # if value is 0 return "0.0"
    if v == 0:
        return "0.0"

    # keep 15 decimals but avoid scientific drift, then strip zeros
    s = f"{v:.15f}".rstrip('0').rstrip('.')

    # ensure at least one decimal for consistency
    return s if '.' in s else s + ".0"


def get_clean_dataframes_by_experiment(config):
    """
    Generates structured, deduplicated DataFrames of simulation parameters 
    for different experiments based on a config file.

    It creates two types of parameter sets per experiment:
    1. 'fixed': Models fixed values of selection coefficients (s).
    2. 'normal': Models normally distributed selection coefficient, varying by standard deviation (sigma).

    Parameters:
    -----------
    config : dict
        A nested configuration dictionary containing global constants and specific experiment parameters.

    Returns:
    --------
    dict
        A dictionary structured as {experiment_name: {"fixed": df_fixed, "normal": df_normal}}
    """
    # --- 1. Extract Global Constants ---
    # These control the range scale and the number of steps/samples to generate for linspace
    sigma_range = config["constants"]["SIGMA_RANGE"]
    n_sigma = config["constants"]["N_SIGMA"]
    n_sel_coef = config["constants"]["N_SEL_COEF"]

    # Dictionary to hold the resulting dataframes for each experiment
    experiment_dfs = {}

    # --- 2. Iterate Through Each Configured Experiment ---
    for exp_name, exp in config["experiments"].items():
        pop_size = exp["pop_size"][0]  # Extract base population size (N)
        fixed_rows = []                # Holds rows for the fixed effects dataframe
        normal_rows = []               # Holds rows for the normally distributed effects dataframe

        # --- Experiment 1: Effective Selection Coefficient ("s_eff") ---
        if exp_name == "s_eff":
            mut_rate = exp["mut_rate"][0] # Mutation rate stays constant for this experiment
            
            for s in exp["s_normal"]:
                # Generate an array of evenly spaced selection coefficients (s_values) across the specified range
                s_values = np.linspace(exp["s_range"][0] * s, s, n_sel_coef+1)
                for sel_coef in s_values:
                    # Format and append parameters for the fixed model: [Population Size, Mutation Rate, Selection Coefficient]
                    fixed_rows.append([pop_size, format_fs_flat(mut_rate), format_fs_flat(sel_coef)])

                # Generate an array of evenly spaced standard deviations (sigma) around 's', rounded to prevent float inaccuracies
                sigma_values = np.round(
                    np.linspace(sigma_range[0] * s, sigma_range[1] * s, n_sigma+1), 15
                )
                for sigma in sigma_values:
                    # Append parameters for the normal distribution model (adds the sigma column)
                    normal_rows.append([pop_size, format_fs_flat(mut_rate), format_fs_flat(s), format_fs_flat(sigma)])

        # --- Experiment 2: Effective Mutation Rate ("U_eff") ---
        elif exp_name == "U_eff":
            sel_coef = exp["sel_coef"][0] # Selection coefficient stays constant for this experiment
            
            for U in exp["mut_rate_normal"]:
                # Generate an array of evenly spaced mutation rates (U_values) symmetrically distributed around U
                U_values = np.linspace(
                    exp["mut_rate_range"][0] * U,
                    U,
                    n_sel_coef+1,
                )
                for mut_rate in U_values:
                    # Format and append parameters for the fixed model: [Population Size, Mutation Rate, Selection Coefficient]
                    fixed_rows.append([pop_size, format_fs_flat(mut_rate), format_fs_flat(sel_coef)])

                # Generate an array of evenly spaced standard deviations (sigma) around 'sel_coef'
                sigma_values = np.round(
                    np.linspace(sigma_range[0] * sel_coef, sigma_range[1] * sel_coef, n_sigma+1), 15
                )
                for sigma in sigma_values:
                    # Append parameters for the normal distribution model (adds the sigma column)
                    normal_rows.append([pop_size, format_fs_flat(U), format_fs_flat(sel_coef), format_fs_flat(sigma)])

        # --- 3. Create, Clean, and Store DataFrames per Experiment ---
        # Convert lists of rows into Pandas DataFrames with appropriate column headers, remove potential duplicates
        df_fixed = pd.DataFrame(fixed_rows, columns=["N", "U", "s"]).drop_duplicates().reset_index(drop=True)
        df_normal = pd.DataFrame(normal_rows, columns=["N", "U", "s", "sigma"]).drop_duplicates().reset_index(drop=True)

        # Map the cleaned dataframes back to the experiment name
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