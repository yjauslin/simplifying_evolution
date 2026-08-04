# Simplifying evolution

## Description

The effect of many weakly deleterious mutations can be summarized by a few strongly deleterious variants (Good, 2012). This effect cannot be simply explained by a rescaled population size and additionally an effective selection coefficient and mutation rate are needed. Above concept, applies under weak selection in non-recombining regions, but was never tested under other selection strength regimes such as intermediate or strong selection. This project addresses this gap using the extended structured coalescent (ESC) by Strütt et al. (2025) and adding a normally distributed selection coefficient into its framework. We derive the effective selection coefficient and mutation rate under different degrees of the standard deviation and for three distinct selection regimes (weak, intermediate and strong). In addition, we also examine under what distribution of fitness effects the ESC's assumption of a fixed selection coefficient is justified.

## Structure

- config
  - [environment.yml](config/environment.yml): Config-file listing all the necessary dependencies of this project
  - [config.yml](config/config.yml): Config-file setting up the different parameters used in the simulations.

- resources/
  - [escsim-0.1.1/](resources/escsim-0.1.1/): Pip-package implementing the forward simulation approach to simulate mutational burden profile.

- results/
  - [coalescent_densities](results/coalescent_densities/): Folder containing coalescent density estimates and empirical coalescent density data to produce Figure 12 and 13.
    - [fixed](results/coalescent_densities/fixed/): Folder containing coalescent density estimates and empirical coalescent density data with the fixed selection coefficient.
    - [normal](results/coalescent_densities/normal/): Folder containing coalescent density estimates and empirical coalescent density data with the Gaussian distribution of mutational effects.
  - [escsim](results/escsim/): Folder containing relative click rate and mutational burden profile data obtained through SLiM simulations used for figure 6 and 7.
    - [fixed](results/escsim/fixed/): Folder containing relative click rate and mutational burden profile data obtained through SLiM simulations with a fixed selection coefficient. This data is a precursor to obtain the coalescent densities for Figure 12 and 13.
    - [normal](results/escsim/normal/): Folder containing relative click rate and mutational burden profile data obtained through SLiM simulations with the Gaussian distribution of mutational effects. This data was used to produce Figure 6 and 7.
  - [escsim_figures/](results/escsim_figures/): Folder containing visualizations of the simulation results.
      - [fixed](results/escsim_figures/fixed/): Folder containing summary pdf for the fixed simulations. 
      - [normal](results/escsim_figures/normal/): Folder containing summary pdf for the normally distributed simulations.
      - [coalescent_density_s_eff.jpg](results/escsim_figures/coalescent_density_s_eff.jpg): Figure 12, Coalescent densities over time for the fixed effective selection coefficient estimates, normally distributed selection coefficient estimates and Wright-Fisher forward simulations. 
      - [coalescent_density_U_eff.jpg](results/escsim_figures/coalescent_density_u_eff.jpg): Figure 13, Coalescent densities over time for the effective mutation rate estimates with a fixed selection coefficient, normally distributed selection coefficient estimates and Wright-Fisher forward simulations. 
      - [combined_effective_mutation_rate.jpg](results/escsim_figures/combined_effective_mutation_rate.jpg): Figure 10, Effective mutation rate and relative effective mutation rate under increasing standard variation.
      - [combined_effective_selection_coefficient.jpg](results/escsim_figures/combined_effective_selection_coefficient.jpg): Figure 8, Effective selection coefficient and relative effective selection coefficient under increasing standard variation.
      - [combined_velocity.jpg](results/escsim_figures/combined_velocity.jpg): Figure 6, Relative click rate under increasing standard variation for different selection coefficients and different mutation rates.
      - [effective_mutation_rate.jpg](results/escsim_figures/effective_mutation_rate.jpg): Subplot 10a, Effective mutation rate under increasing standard deviation.
      - [effective_selection_coefficient](results/escsim_figures/effective_selection_coefficient.jpg): Subplot 8a, Effective selection coefficient under increasing standard deviation.
      - [mut_burden_dist](results/escsim_figures/mut_burden_dist.jpg): Figure 7, Relative mutational burden classes for different selection coefficients and mutation rates under increasing standard variation. 
      - [relative_effective_mutation_rate.jpg](results/escsim_figures/relative_effective_mutation_rate.jpg): Subplot 10b, Relative effective mutation rate under increasing standard deviation.
      - [relative_effective_selection_coefficient.jpg](results/escsim_figures/relative_effective_selection_coefficient.jpg): Subplot 8b, Relative effective selection coefficient under increasing standard deviation.
      - [s_eff_Kolmogorov.jpg](results/escsim_figures/s_eff_Kolmogorov.jpg): Figure 9, Kolmogorov-Smirnov D-Statistic under increasing standard variation for three different selection coefficients.
      - [s_velocity.jpg](results/escsim_figures/s_velocity.jpg): Subplot 6a, Relative click rate under increasing standard deviation for different selection coefficients.
      - [U_eff_Kolmogorov.jpg](results/escsim_figures/U_eff_Kolmogorov.jpg): Figure 11, Kolmogorov-Smirnov D-Statistic under increasing standard variation for three different mutation rates.
      - [u_velocity.jpg](results/escsim_figures/u_velocity.jpg): Subplot 6b, Relative click rate under increasing standard deviation for different mutation rates.
  - [min_values](results/min_values/): Folder containing data on the selection coefficient or mutation rate minimizing the Kolmogorov-Smirnov D-statistic per value of standard deviation.
    - [s_eff](results/min_values/s_eff/): Folder containing data on the selection coefficient minimizing the Kolmogorov-Smirnov D-statistic per value of standard deviation.
    - [U_eff](results/min_values/U_eff/): Folder containing data on the mutation rate minimizing the Kolmogorov-Smirnov D-statistic per value of standard deviation.


- workflow/
  - [rules/](workflow/rules/)
    - [calc_density.smk](workflow/rules/calc_density.smk): Defines the rule to calculate the coalescent density.
    - [comparison.smk](workflow/rules/comparison.smk): Defines the rules that calculate the Kolmogorov-Smirnov D-Statistic and RMSE along the two grids. Defines the rule that picks the minimum value for the D-Statistic for every sigma along the grid.
    - [plots.smk](workflow/rules/plots.smk): Defines the rules for creating wave and result plots
    - [simulations.smk](workflow/rules/simulations.smk): Defines the rules for forward simulations including summary figures.
    - [verify_density.smk](workflow/rules/verify_density.smk): Defines the rules that verify the estimates using Wright-Fisher forward simulations.
  - [scripts/](workflow/scripts/)
    - [calc_coalescent_densities.py](workflow/scripts/calc_coalescent_densitiy.py): Calculates the coalescent densities and writes them into result files based upon the simulation results from the forward-time simulations.
    - [calculate_tmrca.py](workflow/scripts/calculate_tmrca.py): Script that samples 500 lineage pairs from Wright-Fisher forward simulations and calculates their time to most recent common ancestor (TMRCA).
    - [compare_distributions.py](workflow/scripts/compare_distributions.py): converts the coalescent densities probability distributions for fixed and normal coalescent density files into cumulative distribution functions and then calculates the root mean squared error (RMSE) and the Kolmogorov-Smirnov D-statistic.
    - [get_min_kolmogorov.py](workflow/scripts/get_min_kolmogorov.py): picks the minimum value of the D-statistic for every value along the sigma grid.
    - [tree_normal.slim](workflow/scripts/tree_normal.slim): SLiM simulation file that produces tree-files with a normally distributed selection coefficient to verify the coalescent density estimates.
    - [visualize_densities.py](workflow/scripts/visualize_densities.py): Plots the coalescent density estimates and the simulated densities for different values of the standard deviation.
    - [visualize_mut_burden_dist.py](workflow/scripts/visualize_mut_burden_dist.py): Plots the relative mutational burden distribution for different values of the standard deviation.
    - [visualize_rmse.py](workflow/scripts/visualize_rmse.py): Plots the RMSE and D-statistic along the grid of s_fixed or U_fixed.
    - [visualize_sd_vs_ks.py](workflow/scripts/visualize_sd_vs_ks.py): Plots the minimum value of the D-statistic versus the standard deviation.
    - [visualize_sd_vs_s.py](workflow/scripts/visualize_sd_vs_s.py): Plots the effective selection coefficient or mutation rate versus the standard deviation.
    - [visualize_velocity_vs_sd.py](workflow/scripts/visualize_velocity_vs_sd.py): Plots the mean velocity versus the standard deviation.
    - [visualize_wave.py](workflow/scripts/visualize_wave.py): Plots individual waves to check whether they are equilibrated
  - [params.smk](workflow/params.smk): extracts the parameter values from the config and produces datafranes that can be used to produce input or output for individual rules.
  - [snakefile](workflow/snakefile): File defining the necessary output of the pipeline
- [deploy_snakemake.sh](deploy_snakemake.sh): Bash-script that deploys the pipeline on a HPC
- [model_figure.ipynb](model_figure.ipynb): Jupyter-notebook that produces the model figures plus additional figures used in a lab-meeting presentation.


