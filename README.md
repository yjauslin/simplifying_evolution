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
  - [escsim_figures/](results/escsim_figures/): Folder containing visualizations of the simulation results.
      - [comparison_plots](results/escsim_figures/comparison_plots/): Folder containing visualizations of the Kolmogorov-Smirnov and RMSE along the fixed grid of s_fixed (s_eff) and U_fixed (U_eff)
      - [fixed](results/escsim_figures/fixed/): Folder containing summary pdf for the fixed simulations and wave visualizations 
      - [normal](results/escsim_figures/normal/): Folder containing summary pdf for the normally distributed simulations and wave visualizations 

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


