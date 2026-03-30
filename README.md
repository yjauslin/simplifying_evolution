# Simplifying evolution

## Description

This project investigates the distribution of fitness effects in non-recombining regions. It does so by implementing the extended structured coalescent by Strütt et al. (2025) and extending it to include different distributions of fitness effects for the selection coefficient.
## Structure

- config
  - [environment.yml](config/environment.yml): Config-file listing all the necessary dependencies of this project

- resources/
  - [escsim-0.1.1/](resources/escsim-0.1.1/): Pip-package implementing the forward simulation approach to simulate mutational burden profile.

- results/
  - [escsim_figures/](results/escsim_figures/): Folder containing simulation results.

- workflow/
  - [rules/](workflow/rules/)
    - [calc_density.smk](workflow/rules/calc_density.smk): Defines the rule to calculate the coalescent density.
    - [plots.smk](workflow/rules/plots.smk): Defines the rules for creating wave and result plots
    - [simulations.smk]((workflow/rules/simulations.smk)): Defines the rules for forward simulations including summary figures.
  - [scripts/](workflow/scripts/)
    - [calc_coalescent_densities.py](workflow/scripts/calc_coalescent_densities.py): Calculates the coalescent densities and writes the waves into result files based upon the simulation results from the forward-time simulations.
    - [visualize_coalescent_densities.py](workflow/scripts/visualize_coalescent_densities.py): Plots the calculated densities plus the effective population size.
    - [visualize_wave.py](workflow/scripts/visualize_wave.py): Plots individual waves to verify whether they are equilibrated
  - [snakefile](workflow/snakefile): File listing all the necessary rules to produce results.


