# Simplifying evolution

## Description

This project investigates the distribution of fitness effects in non-recombining regions. It does so by implementing the extended structured coalescent by Strütt et al. (2025) and adding an efficient selection coefficient ($s_e$) to it. 

## Structure

- config
  - [environment.yml](config/environment.yml): Config-file listing all the necessary dependencies of this project

- resources/
  - [escsim-0.1.1/](resources/escsim-0.1.1/): Pip-package implementing the forward simulation approach to simulate mutational burden profile.

- results/
  - [escsim/](results/escsim/): Folder containing simulation results of the esc model in forward time.

- workflow/
  - scripts/
    - calc_coalescent_densities.py: Calculates the coalescent densities and writes the waves into result files based upon the simulation results from the forward-time simulations.
    - visualize_coalescent_densities.py: Plots the calculated densities plus the effective population size.
    - visualize_wave.py: Plots individual waves to verify whether they are equilibrated
  - [snakefile](workflow/snakefile): File listing all the necessary rules to produce results.


