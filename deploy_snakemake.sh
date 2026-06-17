#!/usr/bin/env bash

# ---------------------------------------------------------------------
# SLURM DIRECTIVES (For the master Snakemake process itself)
# ---------------------------------------------------------------------
#SBATCH --job-name=snakemake
#SBATCH --output=logs/snakemake_%j.out
#SBATCH --error=logs/snakemake_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=7-00:00:00
#SBATCH --partition="pcmpg_el8"

# Exit immediately if a command exits with a non-zero status
set -euo pipefail

# ---------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------
ENV_PATH="/home/yjauslin/.conda/envs/coalescent_simulations"
SNAKEFILE="Snakefile"
MAX_JOBS=50  # Maximum number of concurrent cluster jobs to submit
PROJECT_DIR=/data/projects/p2026-0005_simplify_evolution/Yannick/simplifying_evolution/

echo "====================================================================="
echo "           Starting Snakemake HPC Deployment Pipeline"
echo "====================================================================="
date

cd $PROJECT_DIR

# 1. Create logs directory if it doesn't exist
mkdir -p logs

# 2. Load Conda/Mamba module (Adapt this to your cluster's module names)
echo "Loading conda module..."
module load Anaconda3/2022.05 || echo "Conda module not found, relying on shell defaults."

export PATH="$ENV_PATH/bin:$PATH"

# 4. Run Snakemake with the Slurm Executor
echo "Executing Snakemake pipeline via Slurm..."
echo "---------------------------------------------------------------------"

snakemake --profile profiles/slurm

echo "---------------------------------------------------------------------"
echo "Master Snakemake script finished."
date